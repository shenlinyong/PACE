#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PACE Predictor Module

Core prediction functions for enhancer-gene regulatory interactions.

Author: Linyong Shen @ Northwest A&F University
"""

import os
from pathlib import Path
import numpy as np
import pandas as pd
from typing import Optional, List, Dict, Tuple, Union

from tools import (
    calculate_distance, power_law_contact, safe_divide,
    create_enhancer_gene_pairs, assign_enhancer_class, logger
)
from hic import ContactEstimator
from pace_core import score_pairs, ScoreConfig


class PACEPredictor:
    """
    PACE predictor for enhancer-gene interactions.
    """
    
    def __init__(self,
                 max_distance: int = 5000000,
                 contact_method: str = 'powerlaw',
                 hic_file: Optional[str] = None,
                 hic_type: str = 'hic',
                 hic_resolution: int = 5000,
                 hic_gamma: float = 1.024238616787792,
                 hic_scale: float = 5.9594510043736655,
                 score_column: str = 'PACE.Score',
                 include_self_promoter: bool = True):
        """
        Initialize predictor.
        
        Args:
            max_distance: Maximum enhancer-gene distance
            contact_method: Contact estimation method (powerlaw, hic)
            hic_file: Path to Hi-C file
            hic_type: Type of Hi-C file
            hic_resolution: Hi-C resolution
            hic_gamma: Power-law exponent
            hic_scale: Power-law scale
            score_column: Name for score column
            include_self_promoter: Include self-promoter interactions
        """
        self.max_distance = max_distance
        self.score_column = score_column
        self.include_self_promoter = include_self_promoter
        
        # Initialize contact estimator
        self.contact_estimator = ContactEstimator(
            method=contact_method,
            hic_file=hic_file,
            hic_type=hic_type,
            resolution=hic_resolution,
            hic_gamma=hic_gamma,
            hic_scale=hic_scale
        )
    
    def create_pairs(self, enhancers, genes):
        enhancers = enhancers.copy()
        genes = genes.copy()
        enhancers['chr'] = enhancers['chr'].astype(str)
        genes['chr'] = genes['chr'].astype(str)
        if 'gene_id' not in genes or genes.gene_id.isna().any():
            raise ValueError('Stable gene_id is required')
        if 'TSS' not in genes:
            if 'tss' in genes:
                genes['TSS'] = genes.tss
            elif {'start', 'end', 'strand'}.issubset(genes):
                genes['TSS'] = np.where(genes.strand == '-', genes.end - 1, genes.start)
            else:
                raise ValueError('Provide BED0 TSS or start/end/strand')
        tss_keys = ['chr', 'gene_id', 'TSS']
        for column in ['tss_weight', 'tss_quality', 'catalogue_quality',
                       'unassigned_mass', 'strand']:
            if column in genes and (genes.groupby(tss_keys)[column].nunique(dropna=False) > 1).any():
                raise ValueError(f'Conflicting {column} for duplicate gene/TSS records')
        genes = genes.drop_duplicates(tss_keys)
        if 'tss_weight' not in genes:
            genes['tss_weight'] = 1 / genes.groupby(['chr', 'gene_id']).TSS.transform('size')
        records = []
        for chrom, ec in enhancers.groupby('chr'):
            gc = genes[genes.chr == chrom].sort_values('TSS')
            positions = gc.TSS.to_numpy()
            for _, e in ec.iterrows():
                # Match the analysis adapter exactly, including odd-width BED
                # intervals and candidates exactly at the open window boundary.
                mid = (e.start + e.end) / 2
                lo = np.searchsorted(positions, mid-self.max_distance, side='right')
                hi = np.searchsorted(positions, mid+self.max_distance, side='left')
                for _, g in gc.iloc[lo:hi].iterrows():
                    dist = abs(mid-g.TSS)
                    row = dict(e)
                    row.update(TargetGene=g.get('gene_name', g.get('name', g.gene_id)),
                               TargetGeneEnsemblID=g.gene_id, TargetGeneTSS=g.TSS,
                               TargetGeneStrand=g.get('strand', '+'), distance=dist,
                               enhancer_mid=mid, tss_weight=g.tss_weight)
                    for c in ['tss_quality', 'catalogue_quality', 'unassigned_mass']:
                        if c in g: row[c] = g[c]
                    row['class'] = assign_enhancer_class(dist)
                    if self.include_self_promoter or row['class'] != 'promoter':
                        records.append(row)
        return pd.DataFrame(records)

    def estimate_contacts(self, pairs):
        pairs = pairs.copy()
        pairs['contact_prior'] = power_law_contact(
            pairs.distance.to_numpy(), hic_gamma=self.contact_estimator.hic_gamma,
            hic_scale=self.contact_estimator.hic_scale)
        pairs['contact_source'] = 'distance_prior'
        pairs['contact_observed'] = np.nan
        processor = self.contact_estimator.hic_processor
        if processor is not None:
            observations = []
            for _, r in pairs.iterrows():
                observations.append(processor.get_contact(
                    str(r['chr']), r.enhancer_mid, r.TargetGeneTSS,
                    use_powerlaw_fallback=False))
            pairs['contact_observed'] = observations
            # A readable contact file does not establish matching tissue or
            # sample provenance. Matching must be declared in contact metadata.
            pairs['contact_source'] = 'unknown'
        elif self.contact_estimator.avg_contacts:
            pairs['contact_observed'] = self.contact_estimator.estimate_batch(
                pairs, pos1_col='enhancer_mid', pos2_col='TargetGeneTSS')
            pairs['contact_source'] = 'surrogate'
        # Raw processed contacts without an expected curve and QC remain
        # visible but unscaled. They cannot silently replace the distance prior.
        return pairs

    def calculate_abc_score(self, pairs, activity_column='activity', contact_column='contact'):
        """Legacy method name; delegates to the sole current gene-level kernel."""
        if activity_column != 'activity':
            raise ValueError('Expression-weighted activity is not part of PACE')
        return score_pairs(pairs)

    def predict(self, enhancers, genes, expression=None, expression_weight=False,
                weight_method='log', min_expression=1.0, contact_metadata=None):
        enhancers = enhancers.copy()
        if 'activity' not in enhancers:
            raise ValueError('Provide quantified activity; missing activity is not 1')
        pairs = self.create_pairs(enhancers, genes)
        if pairs.empty:
            return pd.DataFrame(columns=['PACE.Score', 'evidence_status'])
        pairs = self.estimate_contacts(pairs)
        if contact_metadata is not None:
            keys = ['chr', 'start', 'end', 'TargetGeneEnsemblID', 'TargetGeneTSS']
            allowed = ['contact_observed', 'contact_expected', 'contact_reliability', 'contact_source']
            md = contact_metadata[keys + [c for c in allowed if c in contact_metadata]].copy()
            md['chr'] = md['chr'].astype(str)
            if md.duplicated(keys).any(): raise ValueError('Duplicate contact metadata keys')
            pairs = pairs.merge(md, on=keys, how='left', suffixes=('', '_qc'), validate='one_to_one')
            for c in allowed:
                if c + '_qc' in pairs:
                    pairs[c] = pairs[c + '_qc'].combine_first(pairs[c])
                    pairs = pairs.drop(columns=[c + '_qc'])
        predictions = score_pairs(pairs)
        if expression_weight:
            logger.warning('PACE retains RNA as context; expression multiplication is disabled')
        if expression is not None:
            idcol = next((c for c in ['gene_id', 'TargetGeneEnsemblID'] if c in expression), None)
            if idcol is None or 'TPM' not in expression:
                raise ValueError('RNA context requires stable gene_id and TPM columns')
            expr = expression[[idcol, 'TPM']].drop_duplicates()
            if expr.duplicated(idcol).any(): raise ValueError('Conflicting gene expression records')
            predictions['Expression'] = predictions.TargetGeneEnsemblID.map(expr.set_index(idcol).TPM)
            predictions['expression_status'] = np.where(predictions.Expression.isna(), 'unknown',
                np.where(predictions.Expression >= min_expression, 'detected', 'below_threshold'))
            predictions['isExpressed'] = predictions.Expression >= min_expression
        return predictions

    def filter_predictions(self,
                          predictions: pd.DataFrame,
                          threshold: float = 0.02,
                          only_expressed: bool = False) -> pd.DataFrame:
        """
        Filter predictions by threshold.
        
        Args:
            predictions: DataFrame with predictions
            threshold: Score threshold
            only_expressed: Only keep expressed genes
        
        Returns:
            Filtered predictions
        """
        filtered = predictions.copy()
        
        # Filter by score
        filtered = filtered[filtered[self.score_column] >= threshold]
        
        # Filter to expressed genes
        if only_expressed and 'isExpressed' in filtered.columns:
            filtered = filtered[filtered['isExpressed']]
        
        logger.info(f"Filtered to {len(filtered)} predictions (threshold={threshold})")
        
        return filtered


def run_predictions(enhancer_file: str,
                   gene_file: str,
                   output_file: str,
                   hic_file: Optional[str] = None,
                   hic_type: str = 'hic',
                   hic_resolution: int = 5000,
                   expression_file: Optional[str] = None,
                   use_expression_weight: bool = False,
                   weight_method: str = 'log',
                   min_expression: float = 1.0,
                   max_distance: int = 5000000,
                   hic_gamma: float = 1.024238616787792,
                   hic_scale: float = 5.9594510043736655,
                   score_threshold: float = 0.02,
                   include_self_promoter: bool = True,
                   contact_metadata_file: Optional[str] = None) -> str:
    """
    Run prediction pipeline.
    
    Args:
        enhancer_file: Path to enhancer list file
        gene_file: Path to gene list file
        output_file: Output predictions file
        hic_file: Optional Hi-C file
        hic_type: Type of Hi-C file
        hic_resolution: Hi-C resolution
        expression_file: Optional expression file
        use_expression_weight: Whether to use expression weight
        weight_method: Expression weight method
        min_expression: Minimum expression
        max_distance: Maximum E-G distance
        hic_gamma: Power-law exponent
        hic_scale: Power-law scale
        score_threshold: Score threshold for filtering
        include_self_promoter: Include self-promoter interactions
    
    Returns:
        Path to output file
    """
    # Load data
    logger.info("Loading enhancers and genes")
    
    enhancers = pd.read_csv(enhancer_file, sep='\t')
    genes = pd.read_csv(gene_file, sep='\t')
    
    # Load expression if provided
    expression = None
    if expression_file:
        expression = pd.read_csv(expression_file, sep='\t')
    
    # Initialize predictor
    contact_method = ('avg' if hic_type == 'avg' else 'hic') if hic_file else 'powerlaw'
    
    predictor = PACEPredictor(
        max_distance=max_distance,
        contact_method=contact_method,
        hic_file=hic_file,
        hic_type=hic_type,
        hic_resolution=hic_resolution,
        hic_gamma=hic_gamma,
        hic_scale=hic_scale,
        include_self_promoter=include_self_promoter
    )
    
    # Run predictions
    predictions = predictor.predict(
        enhancers=enhancers,
        genes=genes,
        expression=expression,
        expression_weight=use_expression_weight,
        weight_method=weight_method,
        min_expression=min_expression,
        contact_metadata=pd.read_csv(contact_metadata_file, sep='\t') if contact_metadata_file else None
    )
    
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    # Write all predictions (gzipped)
    if output_file.endswith('.gz'):
        predictions.to_csv(output_file, sep='\t', index=False, compression='gzip')
    else:
        predictions.to_csv(output_file, sep='\t', index=False)
    
    logger.info(f"Wrote {len(predictions)} predictions to {output_file}")
    
    # Also write filtered predictions
    destination = Path(output_file)
    filtered_file = str(destination.with_name(destination.name.replace('AllPutative', 'Filtered')))
    if filtered_file == output_file:
        filtered_file = output_file + '.filtered.tsv'
    # Normalize the extension to a single .tsv (avoid '.tsv.tsv').
    if filtered_file.endswith('.tsv.gz'):
        filtered_file = filtered_file[:-len('.tsv.gz')] + '.tsv'
    elif filtered_file.endswith('.gz'):
        filtered_file = filtered_file[:-len('.gz')]

    filtered = predictor.filter_predictions(predictions, threshold=score_threshold)
    filtered.to_csv(filtered_file, sep='\t', index=False)
    
    logger.info(f"Wrote {len(filtered)} filtered predictions to {filtered_file}")
    
    return output_file
