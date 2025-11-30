#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PACE Predictor Module

Core prediction functions for enhancer-gene regulatory interactions.

Author: Linyong Shen @ Northwest A&F University
"""

import os
import numpy as np
import pandas as pd
from typing import Optional, List, Dict, Tuple, Union

from tools import (
    calculate_distance, power_law_contact, safe_divide,
    create_enhancer_gene_pairs, assign_enhancer_class, logger
)
from hic import ContactEstimator


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
                 score_column: str = 'ABC.Score',
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
    
    def create_pairs(self,
                    enhancers: pd.DataFrame,
                    genes: pd.DataFrame) -> pd.DataFrame:
        """
        Create all enhancer-gene pairs within distance threshold.
        
        Args:
            enhancers: DataFrame with enhancer regions
            genes: DataFrame with gene info
        
        Returns:
            DataFrame with E-G pairs
        """
        logger.info(f"Creating E-G pairs within {self.max_distance/1e6:.1f} Mb")
        
        pairs = []
        
        # Ensure required columns
        if 'TSS' not in genes.columns:
            if 'start' in genes.columns:
                genes['TSS'] = genes['start']
            else:
                logger.error("Genes must have TSS or start column")
                return pd.DataFrame()
        
        # Group by chromosome for efficiency
        for chrom in enhancers['chr'].unique():
            enh_chrom = enhancers[enhancers['chr'] == chrom]
            gene_chrom = genes[genes['chr'] == chrom]
            
            for _, enh in enh_chrom.iterrows():
                enh_mid = (enh['start'] + enh['end']) // 2
                
                for _, gene in gene_chrom.iterrows():
                    tss = gene['TSS']
                    dist = abs(enh_mid - tss)
                    
                    if dist <= self.max_distance:
                        pairs.append({
                            'chr': chrom,
                            'start': enh['start'],
                            'end': enh['end'],
                            'name': enh.get('name', f"{chrom}:{enh['start']}-{enh['end']}"),
                            'TargetGene': gene.get('gene_name', gene.get('name', '')),
                            'TargetGeneEnsemblID': gene.get('gene_id', ''),
                            'TargetGeneTSS': tss,
                            'TargetGeneStrand': gene.get('strand', '+'),
                            'distance': dist,
                            'enhancer_mid': enh_mid,
                            # Copy activity and signals if available
                            'activity': enh.get('activity', 1.0),
                        })
        
        df = pd.DataFrame(pairs)
        
        # Add enhancer class
        if len(df) > 0:
            df['class'] = df.apply(
                lambda x: assign_enhancer_class(x['distance']),
                axis=1
            )
        
        logger.info(f"Created {len(df)} E-G pairs")
        
        return df
    
    def estimate_contacts(self, pairs: pd.DataFrame) -> pd.DataFrame:
        """
        Estimate 3D contacts for all pairs.
        
        Args:
            pairs: DataFrame with E-G pairs
        
        Returns:
            DataFrame with contact column added
        """
        logger.info("Estimating 3D contacts")
        
        pairs = pairs.copy()
        
        contacts = self.contact_estimator.estimate_batch(
            pairs,
            chrom_col='chr',
            pos1_col='enhancer_mid',
            pos2_col='TargetGeneTSS'
        )
        
        pairs['contact'] = contacts
        
        return pairs
    
    def calculate_abc_score(self,
                           pairs: pd.DataFrame,
                           activity_column: str = 'activity',
                           contact_column: str = 'contact') -> pd.DataFrame:
        """
        Calculate ABC/PACE score.
        
        Args:
            pairs: DataFrame with E-G pairs
            activity_column: Column with activity values
            contact_column: Column with contact values
        
        Returns:
            DataFrame with score column added
        """
        logger.info("Calculating prediction scores")
        
        pairs = pairs.copy()
        
        # Calculate A × C product
        pairs['AxC'] = pairs[activity_column] * pairs[contact_column]
        
        # Normalize by sum per gene
        gene_sums = pairs.groupby('TargetGene')['AxC'].transform('sum')
        
        pairs[self.score_column] = safe_divide(pairs['AxC'], gene_sums, default=0.0)
        
        # Handle self-promoters
        if not self.include_self_promoter:
            pairs.loc[pairs['class'] == 'promoter', self.score_column] = 0.0
        
        # Clean up
        pairs.drop('AxC', axis=1, inplace=True)
        
        logger.info(f"Score range: {pairs[self.score_column].min():.4f} - {pairs[self.score_column].max():.4f}")
        
        return pairs
    
    def add_expression_weight(self,
                             pairs: pd.DataFrame,
                             expression: pd.DataFrame,
                             gene_column: str = 'TargetGene',
                             expr_column: str = 'TPM',
                             weight_method: str = 'log',
                             min_expression: float = 1.0) -> pd.DataFrame:
        """
        Add expression weight to predictions.
        
        Args:
            pairs: DataFrame with E-G pairs
            expression: DataFrame with gene expression
            gene_column: Column with gene names
            expr_column: Column with expression values
            weight_method: Weight method (binary, linear, log)
            min_expression: Minimum expression threshold
        
        Returns:
            DataFrame with expression weight
        """
        logger.info("Adding expression weights")
        
        pairs = pairs.copy()
        
        # Create expression dictionary
        expr_dict = dict(zip(expression.iloc[:, 0], expression[expr_column]))
        
        # Map expression to genes
        pairs['Expression'] = pairs[gene_column].map(expr_dict).fillna(0)
        
        # Filter to expressed genes
        pairs['isExpressed'] = pairs['Expression'] >= min_expression
        
        # Calculate expression weight
        if weight_method == 'binary':
            pairs['ExpressionWeight'] = pairs['isExpressed'].astype(float)
        
        elif weight_method == 'linear':
            max_expr = pairs['Expression'].max()
            pairs['ExpressionWeight'] = pairs['Expression'] / (max_expr + 1e-10)
        
        elif weight_method == 'log':
            log_expr = np.log2(pairs['Expression'] + 1)
            max_log = log_expr.max()
            pairs['ExpressionWeight'] = log_expr / (max_log + 1e-10)
        
        else:
            pairs['ExpressionWeight'] = 1.0
        
        logger.info(f"Expressed genes: {pairs['isExpressed'].sum()}")
        
        return pairs
    
    def predict(self,
               enhancers: pd.DataFrame,
               genes: pd.DataFrame,
               expression: Optional[pd.DataFrame] = None,
               expression_weight: bool = False,
               weight_method: str = 'log',
               min_expression: float = 1.0) -> pd.DataFrame:
        """
        Run full prediction pipeline.
        
        Args:
            enhancers: DataFrame with enhancers (must have 'activity' column)
            genes: DataFrame with genes
            expression: Optional expression DataFrame
            expression_weight: Whether to weight by expression
            weight_method: Expression weight method
            min_expression: Minimum expression threshold
        
        Returns:
            DataFrame with predictions
        """
        # Create pairs
        pairs = self.create_pairs(enhancers, genes)
        
        if len(pairs) == 0:
            logger.warning("No E-G pairs created!")
            return pd.DataFrame()
        
        # Copy activity from enhancers (handle different column names)
        activity_col = None
        for col in ['activity', 'activity_base', 'Activity', 'score', 'signalValue']:
            if col in enhancers.columns:
                activity_col = col
                break
        
        if activity_col:
            enh_activity = dict(zip(enhancers['name'], enhancers[activity_col]))
            pairs['activity'] = pairs['name'].map(enh_activity).fillna(1.0)
        else:
            logger.warning("No activity column found in enhancers, using default 1.0")
            pairs['activity'] = 1.0
        
        # Estimate contacts
        pairs = self.estimate_contacts(pairs)
        
        # Add expression weight if requested
        if expression_weight and expression is not None:
            pairs = self.add_expression_weight(
                pairs, expression,
                weight_method=weight_method,
                min_expression=min_expression
            )
            
            # Incorporate expression weight into activity
            pairs['activity_weighted'] = pairs['activity'] * pairs['ExpressionWeight']
            
            # Calculate score with expression weight
            pairs = self.calculate_abc_score(
                pairs,
                activity_column='activity_weighted',
                contact_column='contact'
            )
        else:
            # Calculate score without expression weight
            pairs = self.calculate_abc_score(pairs)
        
        # Sort by score
        pairs = pairs.sort_values(self.score_column, ascending=False)
        
        return pairs
    
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
                   include_self_promoter: bool = True) -> str:
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
    if expression_file and os.path.exists(expression_file):
        expression = pd.read_csv(expression_file, sep='\t')
    
    # Initialize predictor
    contact_method = 'hic' if hic_file else 'powerlaw'
    
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
        min_expression=min_expression
    )
    
    # Write all predictions (gzipped)
    if output_file.endswith('.gz'):
        predictions.to_csv(output_file, sep='\t', index=False, compression='gzip')
    else:
        predictions.to_csv(output_file, sep='\t', index=False)
    
    logger.info(f"Wrote {len(predictions)} predictions to {output_file}")
    
    # Also write filtered predictions
    filtered_file = output_file.replace('AllPutative', 'Filtered')
    filtered_file = filtered_file.replace('.gz', '.tsv')
    
    filtered = predictor.filter_predictions(predictions, threshold=score_threshold)
    filtered.to_csv(filtered_file, sep='\t', index=False)
    
    logger.info(f"Wrote {len(filtered)} filtered predictions to {filtered_file}")
    
    return output_file
