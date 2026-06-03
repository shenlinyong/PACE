#!/usr/bin/env python3
"""
PACE: Enhanced ABC Score Calculator

This script computes the enhanced PACE score by integrating:
1. Chromatin accessibility (ATAC-seq/DNase-seq)
2. Multiple histone modifications
3. 3D chromatin contact (Hi-C or power-law)
4. Gene expression (RNA-seq)
5. DNA methylation (inhibitory)
6. TF binding
7. eQTL validation

PACE Score Formula:
==================

    PACE_Score(E,G) = [Activity(E) × Contact(E,G) × Expr_Weight(G) × Reg_Factor(E)] / Σ(...)

Where:
    Activity(E) = Aggregation of chromatin signals at enhancer E
    Contact(E,G) = 3D contact frequency between E and gene G promoter
    Expr_Weight(G) = Expression-based weight for gene G
    Reg_Factor(E) = Regulatory factor (TF binding, methylation, etc.)

Activity Calculation:
====================

For geometric mean (original ABC):
    Activity = √(Accessibility × H3K27ac)

For weighted geometric mean (PACE):
    Activity = ∏(Signal_i ^ Weight_i) × (1 - Inhibitory_Score)

Where Inhibitory_Score combines repressive marks and methylation.

Author: Linyong Shen @ Northwest A&F University
"""

import argparse
import numpy as np
import pandas as pd
import logging
from typing import Dict, List, Optional, Tuple
import os
import sys

# Import from multiomics module
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from multiomics_activity import (
    MultiOmicsActivityCalculator,
    SignalType,
    AggregationMethod,
    ExpressionFilter,
    eQTLValidator,
    MethylationIntegrator,
)

# Reuse the real Hi-C contact estimator from the workflow implementation so
# that the standalone calculator and the Snakemake pipeline share identical
# contact logic (.hic via hic-straw, .cool via cooler, BEDPE, or power-law).
_WORKFLOW_SCRIPTS = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "workflow", "scripts",
)
if os.path.isdir(_WORKFLOW_SCRIPTS):
    sys.path.insert(0, _WORKFLOW_SCRIPTS)
try:
    from hic import ContactEstimator
    HAS_HIC = True
except ImportError:  # pragma: no cover
    ContactEstimator = None
    HAS_HIC = False

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class PACEScoreCalculator:
    """
    Enhanced ABC score calculator with multi-omics integration.
    
    This is the main class for PACE, implementing flexible
    integration of diverse epigenomic data types.
    """
    
    def __init__(self, config: Dict):
        """
        Initialize PACE calculator.
        
        Args:
            config: Configuration dictionary
        """
        self.config = config
        self.activity_method = config.get('activity_method', 'geometric_mean')
        
        # Initialize components
        self._init_activity_calculator()
        self._init_expression_filter()
        self._init_eqtl_validator()
        
    def _init_activity_calculator(self):
        """Initialize multi-omics activity calculator"""
        method_map = {
            'geometric_mean': AggregationMethod.GEOMETRIC_MEAN,
            'weighted_geometric': AggregationMethod.WEIGHTED_GEOMETRIC,
            'weighted_sum': AggregationMethod.WEIGHTED_SUM,
            'arithmetic_mean': AggregationMethod.ARITHMETIC_MEAN
        }
        method = method_map.get(self.activity_method, AggregationMethod.GEOMETRIC_MEAN)
        self.activity_calc = MultiOmicsActivityCalculator(method=method)
        
    def _init_expression_filter(self):
        """Initialize expression filter if enabled"""
        expr_config = self.config.get('expression', {})
        if expr_config.get('enabled', False):
            self.expression_filter = ExpressionFilter(
                expression_file=expr_config.get('file', ''),
                min_expression=expr_config.get('min_expression', 1.0),
                expression_column=expr_config.get('value_column', 'TPM')
            )
        else:
            self.expression_filter = None
            
    def _init_eqtl_validator(self):
        """Initialize eQTL validator if enabled"""
        eqtl_config = self.config.get('eqtl_validation', {})
        if eqtl_config.get('enabled', False):
            cols = eqtl_config.get('columns', {})
            self.eqtl_validator = eQTLValidator(
                eqtl_file=eqtl_config.get('file', ''),
                chr_col=cols.get('chr', 'chr'),
                pos_col=cols.get('pos', 'pos'),
                gene_col=cols.get('gene', 'gene'),
                pval_col=cols.get('pvalue', 'pvalue'),
                beta_col=cols.get('beta', 'beta'),
            )
        else:
            self.eqtl_validator = None
        self.eqtl_enrichment = None
    
    def add_signals_from_sample(self, sample_config: Dict):
        """
        Add all available signals from sample configuration.
        
        Args:
            sample_config: Sample-specific configuration
        """
        # Accessibility (required)
        if sample_config.get('ATAC'):
            self.activity_calc.add_signal(
                "ATAC", SignalType.ACCESSIBILITY, 
                sample_config['ATAC'],
                weight=self.config.get('accessibility', {}).get('weight', 1.5)
            )
        elif sample_config.get('DHS'):
            self.activity_calc.add_signal(
                "DHS", SignalType.ACCESSIBILITY,
                sample_config['DHS'],
                weight=self.config.get('accessibility', {}).get('weight', 1.5)
            )
        
        # Histone modifications
        histone_config = self.config.get('histone_marks', {})
        
        histone_mapping = {
            'H3K27ac': SignalType.ACTIVE_ENHANCER,
            'H3K4me1': SignalType.ENHANCER_MARK,
            'H3K4me3': SignalType.PROMOTER_MARK,
            'H3K36me3': SignalType.TRANSCRIPTION,
            'H3K9ac': SignalType.ACTIVE_ENHANCER,
            'H3K27me3': SignalType.REPRESSIVE,
            'H3K9me3': SignalType.REPRESSIVE,
        }
        
        for mark, signal_type in histone_mapping.items():
            if sample_config.get(mark):
                mark_config = histone_config.get(mark, {})
                if mark_config.get('enabled', True):
                    # Repressive marks (H3K27me3/H3K9me3) are inhibitory by
                    # default; this can be overridden in the config.
                    default_inhib = signal_type == SignalType.REPRESSIVE
                    self.activity_calc.add_signal(
                        mark, signal_type,
                        sample_config[mark],
                        weight=mark_config.get('weight', 1.0),
                        inhibitory=mark_config.get('inhibitory', default_inhib)
                    )
        
        # Methylation
        if sample_config.get('methylation'):
            meth_config = self.config.get('methylation', {})
            if meth_config.get('enabled', True):
                self.activity_calc.add_signal(
                    "methylation", SignalType.METHYLATION,
                    sample_config['methylation'],
                    weight=meth_config.get('weight', 0.5),
                    inhibitory=True
                )
        
        # TF binding. Transcription factors can be activators or repressors;
        # PACE does not assume a fixed sign. The mode is
        # taken from the per-sample 'TF_modes' column (activator/repressor) or
        # from the 'specific_tfs' config; the default is activator.
        if sample_config.get('TF_binding'):
            tf_files = str(sample_config['TF_binding']).split(',')
            tf_names = str(sample_config.get('TF_names', '')).split(',')
            tf_modes = str(sample_config.get('TF_modes', '')).split(',')

            tf_config = self.config.get('transcription_factors', {})
            default_weight = tf_config.get('default_weight', 0.3)
            specific_tfs = tf_config.get('specific_tfs', {})

            for i, (tf_file, tf_name) in enumerate(zip(tf_files, tf_names)):
                tf_name = tf_name.strip()
                tf_file = tf_file.strip()
                if not tf_file:
                    continue

                spec = specific_tfs.get(tf_name, {})
                weight = spec.get('weight', default_weight)

                # Determine activator vs repressor.
                mode = tf_modes[i].strip().lower() if i < len(tf_modes) else ''
                if mode in ('repressor', 'inhibitory', 'repressive'):
                    is_inhibitory = True
                elif mode in ('activator', 'activating'):
                    is_inhibitory = False
                else:
                    is_inhibitory = bool(spec.get('inhibitory', False))

                self.activity_calc.add_signal(
                    tf_name, SignalType.TF_BINDING,
                    tf_file,
                    weight=weight,
                    inhibitory=is_inhibitory,
                )
    
    def calculate_activity(self, 
                           candidate_regions: pd.DataFrame) -> Tuple[np.ndarray, Dict]:
        """
        Calculate activity scores for candidate enhancer regions.
        
        Args:
            candidate_regions: DataFrame with enhancer regions
            
        Returns:
            Tuple of (activity_scores, component_scores)
        """
        logger.info(f"Calculating activity for {len(candidate_regions)} regions")
        
        result = self.activity_calc.calculate_activity(
            candidate_regions, 
            return_components=True
        )
        
        return result['activity'], result
    
    def calculate_contact(self,
                          enhancers: pd.DataFrame,
                          genes: pd.DataFrame,
                          hic_file: str = None,
                          hic_type: str = None,
                          hic_resolution: int = 5000) -> pd.DataFrame:
        """
        Calculate contact frequencies between enhancers and gene promoters.
        
        Args:
            enhancers: DataFrame with enhancer regions
            genes: DataFrame with gene annotations
            hic_file: Path to Hi-C file (optional)
            hic_type: Type of Hi-C file
            hic_resolution: Hi-C resolution
            
        Returns:
            DataFrame with E-G pairs and contact frequencies
        """
        # Build all candidate E-G pairs within max_distance and attach the
        # genomic distance; contact is then filled either from Hi-C or the
        # power-law model.
        pairs = self._create_pairs(enhancers, genes)
        if len(pairs) == 0:
            return pairs

        predict_config = self.config.get('params_predict', {})
        gamma = predict_config.get('hic_gamma', 1.024238616787792)
        scale = predict_config.get('hic_scale', 5.9594510043736655)
        pseudocount = predict_config.get('hic_pseudocount_distance', 5000)

        if hic_file and HAS_HIC:
            logger.info("Calculating contact from Hi-C data (%s)", hic_type)
            estimator = ContactEstimator(
                method='hic', hic_file=hic_file, hic_type=hic_type,
                resolution=hic_resolution, hic_gamma=gamma, hic_scale=scale,
            )
            pairs['contact'] = estimator.estimate_batch(
                pairs, chrom_col='enhancer_chr',
                pos1_col='enhancer_mid', pos2_col='gene_tss')
        else:
            if hic_file and not HAS_HIC:
                logger.warning("hic module unavailable; using power-law")
            logger.info("Calculating contact using power-law")
            distance = pairs['distance'].to_numpy(dtype=float)
            pairs['contact'] = scale / np.power(distance + pseudocount, gamma)

        return pairs

    def _create_pairs(self,
                      enhancers: pd.DataFrame,
                      genes: pd.DataFrame,
                      max_distance: int = 5000000) -> pd.DataFrame:
        """Create all enhancer-gene pairs within ``max_distance``."""
        pairs = []
        for _, enh in enhancers.iterrows():
            enh_center = (enh['start'] + enh['end']) // 2
            nearby_genes = genes[
                (genes['chr'] == enh['chr']) &
                (abs(genes['tss'] - enh_center) <= max_distance)
            ]
            for _, gene in nearby_genes.iterrows():
                distance = abs(gene['tss'] - enh_center)
                pairs.append({
                    'enhancer_chr': enh['chr'],
                    'enhancer_start': enh['start'],
                    'enhancer_end': enh['end'],
                    'enhancer_mid': enh_center,
                    'enhancer_name': enh.get('name', ''),
                    'gene_id': gene['gene_id'],
                    'gene_name': gene['gene_name'],
                    'gene_tss': gene['tss'],
                    'distance': distance,
                })
        return pd.DataFrame(pairs)
    
    def calculate_pace_scores(self,
                              candidate_regions: pd.DataFrame,
                              genes: pd.DataFrame,
                              sample_config: Dict) -> pd.DataFrame:
        """
        Calculate PACE scores for all E-G pairs.
        
        This is the main method that integrates all components.
        
        Args:
            candidate_regions: Candidate enhancer regions
            genes: Gene annotations
            sample_config: Sample-specific configuration
            
        Returns:
            DataFrame with PACE predictions
        """
        logger.info("=" * 60)
        logger.info("PACE Score Calculation")
        logger.info("=" * 60)
        
        # Add signals from sample
        self.add_signals_from_sample(sample_config)

        # Ensure every candidate region has a unique name (used for the
        # contact pairing and the component-score broadcast below).
        candidate_regions = candidate_regions.copy()
        if 'name' not in candidate_regions.columns:
            candidate_regions['name'] = [
                f"{r['chr']}:{r['start']}-{r['end']}"
                for _, r in candidate_regions.iterrows()]

        # Calculate activity
        logger.info("\n[1/5] Calculating enhancer activity...")
        activity, activity_components = self.calculate_activity(candidate_regions)
        candidate_regions['activity'] = activity
        
        # Calculate contact
        logger.info("\n[2/5] Calculating E-G contact frequencies...")
        eg_pairs = self.calculate_contact(
            candidate_regions, genes,
            hic_file=sample_config.get('HiC_file'),
            hic_type=sample_config.get('HiC_type'),
            hic_resolution=sample_config.get('HiC_resolution', 5000)
        )
        
        # Merge activity with E-G pairs
        eg_pairs = eg_pairs.merge(
            candidate_regions[['chr', 'start', 'end', 'activity']],
            left_on=['enhancer_chr', 'enhancer_start', 'enhancer_end'],
            right_on=['chr', 'start', 'end'],
            how='left'
        )
        
        # Calculate raw ABC score
        logger.info("\n[3/5] Calculating raw ABC scores...")
        eg_pairs['activity_x_contact'] = eg_pairs['activity'] * eg_pairs['contact']
        
        # Normalize by sum per gene
        gene_sums = eg_pairs.groupby('gene_id')['activity_x_contact'].sum()
        eg_pairs['gene_sum'] = eg_pairs['gene_id'].map(gene_sums)
        eg_pairs['ABC.Score'] = eg_pairs['activity_x_contact'] / eg_pairs['gene_sum']
        
        # Apply expression filter/weight
        logger.info("\n[4/5] Applying expression filter/weight...")
        if self.expression_filter and sample_config.get('RNA_seq'):
            self.expression_filter.expression_file = sample_config['RNA_seq']
            eg_pairs = self.expression_filter.filter_predictions(eg_pairs)
            
            weight_method = self.config.get('expression', {}).get('weight_method', 'none')
            if weight_method != 'none':
                eg_pairs = self.expression_filter.add_expression_weight(eg_pairs, weight_method)
                eg_pairs['PACE.Score'] = eg_pairs['ABC.Score'] * eg_pairs['expression_weight']
            else:
                eg_pairs['PACE.Score'] = eg_pairs['ABC.Score']
        else:
            eg_pairs['PACE.Score'] = eg_pairs['ABC.Score']
        
        # Add eQTL validation
        logger.info("\n[5/5] Adding eQTL validation...")
        if self.eqtl_validator and sample_config.get('eQTL_file'):
            eqtl_cfg = self.config.get('eqtl_validation', {})
            self.eqtl_validator.eqtl_file = sample_config['eQTL_file']
            eg_pairs = self.eqtl_validator.validate_predictions(
                eg_pairs,
                pvalue_threshold=eqtl_cfg.get('pvalue_threshold', 1e-5),
                window=eqtl_cfg.get('window', 1000),
            )
            try:
                self.eqtl_enrichment = self.eqtl_validator.calculate_enrichment(
                    eg_pairs, score_column='PACE.Score',
                    threshold=self.config.get('params_filter_predictions', {})
                    .get('threshold', 0.02) or 0.02)
                logger.info("eQTL enrichment OR=%.2f (P=%.2e)",
                            self.eqtl_enrichment['odds_ratio'],
                            self.eqtl_enrichment['pvalue'])
            except Exception as exc:  # pragma: no cover
                logger.warning("Could not compute eQTL enrichment: %s", exc)

        # Attach per-signal component scores to each E-G pair (one value per
        # enhancer, broadcast to all of that enhancer's pairs).
        if self.config.get('output_options', {}).get('output_signal_scores', True):
            comp = activity_components.get('activating_signals', {})
            comp.update(activity_components.get('inhibitory_signals', {}))
            for signal_name, values in comp.items():
                sig_map = dict(zip(candidate_regions['name'], values))
                eg_pairs[signal_name] = eg_pairs['enhancer_name'].map(sig_map)

        # Standardize to the canonical PACE output schema so that the output
        # of the standalone calculator matches the Snakemake pipeline and is
        # directly consumable by the ML and validation modules.
        eg_pairs = eg_pairs.rename(columns={
            'enhancer_name': 'name',
            'gene_name': 'TargetGene',
            'gene_id': 'TargetGeneEnsemblID',
            'gene_tss': 'TargetGeneTSS',
        })
        if 'class' not in eg_pairs.columns:
            from tools import assign_enhancer_class
            eg_pairs['class'] = eg_pairs['distance'].apply(assign_enhancer_class)

        # Final cleanup and sorting
        eg_pairs = eg_pairs.sort_values('PACE.Score', ascending=False)

        logger.info("\n" + "=" * 60)
        logger.info(f"Completed! Generated {len(eg_pairs)} E-G predictions")
        logger.info("=" * 60)

        return eg_pairs


def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description='PACE: Enhanced ABC Score Calculator',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Basic usage
    python calculate_pace_score.py \\
        --config config/config_v1.yaml \\
        --sample Pig_Liver \\
        --output results/pace_predictions.tsv
        
    # With all options
    python calculate_pace_score.py \\
        --config config/config_v1.yaml \\
        --sample Pig_Liver \\
        --output results/pace_predictions.tsv \\
        --threshold 0.02 \\
        --output_components
        """
    )
    
    parser.add_argument('--config', required=True, help='Configuration file')
    parser.add_argument('--sample', required=True, help='Sample name')
    parser.add_argument('--candidates', required=True, help='Candidate regions BED file')
    parser.add_argument('--genes', required=True, help='Genes BED file')
    parser.add_argument('--output', required=True, help='Output file')
    parser.add_argument('--threshold', type=float, default=0.02, help='Score threshold')
    parser.add_argument('--output_components', action='store_true', 
                        help='Output individual signal scores')
    
    args = parser.parse_args()
    
    # Load configuration
    import yaml
    with open(args.config) as f:
        config = yaml.safe_load(f)
    
    # Load candidate regions
    candidates = pd.read_csv(args.candidates, sep='\t', 
                             names=['chr', 'start', 'end', 'name', 'score', 'strand'])
    
    # Load genes
    genes = pd.read_csv(args.genes, sep='\t',
                        names=['chr', 'start', 'end', 'name', 'score', 'strand', 
                               'gene_id', 'gene_type'])
    genes['gene_name'] = genes['name'].str.split(';').str[0]
    genes['tss'] = np.where(genes['strand'] == '+', genes['start'], genes['end'])
    
    # Get sample configuration
    biosamples = pd.read_csv(config['biosamplesTable'], sep='\t', comment='#')
    sample_config = biosamples[biosamples['biosample'] == args.sample].iloc[0].to_dict()
    
    # Calculate PACE scores
    calculator = PACEScoreCalculator(config)
    predictions = calculator.calculate_pace_scores(candidates, genes, sample_config)
    
    # Filter and save
    filtered = predictions[predictions['PACE.Score'] >= args.threshold]
    filtered.to_csv(args.output, sep='\t', index=False)
    
    logger.info(f"Saved {len(filtered)} predictions to {args.output}")


if __name__ == '__main__':
    main()
