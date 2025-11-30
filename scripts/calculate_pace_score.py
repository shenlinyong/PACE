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
try:
    from multiomics_activity import (
        MultiOmicsActivityCalculator,
        SignalType,
        AggregationMethod,
        ExpressionFilter,
        eQTLValidator,
        MethylationIntegrator
    )
except ImportError:
    # Fallback if running standalone
    pass

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
            self.eqtl_validator = eQTLValidator(
                eqtl_file=eqtl_config.get('file', '')
            )
        else:
            self.eqtl_validator = None
    
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
                    self.activity_calc.add_signal(
                        mark, signal_type,
                        sample_config[mark],
                        weight=mark_config.get('weight', 1.0),
                        inhibitory=mark_config.get('inhibitory', False)
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
        
        # TF binding
        if sample_config.get('TF_binding'):
            tf_files = sample_config['TF_binding'].split(',')
            tf_names = sample_config.get('TF_names', '').split(',')
            
            tf_config = self.config.get('transcription_factors', {})
            default_weight = tf_config.get('default_weight', 0.3)
            specific_tfs = tf_config.get('specific_tfs', {})
            
            for tf_file, tf_name in zip(tf_files, tf_names):
                tf_name = tf_name.strip()
                tf_file = tf_file.strip()
                
                weight = specific_tfs.get(tf_name, {}).get('weight', default_weight)
                
                self.activity_calc.add_signal(
                    tf_name, SignalType.TF_BINDING,
                    tf_file,
                    weight=weight
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
        contact_config = self.config.get('params_contact', {})
        
        if hic_file and contact_config.get('method') == 'hic':
            logger.info("Calculating contact from Hi-C data")
            # Implementation would use straw or cooler
            contact = self._calculate_hic_contact(
                enhancers, genes, hic_file, hic_type, hic_resolution
            )
        else:
            logger.info("Calculating contact using power-law")
            contact = self._calculate_powerlaw_contact(enhancers, genes)
        
        return contact
    
    def _calculate_powerlaw_contact(self,
                                    enhancers: pd.DataFrame,
                                    genes: pd.DataFrame) -> pd.DataFrame:
        """Calculate contact using power-law distance decay"""
        contact_config = self.config.get('params_contact', {})
        gamma = contact_config.get('powerlaw_gamma', 1.024238616787792)
        scale = contact_config.get('powerlaw_scale', 5.9594510043736655)
        pseudocount = contact_config.get('pseudocount_distance', 5000)
        
        # Create E-G pairs within 5Mb
        pairs = []
        
        for _, enh in enhancers.iterrows():
            enh_center = (enh['start'] + enh['end']) // 2
            
            # Find genes within 5Mb
            nearby_genes = genes[
                (genes['chr'] == enh['chr']) &
                (abs(genes['tss'] - enh_center) <= 5000000)
            ]
            
            for _, gene in nearby_genes.iterrows():
                distance = abs(gene['tss'] - enh_center)
                
                # Power-law contact: C = scale / (distance + pseudocount)^gamma
                contact = scale / ((distance + pseudocount) ** gamma)
                
                pairs.append({
                    'enhancer_chr': enh['chr'],
                    'enhancer_start': enh['start'],
                    'enhancer_end': enh['end'],
                    'enhancer_name': enh.get('name', ''),
                    'gene_id': gene['gene_id'],
                    'gene_name': gene['gene_name'],
                    'gene_tss': gene['tss'],
                    'distance': distance,
                    'contact': contact
                })
        
        return pd.DataFrame(pairs)
    
    def _calculate_hic_contact(self,
                               enhancers: pd.DataFrame,
                               genes: pd.DataFrame,
                               hic_file: str,
                               hic_type: str,
                               resolution: int) -> pd.DataFrame:
        """Calculate contact from Hi-C data"""
        # Implementation would use hic-straw or cooler
        # Placeholder for now
        logger.warning("Hi-C contact calculation not fully implemented, using power-law")
        return self._calculate_powerlaw_contact(enhancers, genes)
    
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
        
        # Calculate activity
        logger.info("\n[1/5] Calculating enhancer activity...")
        activity, activity_components = self.calculate_activity(candidate_regions)
        candidate_regions = candidate_regions.copy()
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
            self.eqtl_validator.eqtl_file = sample_config['eQTL_file']
            eg_pairs = self.eqtl_validator.validate_predictions(eg_pairs)
        
        # Add component scores if requested
        if self.config.get('output_options', {}).get('output_signal_scores', True):
            for signal_name, values in activity_components.get('activating_signals', {}).items():
                # Would need to properly merge these
                pass
        
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
