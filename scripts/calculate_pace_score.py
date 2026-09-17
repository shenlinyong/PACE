#!/usr/bin/env python3
"""
PACE signal-file scoring interface.

Quantifies accessibility and configured activating signals, then scores
enhancer-gene candidates using the shared activity-contact kernel.

PACE uses the single kernel in workflow/scripts/pace_core.py.
RNA is context only; see docs/FORMULA.md for the versioned equations.

Author: 申林用 (Linyong Shen) @ Northwest A&F University
"""

import argparse
import numpy as np
import pandas as pd
import logging
from typing import Dict, List, Optional, Tuple
import os
import sys
from pathlib import Path

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
from tools import read_candidate_regions
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
        self.activity_method = config.get('activity_method', 'missing_geometric')
        
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
        if self.activity_method != 'missing_geometric':
            raise ValueError('PACE requires activity_method: missing_geometric; use archived v1 for old modes')
        method = AggregationMethod.WEIGHTED_GEOMETRIC
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
                weight=self.config.get('accessibility', {}).get('weight', 1.0)
            )
        elif sample_config.get('DHS'):
            self.activity_calc.add_signal(
                "DHS", SignalType.ACCESSIBILITY,
                sample_config['DHS'],
                weight=self.config.get('accessibility', {}).get('weight', 1.0)
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
                if mark_config.get('enabled', mark == 'H3K27ac'):
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
            if meth_config.get('enabled', False):
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
        if sample_config.get('TF_binding') and self.config.get('transcription_factors', {}).get('enabled', False):
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
            enh_center = (enh['start'] + enh['end']) / 2
            nearby_genes = genes[
                (genes['chr'] == enh['chr']) &
                (abs(genes['tss'] - enh_center) < max_distance)
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
    
    def calculate_pace_scores(self, candidate_regions, genes, sample_config):
        from predictor import PACEPredictor
        self.activity_calc.signals.clear()
        # Empty TSV cells are absent assays, not float NaN paths.
        sample_config = {k: v for k, v in sample_config.items()
                         if v is not None and not (isinstance(v, float) and np.isnan(v))}
        self.add_signals_from_sample(sample_config)
        activity, components = self.calculate_activity(candidate_regions)
        enh = candidate_regions.copy()
        enh['activity'] = activity
        if hasattr(self.activity_calc, 'last_qc'):
            for c in self.activity_calc.last_qc:
                if c != 'activity': enh[c] = self.activity_calc.last_qc[c].to_numpy()
        g = genes.copy()
        if 'TSS' not in g and 'tss' in g: g['TSS'] = g.tss
        hic = sample_config.get('HiC_file')
        params = self.config.get('params_predict', {})
        predictor = PACEPredictor(
            max_distance=params.get('window', 5000000),
            contact_method=('avg' if sample_config.get('HiC_type') == 'avg' else 'hic') if hic else 'powerlaw',
            hic_file=hic, hic_type=sample_config.get('HiC_type', 'hic'),
            hic_resolution=int(sample_config.get('HiC_resolution', 5000)),
            hic_gamma=params.get('hic_gamma', 1.024238616787792),
            hic_scale=params.get('hic_scale', 5.9594510043736655))
        expr_path = sample_config.get('RNA_seq')
        expression = pd.read_csv(expr_path, sep='\t') if expr_path else None
        qc_path = sample_config.get('contact_metadata')
        qc = pd.read_csv(qc_path, sep='\t') if qc_path else None
        return predictor.predict(enh, g, expression=expression, contact_metadata=qc)


def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description='PACE: score enhancer-gene candidates from signal files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Basic usage
    python calculate_pace_score.py \\
        --config config/config.yaml \\
        --sample Pig_Liver \\
        --candidates candidates.bed --genes genes.tsv \\
        --output results/pace_predictions.tsv
        
    # With all options
    python calculate_pace_score.py \\
        --config config/config.yaml \\
        --sample Pig_Liver \\
        --candidates candidates.bed --genes genes.tsv \\
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
    candidates = read_candidate_regions(args.candidates)
    
    # Accept the current all-TSS TSV as well as historical extended BED.
    with open(args.genes) as handle:
        header = handle.readline().split('\t')[0].lower()
    if header in ('chr', 'chrom', '#chr'):
        genes = pd.read_csv(args.genes, sep='\t')
    else:
        genes = pd.read_csv(args.genes, sep='\t',
                            names=['chr','start','end','name','score','strand','gene_id','gene_type'])
    if 'gene_name' not in genes:
        genes['gene_name'] = genes['name'].str.split(';').str[0] if 'name' in genes else genes['gene_id']
    if 'TSS' not in genes and 'tss' not in genes:
        genes['TSS'] = np.where(genes['strand'] == '+', genes['start'], genes['end'] - 1)

    # Get sample configuration
    biosamples = pd.read_csv(config['biosamplesTable'], sep='\t', comment='#')
    selected = biosamples[biosamples['biosample'] == args.sample]
    if len(selected) != 1:
        parser.error(f'Sample {args.sample!r} must occur exactly once in biosamplesTable')
    sample_config = selected.iloc[0].to_dict()
    
    # Calculate PACE scores
    calculator = PACEScoreCalculator(config)
    predictions = calculator.calculate_pace_scores(candidates, genes, sample_config)
    
    # Filter and save
    filtered = predictions[predictions['PACE.Score'] >= args.threshold]
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    predictions.to_csv(args.output, sep='\t', index=False)
    filtered.to_csv(args.output + '.filtered.tsv', sep='\t', index=False)
    
    logger.info(f"Saved {len(predictions)} predictions to {args.output}; {len(filtered)} to the filtered file")


if __name__ == '__main__':
    main()
