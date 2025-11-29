#!/usr/bin/env python3
"""
PACE: Multi-omics Activity Calculator

This module provides flexible activity score calculation by integrating
multiple epigenomic signals with configurable weights.

Supported data types:
- Chromatin accessibility: ATAC-seq, DNase-seq
- Histone modifications: H3K27ac, H3K4me1, H3K4me3, H3K9ac, H3K36me3, etc.
- Transcription factors: TF ChIP-seq, motif scores
- DNA methylation: WGBS, RRBS
- Gene expression: RNA-seq

Author: Linyong Shen @ Northwest A&F University
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Union
from dataclasses import dataclass, field
from enum import Enum
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class AggregationMethod(Enum):
    """Methods for aggregating multiple signals"""
    GEOMETRIC_MEAN = "geometric_mean"      # Original ABC method
    ARITHMETIC_MEAN = "arithmetic_mean"
    WEIGHTED_SUM = "weighted_sum"
    WEIGHTED_GEOMETRIC = "weighted_geometric"
    MAX = "max"
    PRODUCT = "product"


class SignalType(Enum):
    """Types of epigenomic signals"""
    ACCESSIBILITY = "accessibility"        # ATAC, DNase
    ACTIVE_ENHANCER = "active_enhancer"   # H3K27ac
    ENHANCER_MARK = "enhancer_mark"       # H3K4me1
    PROMOTER_MARK = "promoter_mark"       # H3K4me3
    TRANSCRIPTION = "transcription"       # H3K36me3
    REPRESSIVE = "repressive"             # H3K27me3, H3K9me3
    TF_BINDING = "tf_binding"             # TF ChIP-seq
    METHYLATION = "methylation"           # DNA methylation (inhibitory)
    EXPRESSION = "expression"             # RNA-seq, eRNA


@dataclass
class SignalConfig:
    """Configuration for a single signal type"""
    name: str
    signal_type: SignalType
    file_path: str
    weight: float = 1.0
    is_inhibitory: bool = False  # True for repressive marks, methylation
    normalize: bool = True
    log_transform: bool = False
    pseudocount: float = 1.0


@dataclass 
class ActivityConfig:
    """Configuration for activity calculation"""
    signals: List[SignalConfig] = field(default_factory=list)
    aggregation_method: AggregationMethod = AggregationMethod.GEOMETRIC_MEAN
    min_signals_required: int = 1
    accessibility_required: bool = True
    output_individual_scores: bool = True


class MultiOmicsActivityCalculator:
    """
    Calculate enhancer activity from multiple epigenomic signals.
    
    This is the core innovation of PACE, allowing flexible integration
    of diverse data types with configurable weights and aggregation methods.
    
    Example usage:
    ```python
    # Basic usage (like original ABC)
    calc = MultiOmicsActivityCalculator()
    calc.add_signal("ATAC", SignalType.ACCESSIBILITY, "atac.bw", weight=1.0)
    calc.add_signal("H3K27ac", SignalType.ACTIVE_ENHANCER, "h3k27ac.bw", weight=1.0)
    activity = calc.calculate_activity(regions)
    
    # Advanced usage (multi-histone + TF + methylation)
    calc = MultiOmicsActivityCalculator(method=AggregationMethod.WEIGHTED_GEOMETRIC)
    calc.add_signal("ATAC", SignalType.ACCESSIBILITY, "atac.bw", weight=1.5)
    calc.add_signal("H3K27ac", SignalType.ACTIVE_ENHANCER, "h3k27ac.bw", weight=1.0)
    calc.add_signal("H3K4me1", SignalType.ENHANCER_MARK, "h3k4me1.bw", weight=0.8)
    calc.add_signal("H3K4me3", SignalType.PROMOTER_MARK, "h3k4me3.bw", weight=0.5)
    calc.add_signal("CTCF", SignalType.TF_BINDING, "ctcf.bw", weight=0.3)
    calc.add_signal("mCG", SignalType.METHYLATION, "methylation.bw", weight=0.5, inhibitory=True)
    activity = calc.calculate_activity(regions)
    ```
    """
    
    def __init__(self, 
                 method: AggregationMethod = AggregationMethod.GEOMETRIC_MEAN,
                 accessibility_required: bool = True):
        """
        Initialize the calculator.
        
        Args:
            method: Aggregation method for combining signals
            accessibility_required: Whether accessibility data is required
        """
        self.signals: List[SignalConfig] = []
        self.method = method
        self.accessibility_required = accessibility_required
        self._signal_data: Dict[str, np.ndarray] = {}
        
    def add_signal(self, 
                   name: str,
                   signal_type: SignalType,
                   file_path: str,
                   weight: float = 1.0,
                   inhibitory: bool = False,
                   normalize: bool = True,
                   log_transform: bool = False) -> None:
        """
        Add a signal to the activity calculation.
        
        Args:
            name: Signal name (e.g., "H3K27ac")
            signal_type: Type of signal
            file_path: Path to signal file (bigWig/BED/etc.)
            weight: Weight for this signal (default 1.0)
            inhibitory: Whether this is an inhibitory signal
            normalize: Whether to normalize the signal
            log_transform: Whether to log-transform the signal
        """
        config = SignalConfig(
            name=name,
            signal_type=signal_type,
            file_path=file_path,
            weight=weight,
            is_inhibitory=inhibitory,
            normalize=normalize,
            log_transform=log_transform
        )
        self.signals.append(config)
        logger.info(f"Added signal: {name} (type={signal_type.value}, weight={weight})")
        
    def _load_signal(self, config: SignalConfig, regions: pd.DataFrame) -> np.ndarray:
        """Load and process signal for given regions"""
        # Implementation depends on file format
        # This is a placeholder - actual implementation would use pyBigWig, etc.
        logger.info(f"Loading signal: {config.name} from {config.file_path}")
        
        # Placeholder: return random values for demonstration
        n_regions = len(regions)
        values = np.random.rand(n_regions) * 100
        
        if config.log_transform:
            values = np.log2(values + config.pseudocount)
            
        if config.normalize:
            values = (values - values.mean()) / (values.std() + 1e-8)
            values = (values - values.min()) / (values.max() - values.min() + 1e-8)
            
        return values
    
    def _aggregate_signals(self, 
                           activating: Dict[str, np.ndarray],
                           inhibitory: Dict[str, np.ndarray],
                           weights_act: Dict[str, float],
                           weights_inh: Dict[str, float]) -> np.ndarray:
        """
        Aggregate multiple signals into a single activity score.
        
        The formula depends on the aggregation method:
        
        GEOMETRIC_MEAN (original ABC):
            Activity = (∏ signal_i)^(1/n)
            
        WEIGHTED_GEOMETRIC:
            Activity = ∏ (signal_i ^ weight_i)
            
        WEIGHTED_SUM:
            Activity = Σ (weight_i × signal_i)
            
        With inhibitory signals:
            Final_Activity = Activating_Score × (1 - Inhibitory_Score)
        """
        n_regions = len(list(activating.values())[0])
        
        # Calculate activating score
        if self.method == AggregationMethod.GEOMETRIC_MEAN:
            # Original ABC method: geometric mean
            product = np.ones(n_regions)
            for name, values in activating.items():
                product *= np.maximum(values, 1e-10)
            activating_score = np.power(product, 1.0 / len(activating))
            
        elif self.method == AggregationMethod.WEIGHTED_GEOMETRIC:
            # Weighted geometric mean
            log_sum = np.zeros(n_regions)
            weight_sum = 0
            for name, values in activating.items():
                w = weights_act[name]
                log_sum += w * np.log(np.maximum(values, 1e-10))
                weight_sum += w
            activating_score = np.exp(log_sum / weight_sum)
            
        elif self.method == AggregationMethod.WEIGHTED_SUM:
            # Weighted sum
            activating_score = np.zeros(n_regions)
            for name, values in activating.items():
                activating_score += weights_act[name] * values
                
        elif self.method == AggregationMethod.ARITHMETIC_MEAN:
            # Simple arithmetic mean
            activating_score = np.mean(list(activating.values()), axis=0)
            
        elif self.method == AggregationMethod.MAX:
            # Maximum value
            activating_score = np.max(list(activating.values()), axis=0)
            
        elif self.method == AggregationMethod.PRODUCT:
            # Product of all signals
            activating_score = np.ones(n_regions)
            for values in activating.values():
                activating_score *= values
        else:
            raise ValueError(f"Unknown aggregation method: {self.method}")
        
        # Apply inhibitory signals
        if inhibitory:
            inhibitory_score = np.zeros(n_regions)
            for name, values in inhibitory.items():
                inhibitory_score += weights_inh[name] * values
            inhibitory_score = np.clip(inhibitory_score, 0, 1)
            
            # Final activity = activating × (1 - inhibitory)
            final_activity = activating_score * (1 - inhibitory_score)
        else:
            final_activity = activating_score
            
        return final_activity
    
    def calculate_activity(self, 
                           regions: pd.DataFrame,
                           return_components: bool = False) -> Union[np.ndarray, Dict]:
        """
        Calculate activity scores for given regions.
        
        Args:
            regions: DataFrame with chr, start, end columns
            return_components: Whether to return individual signal scores
            
        Returns:
            Activity scores (and optionally component scores)
        """
        if self.accessibility_required:
            has_accessibility = any(
                s.signal_type == SignalType.ACCESSIBILITY for s in self.signals
            )
            if not has_accessibility:
                raise ValueError("Accessibility signal is required but not provided")
        
        # Load all signals
        activating = {}
        inhibitory = {}
        weights_act = {}
        weights_inh = {}
        
        for config in self.signals:
            values = self._load_signal(config, regions)
            
            if config.is_inhibitory:
                inhibitory[config.name] = values
                weights_inh[config.name] = config.weight
            else:
                activating[config.name] = values
                weights_act[config.name] = config.weight
        
        # Calculate aggregated activity
        activity = self._aggregate_signals(activating, inhibitory, weights_act, weights_inh)
        
        if return_components:
            return {
                'activity': activity,
                'activating_signals': activating,
                'inhibitory_signals': inhibitory,
                'weights': {'activating': weights_act, 'inhibitory': weights_inh}
            }
        else:
            return activity


class ExpressionFilter:
    """
    Filter enhancer-gene predictions based on gene expression.
    
    Only genes with sufficient expression are considered as potential
    targets, as unexpressed genes cannot be regulated.
    """
    
    def __init__(self, 
                 expression_file: str,
                 min_expression: float = 1.0,
                 expression_column: str = "TPM"):
        """
        Args:
            expression_file: Path to expression file (TSV/CSV)
            min_expression: Minimum expression threshold
            expression_column: Column name for expression values
        """
        self.expression_file = expression_file
        self.min_expression = min_expression
        self.expression_column = expression_column
        self._expression_data = None
        
    def load_expression(self) -> pd.DataFrame:
        """Load expression data"""
        logger.info(f"Loading expression data from {self.expression_file}")
        self._expression_data = pd.read_csv(self.expression_file, sep='\t')
        return self._expression_data
    
    def get_expressed_genes(self) -> List[str]:
        """Get list of expressed genes"""
        if self._expression_data is None:
            self.load_expression()
            
        expressed = self._expression_data[
            self._expression_data[self.expression_column] >= self.min_expression
        ]
        return expressed['gene_id'].tolist()
    
    def filter_predictions(self, predictions: pd.DataFrame) -> pd.DataFrame:
        """Filter predictions to only include expressed genes"""
        expressed_genes = set(self.get_expressed_genes())
        
        filtered = predictions[
            predictions['TargetGene'].isin(expressed_genes)
        ]
        
        logger.info(f"Filtered {len(predictions)} -> {len(filtered)} predictions "
                   f"({len(expressed_genes)} expressed genes)")
        return filtered
    
    def add_expression_weight(self, 
                              predictions: pd.DataFrame,
                              weight_method: str = "log") -> pd.DataFrame:
        """
        Add expression-based weight to predictions.
        
        Args:
            predictions: Prediction DataFrame
            weight_method: "log", "linear", or "binary"
        """
        if self._expression_data is None:
            self.load_expression()
            
        expr_dict = dict(zip(
            self._expression_data['gene_id'],
            self._expression_data[self.expression_column]
        ))
        
        predictions = predictions.copy()
        predictions['expression'] = predictions['TargetGene'].map(expr_dict).fillna(0)
        
        if weight_method == "log":
            predictions['expression_weight'] = np.log2(predictions['expression'] + 1)
        elif weight_method == "linear":
            max_expr = predictions['expression'].max()
            predictions['expression_weight'] = predictions['expression'] / max_expr
        elif weight_method == "binary":
            predictions['expression_weight'] = (predictions['expression'] >= self.min_expression).astype(float)
        else:
            raise ValueError(f"Unknown weight method: {weight_method}")
            
        return predictions


class eQTLValidator:
    """
    Validate enhancer-gene predictions using eQTL data.
    
    eQTL (expression Quantitative Trait Loci) provide genetic evidence
    for regulatory relationships between genomic variants and gene expression.
    """
    
    def __init__(self, eqtl_file: str):
        """
        Args:
            eqtl_file: Path to eQTL summary statistics file
        """
        self.eqtl_file = eqtl_file
        self._eqtl_data = None
        
    def load_eqtl(self) -> pd.DataFrame:
        """Load eQTL data"""
        logger.info(f"Loading eQTL data from {self.eqtl_file}")
        self._eqtl_data = pd.read_csv(self.eqtl_file, sep='\t')
        return self._eqtl_data
    
    def validate_predictions(self, 
                             predictions: pd.DataFrame,
                             pvalue_threshold: float = 1e-5,
                             window: int = 1000) -> pd.DataFrame:
        """
        Add eQTL validation evidence to predictions.
        
        Args:
            predictions: Prediction DataFrame with enhancer coordinates
            pvalue_threshold: P-value threshold for significant eQTLs
            window: Window size for matching enhancers to eQTL variants
            
        Returns:
            Predictions with eQTL evidence columns added
        """
        if self._eqtl_data is None:
            self.load_eqtl()
            
        predictions = predictions.copy()
        
        # Find eQTLs overlapping with enhancers
        # This is a simplified implementation
        predictions['eQTL_support'] = False
        predictions['eQTL_pvalue'] = np.nan
        predictions['eQTL_beta'] = np.nan
        
        # Implementation would match eQTL variants to enhancer regions
        # and check if the eQTL target gene matches the predicted target
        
        logger.info("eQTL validation completed")
        return predictions
    
    def calculate_enrichment(self, predictions: pd.DataFrame) -> Dict:
        """
        Calculate enrichment of eQTLs in predicted enhancers.
        
        Returns statistics on how well predictions are supported by eQTL evidence.
        """
        # Calculate fold enrichment of eQTLs in high-score enhancers
        # vs random genomic regions
        
        enrichment_stats = {
            'n_predictions': len(predictions),
            'n_with_eqtl': predictions['eQTL_support'].sum(),
            'fraction_supported': predictions['eQTL_support'].mean(),
            'fold_enrichment': np.nan,  # Would need background rate
            'pvalue': np.nan
        }
        
        return enrichment_stats


class TFEnrichmentAnalyzer:
    """
    Analyze transcription factor binding enrichment at enhancers.
    
    Identifies which TFs are enriched at predicted enhancers,
    providing insight into regulatory mechanisms.
    """
    
    def __init__(self, tf_binding_dir: str = None, motif_file: str = None):
        """
        Args:
            tf_binding_dir: Directory containing TF ChIP-seq peak files
            motif_file: File containing TF motif matches
        """
        self.tf_binding_dir = tf_binding_dir
        self.motif_file = motif_file
        
    def calculate_tf_activity(self, 
                              regions: pd.DataFrame,
                              tf_list: Optional[List[str]] = None) -> pd.DataFrame:
        """
        Calculate TF binding activity at enhancer regions.
        
        Args:
            regions: DataFrame with enhancer coordinates
            tf_list: List of TFs to analyze (None = all available)
            
        Returns:
            DataFrame with TF binding scores for each region
        """
        logger.info("Calculating TF binding activity")
        
        # Implementation would:
        # 1. Load TF ChIP-seq peaks
        # 2. Count overlaps with enhancer regions
        # 3. Calculate binding scores
        
        tf_scores = regions.copy()
        # Add TF binding columns...
        
        return tf_scores
    
    def find_enriched_tfs(self, 
                          enhancers: pd.DataFrame,
                          background: pd.DataFrame,
                          min_fold_change: float = 2.0) -> pd.DataFrame:
        """
        Find TFs enriched at enhancers compared to background.
        
        Returns DataFrame with TF enrichment statistics.
        """
        # Fisher's exact test or hypergeometric test
        # for each TF
        
        enrichment_results = []
        # Implementation...
        
        return pd.DataFrame(enrichment_results)


class MethylationIntegrator:
    """
    Integrate DNA methylation data as regulatory inhibitor.
    
    High DNA methylation at enhancers typically indicates
    repressed/inactive state.
    """
    
    def __init__(self, 
                 methylation_file: str,
                 context: str = "CG"):
        """
        Args:
            methylation_file: Path to methylation file (bedGraph/bigWig)
            context: Methylation context (CG, CHG, CHH)
        """
        self.methylation_file = methylation_file
        self.context = context
        
    def get_methylation_scores(self, regions: pd.DataFrame) -> np.ndarray:
        """
        Get methylation levels for given regions.
        
        Returns values between 0-1 representing methylation fraction.
        """
        logger.info(f"Loading methylation data from {self.methylation_file}")
        
        # Implementation would load methylation data
        # and calculate average methylation per region
        
        n_regions = len(regions)
        methylation = np.random.rand(n_regions)  # Placeholder
        
        return methylation
    
    def calculate_inhibitory_factor(self, 
                                    methylation: np.ndarray,
                                    method: str = "linear") -> np.ndarray:
        """
        Convert methylation to inhibitory factor.
        
        Args:
            methylation: Methylation levels (0-1)
            method: "linear" or "sigmoid"
            
        Returns:
            Inhibitory factor (0-1, higher = more inhibition)
        """
        if method == "linear":
            return methylation
        elif method == "sigmoid":
            # Sigmoid transformation for non-linear effect
            return 1 / (1 + np.exp(-10 * (methylation - 0.5)))
        else:
            raise ValueError(f"Unknown method: {method}")


# Convenience function for typical usage
def calculate_pace_score(
    regions: pd.DataFrame,
    atac_file: str,
    h3k27ac_file: str = None,
    h3k4me1_file: str = None,
    h3k4me3_file: str = None,
    tf_files: Dict[str, str] = None,
    methylation_file: str = None,
    method: AggregationMethod = AggregationMethod.GEOMETRIC_MEAN,
    weights: Dict[str, float] = None
) -> np.ndarray:
    """
    Convenience function to calculate PACE activity score.
    
    Args:
        regions: DataFrame with enhancer regions
        atac_file: Path to ATAC-seq signal file (required)
        h3k27ac_file: Path to H3K27ac signal file (optional)
        h3k4me1_file: Path to H3K4me1 signal file (optional)
        h3k4me3_file: Path to H3K4me3 signal file (optional)
        tf_files: Dictionary of TF name -> file path
        methylation_file: Path to methylation file (optional)
        method: Aggregation method
        weights: Custom weights for each signal
        
    Returns:
        Activity scores for each region
    """
    weights = weights or {}
    
    calc = MultiOmicsActivityCalculator(method=method)
    
    # Add accessibility (required)
    calc.add_signal("ATAC", SignalType.ACCESSIBILITY, atac_file, 
                    weight=weights.get("ATAC", 1.0))
    
    # Add histone modifications (optional)
    if h3k27ac_file:
        calc.add_signal("H3K27ac", SignalType.ACTIVE_ENHANCER, h3k27ac_file,
                        weight=weights.get("H3K27ac", 1.0))
    
    if h3k4me1_file:
        calc.add_signal("H3K4me1", SignalType.ENHANCER_MARK, h3k4me1_file,
                        weight=weights.get("H3K4me1", 0.8))
    
    if h3k4me3_file:
        calc.add_signal("H3K4me3", SignalType.PROMOTER_MARK, h3k4me3_file,
                        weight=weights.get("H3K4me3", 0.5))
    
    # Add TF binding (optional)
    if tf_files:
        for tf_name, tf_file in tf_files.items():
            calc.add_signal(tf_name, SignalType.TF_BINDING, tf_file,
                            weight=weights.get(tf_name, 0.3))
    
    # Add methylation as inhibitory signal (optional)
    if methylation_file:
        calc.add_signal("methylation", SignalType.METHYLATION, methylation_file,
                        weight=weights.get("methylation", 0.5), inhibitory=True)
    
    return calc.calculate_activity(regions)


if __name__ == "__main__":
    # Example usage
    print("PACE Multi-omics Activity Calculator")
    print("=" * 50)
    
    # Create example regions
    regions = pd.DataFrame({
        'chr': ['chr1'] * 10,
        'start': range(1000, 11000, 1000),
        'end': range(1500, 11500, 1000)
    })
    
    # Example 1: Basic (like original ABC)
    print("\nExample 1: Basic activity (ATAC + H3K27ac)")
    calc = MultiOmicsActivityCalculator()
    calc.add_signal("ATAC", SignalType.ACCESSIBILITY, "atac.bw")
    calc.add_signal("H3K27ac", SignalType.ACTIVE_ENHANCER, "h3k27ac.bw")
    activity1 = calc.calculate_activity(regions)
    print(f"Activity scores: {activity1[:5]}")
    
    # Example 2: Multi-histone
    print("\nExample 2: Multi-histone activity")
    calc2 = MultiOmicsActivityCalculator(method=AggregationMethod.WEIGHTED_GEOMETRIC)
    calc2.add_signal("ATAC", SignalType.ACCESSIBILITY, "atac.bw", weight=1.5)
    calc2.add_signal("H3K27ac", SignalType.ACTIVE_ENHANCER, "h3k27ac.bw", weight=1.0)
    calc2.add_signal("H3K4me1", SignalType.ENHANCER_MARK, "h3k4me1.bw", weight=0.8)
    calc2.add_signal("H3K4me3", SignalType.PROMOTER_MARK, "h3k4me3.bw", weight=0.5)
    activity2 = calc2.calculate_activity(regions)
    print(f"Activity scores: {activity2[:5]}")
    
    # Example 3: With methylation inhibition
    print("\nExample 3: With methylation inhibition")
    calc3 = MultiOmicsActivityCalculator()
    calc3.add_signal("ATAC", SignalType.ACCESSIBILITY, "atac.bw")
    calc3.add_signal("H3K27ac", SignalType.ACTIVE_ENHANCER, "h3k27ac.bw")
    calc3.add_signal("mCG", SignalType.METHYLATION, "methylation.bw", 
                     weight=0.5, inhibitory=True)
    activity3 = calc3.calculate_activity(regions)
    print(f"Activity scores: {activity3[:5]}")
    
    print("\nDone!")
