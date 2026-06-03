#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PACE Neighborhoods Module

Functions for neighborhood analysis - quantifying chromatin features around enhancers and genes.

Author: Linyong Shen @ Northwest A&F University
"""

import os
import numpy as np
import pandas as pd
from typing import Optional, List, Dict, Tuple, Union
import tempfile

# Optional pyBigWig import
try:
    import pyBigWig
    HAS_PYBIGWIG = True
except ImportError:
    pyBigWig = None
    HAS_PYBIGWIG = False

from tools import (
    read_bed, write_bed, run_command, bedtools_intersect,
    normalize_signal, geometric_mean, weighted_geometric_mean,
    logger
)


class NeighborhoodAnalyzer:
    """
    Analyzer for quantifying chromatin features in enhancer neighborhoods.
    """
    
    def __init__(self,
                 candidate_regions: pd.DataFrame,
                 genes: pd.DataFrame,
                 chrom_sizes: Dict[str, int]):
        """
        Initialize analyzer.
        
        Args:
            candidate_regions: DataFrame with candidate enhancer regions
            genes: DataFrame with gene information
            chrom_sizes: Chromosome sizes dictionary
        """
        self.regions = candidate_regions.copy()
        self.genes = genes.copy()
        self.chrom_sizes = chrom_sizes
        
        # Ensure required columns
        if 'name' not in self.regions.columns:
            self.regions['name'] = [f"region_{i}" for i in range(len(self.regions))]
        
        # Initialize signal columns
        self.signal_columns = []
    
    def quantify_signal_from_bigwig(self,
                                    bigwig_file: str,
                                    column_name: str,
                                    stat: str = 'mean') -> None:
        """
        Quantify signal from bigWig file over regions.
        
        Args:
            bigwig_file: Path to bigWig file
            column_name: Name for output column
            stat: Statistic to compute (mean, max, sum)
        """
        if not HAS_PYBIGWIG:
            raise ImportError("pyBigWig is required for bigWig processing. "
                            "Install with: pip install pyBigWig")
        
        logger.info(f"Quantifying {column_name} from {bigwig_file}")
        
        bw = pyBigWig.open(bigwig_file)
        values = []
        
        for _, row in self.regions.iterrows():
            chrom = row['chr']
            start = int(row['start'])
            end = int(row['end'])
            
            try:
                if stat == 'mean':
                    val = bw.stats(chrom, start, end, type='mean')[0]
                elif stat == 'max':
                    val = bw.stats(chrom, start, end, type='max')[0]
                elif stat == 'sum':
                    val = bw.stats(chrom, start, end, type='sum')[0]
                else:
                    val = bw.stats(chrom, start, end, type='mean')[0]
                
                values.append(val if val is not None else 0)
            except Exception:
                values.append(0)
        
        bw.close()
        
        self.regions[column_name] = values
        self.signal_columns.append(column_name)
        
        logger.info(f"Added signal column: {column_name}")
    
    def quantify_signal_from_bam(self,
                                 bam_file: str,
                                 column_name: str,
                                 normalize: bool = True) -> None:
        """
        Quantify signal from BAM file over regions.
        
        Args:
            bam_file: Path to BAM file
            column_name: Name for output column
            normalize: Whether to normalize to RPM
        """
        logger.info(f"Quantifying {column_name} from {bam_file}")
        
        # Write regions to temp file
        with tempfile.NamedTemporaryFile(suffix='.bed', delete=False, mode='w') as tmp:
            self.regions[['chr', 'start', 'end', 'name']].to_csv(
                tmp, sep='\t', header=False, index=False
            )
            tmp_name = tmp.name
        
        # Count reads using bedtools
        with tempfile.NamedTemporaryFile(suffix='.bed', delete=False) as tmp_out:
            cmd = f"bedtools coverage -a {tmp_name} -b {bam_file} -counts > {tmp_out.name}"
            run_command(cmd)
            
            counts = pd.read_csv(
                tmp_out.name, sep='\t', header=None,
                names=['chr', 'start', 'end', 'name', 'count']
            )
            os.unlink(tmp_out.name)
        
        os.unlink(tmp_name)
        
        # Merge counts
        count_dict = dict(zip(counts['name'], counts['count']))
        self.regions[column_name] = self.regions['name'].map(count_dict).fillna(0)
        
        # Normalize if requested
        if normalize:
            total = self.regions[column_name].sum()
            if total > 0:
                self.regions[column_name] = self.regions[column_name] * 1e6 / total
        
        self.signal_columns.append(column_name)
        logger.info(f"Added signal column: {column_name}")
    
    def quantify_signal_from_tagalign(self,
                                      tagalign_file: str,
                                      column_name: str,
                                      normalize: bool = True) -> None:
        """
        Quantify signal from tagAlign file over regions.
        
        Args:
            tagalign_file: Path to tagAlign file
            column_name: Name for output column
            normalize: Whether to normalize to RPM
        """
        # tagAlign is essentially BED format, use same method as BAM
        self.quantify_signal_from_bam(tagalign_file, column_name, normalize)
    
    def quantify_methylation(self,
                            methylation_file: str,
                            column_name: str = 'methylation',
                            context: str = 'CG') -> None:
        """
        Quantify DNA methylation over regions.
        
        Args:
            methylation_file: Path to methylation file (bigWig or bedGraph)
            column_name: Name for output column
            context: Methylation context (CG, CHG, CHH)
        """
        if methylation_file.endswith('.bw') or methylation_file.endswith('.bigwig'):
            self.quantify_signal_from_bigwig(methylation_file, column_name)
        else:
            # bedGraph format
            logger.info(f"Quantifying methylation from {methylation_file}")
            
            # Load methylation data
            meth = pd.read_csv(
                methylation_file, sep='\t', header=None,
                names=['chr', 'start', 'end', 'value']
            )
            
            # Write regions to temp file
            with tempfile.NamedTemporaryFile(suffix='.bed', delete=False, mode='w') as tmp:
                self.regions[['chr', 'start', 'end', 'name']].to_csv(
                    tmp, sep='\t', header=False, index=False
                )
                tmp_name = tmp.name
            
            # Write methylation to temp file
            with tempfile.NamedTemporaryFile(suffix='.bed', delete=False, mode='w') as tmp_meth:
                meth.to_csv(tmp_meth, sep='\t', header=False, index=False)
                tmp_meth_name = tmp_meth.name
            
            # Use bedtools map to get mean methylation
            with tempfile.NamedTemporaryFile(suffix='.bed', delete=False) as tmp_out:
                cmd = f"bedtools map -a {tmp_name} -b {tmp_meth_name} -c 4 -o mean > {tmp_out.name}"
                run_command(cmd)
                
                result = pd.read_csv(
                    tmp_out.name, sep='\t', header=None,
                    names=['chr', 'start', 'end', 'name', column_name]
                )
                os.unlink(tmp_out.name)
            
            os.unlink(tmp_name)
            os.unlink(tmp_meth_name)
            
            # Handle '.' values (no data)
            result[column_name] = pd.to_numeric(result[column_name], errors='coerce').fillna(0)
            
            # Merge back
            meth_dict = dict(zip(result['name'], result[column_name]))
            self.regions[column_name] = self.regions['name'].map(meth_dict).fillna(0)
            
            self.signal_columns.append(column_name)
    
    def calculate_activity(self,
                          accessibility_column: str,
                          histone_columns: Optional[List[str]] = None,
                          histone_weights: Optional[Dict[str, float]] = None,
                          inhibitory_columns: Optional[List[str]] = None,
                          inhibitory_weights: Optional[Dict[str, float]] = None,
                          method: str = 'geometric_mean',
                          accessibility_weight: float = 1.5) -> None:
        """
        Calculate enhancer activity score.

        Activating signals are aggregated according to ``method`` and then
        scaled down by an inhibitory factor derived from repressive marks /
        DNA methylation:

            geometric_mean:      A = (prod s_i) ** (1/n)              [ABC]
            weighted_geometric:  A = (prod s_i ** w_i) ** (1/sum w_i) [PACE]
            weighted_sum:        A = sum (w_i * s_i)
            A_final = A * (1 - clip(sum_k w_k * minmax(I_k), 0, 1))

        The weighted geometric mean is normalized by the sum of weights so the
        activity scale is independent of the absolute weight magnitudes
        (consistent with scripts/multiomics_activity.py).

        Args:
            accessibility_column: Column with accessibility signal.
            histone_columns: Columns with activating histone signals.
            histone_weights: Weights keyed by column name.
            inhibitory_columns: Columns with inhibitory signals
                (methylation, H3K27me3, H3K9me3).
            inhibitory_weights: Weights for inhibitory signals.
            method: Aggregation method.
            accessibility_weight: Weight for the accessibility signal.
        """
        logger.info(f"Calculating activity using {method}")

        # Prepare activating signals and their weights.
        signal_arrays = [self.regions[accessibility_column].to_numpy(dtype=float)]
        weights = [accessibility_weight]

        if histone_columns:
            for col in histone_columns:
                if col in self.regions.columns:
                    signal_arrays.append(self.regions[col].to_numpy(dtype=float))
                    w = histone_weights.get(col, 1.0) if histone_weights else 1.0
                    weights.append(w)

        # Stack without mutating the source arrays.
        activating = np.clip(np.vstack(signal_arrays), 1e-10, None)
        weights = np.asarray(weights, dtype=float)

        if method == 'geometric_mean':
            activity = np.exp(np.mean(np.log(activating), axis=0))

        elif method == 'weighted_geometric':
            w_sum = weights.sum() if weights.sum() > 0 else 1.0
            log_activity = np.tensordot(weights, np.log(activating), axes=1) / w_sum
            activity = np.exp(log_activity)

        elif method == 'weighted_sum':
            activity = np.tensordot(weights, activating, axes=1)

        else:  # arithmetic_mean
            activity = np.mean(activating, axis=0)

        # Apply inhibitory signals (each scaled to [0, 1] independently).
        if inhibitory_columns:
            inhibitory_score = np.zeros(len(self.regions))
            for col in inhibitory_columns:
                if col in self.regions.columns:
                    w = (inhibitory_weights.get(col, 0.5)
                         if inhibitory_weights else 0.5)
                    inhib = self.regions[col].to_numpy(dtype=float)
                    rng = inhib.max() - inhib.min()
                    inhib_norm = ((inhib - inhib.min()) / rng
                                  if rng > 0 else np.zeros_like(inhib))
                    inhibitory_score += w * inhib_norm
            inhibitory_score = np.clip(inhibitory_score, 0, 1)
            activity = activity * (1 - inhibitory_score)

        self.regions['activity'] = activity
        logger.info(f"Activity range: {activity.min():.4f} - {activity.max():.4f}")
    
    def normalize_signals(self,
                         columns: Optional[List[str]] = None,
                         method: str = 'quantile',
                         reference_file: Optional[str] = None) -> None:
        """
        Normalize signal columns.
        
        Args:
            columns: Columns to normalize (default: all signal columns)
            method: Normalization method
            reference_file: Reference for quantile normalization
        """
        if columns is None:
            columns = self.signal_columns
        
        for col in columns:
            if col in self.regions.columns:
                self.regions[col] = normalize_signal(
                    self.regions[col].values,
                    method=method
                )
        
        logger.info(f"Normalized {len(columns)} signal columns")
    
    def create_enhancer_list(self,
                            output_file: str,
                            include_activity: bool = True) -> pd.DataFrame:
        """
        Create enhancer list file.
        
        Args:
            output_file: Output file path
            include_activity: Whether to include activity score
        
        Returns:
            DataFrame with enhancer list
        """
        output_cols = ['chr', 'start', 'end', 'name']
        
        for col in self.signal_columns:
            if col in self.regions.columns:
                output_cols.append(col)
        
        if include_activity and 'activity' in self.regions.columns:
            output_cols.append('activity')
        
        enhancers = self.regions[output_cols].copy()
        enhancers.to_csv(output_file, sep='\t', index=False)
        
        logger.info(f"Wrote {len(enhancers)} enhancers to {output_file}")
        
        return enhancers
    
    def create_gene_list(self,
                        output_file: str,
                        expression_file: Optional[str] = None,
                        expression_column: str = 'TPM') -> pd.DataFrame:
        """
        Create gene list file with expression data.
        
        Args:
            output_file: Output file path
            expression_file: Optional RNA-seq expression file
            expression_column: Column name in expression file
        
        Returns:
            DataFrame with gene list
        """
        genes = self.genes.copy()
        
        # Add expression if available
        if expression_file:
            expr = pd.read_csv(expression_file, sep='\t')
            
            # Try different gene ID columns
            for id_col in ['gene_id', 'gene_name', 'GeneID', 'Gene']:
                if id_col in expr.columns:
                    gene_id_col = id_col
                    break
            else:
                gene_id_col = expr.columns[0]
            
            # Create expression dictionary
            if expression_column in expr.columns:
                expr_dict = dict(zip(expr[gene_id_col], expr[expression_column]))
            else:
                # Use first numeric column
                numeric_cols = expr.select_dtypes(include=[np.number]).columns
                if len(numeric_cols) > 0:
                    expr_dict = dict(zip(expr[gene_id_col], expr[numeric_cols[0]]))
                else:
                    expr_dict = {}
            
            # Map expression to genes
            for id_col in ['gene_name', 'gene_id', 'name']:
                if id_col in genes.columns:
                    genes['Expression'] = genes[id_col].map(expr_dict).fillna(0)
                    break
            else:
                genes['Expression'] = 0
        
        genes.to_csv(output_file, sep='\t', index=False)
        logger.info(f"Wrote {len(genes)} genes to {output_file}")
        
        return genes


def run_neighborhoods(candidate_regions_file: str,
                     genes_file: str,
                     chrom_sizes_file: str,
                     output_dir: str,
                     accessibility_file: str,
                     accessibility_type: str = 'ATAC',
                     histone_files: Optional[Dict[str, str]] = None,
                     methylation_file: Optional[str] = None,
                     expression_file: Optional[str] = None,
                     activity_method: str = 'geometric_mean',
                     histone_weights: Optional[Dict[str, float]] = None,
                     inhibitory_histone_files: Optional[Dict[str, str]] = None,
                     inhibitory_weights: Optional[Dict[str, float]] = None,
                     accessibility_weight: float = 1.5) -> Tuple[str, str]:
    """
    Run neighborhood analysis.

    This is the main entry point for neighborhood analysis.

    Args:
        candidate_regions_file: Path to candidate regions BED
        genes_file: Path to genes BED
        chrom_sizes_file: Path to chromosome sizes
        output_dir: Output directory
        accessibility_file: ATAC or DNase file
        accessibility_type: Type of accessibility data
        histone_files: Dictionary of activating histone mark files
        methylation_file: Optional methylation file (inhibitory)
        expression_file: Optional RNA-seq expression file
        activity_method: Activity calculation method
        histone_weights: Weights for histone marks (keyed by '<mark>_signal')
        inhibitory_histone_files: Repressive marks (H3K27me3, H3K9me3) used as
            inhibitory signals.
        inhibitory_weights: Weights for inhibitory signals (keyed by
            '<mark>_signal' / 'methylation_signal').
        accessibility_weight: Weight for the accessibility signal.

    Returns:
        Tuple of (enhancer_list_file, gene_list_file)
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Load data
    logger.info("Loading candidate regions and genes")
    
    regions = pd.read_csv(candidate_regions_file, sep='\t', header=None,
                         names=['chr', 'start', 'end', 'name'])
    
    # Read the gene annotation robustly. PACE gene BEDs may be BED6 or an
    # extended BED with gene_id / gene_type columns (as produced by
    # prepare_reference.py). Passing too few column names to pandas silently
    # corrupts the table (extra leading columns become an index), so we size
    # the names to the actual column count.
    with open(genes_file) as f:
        first_line = f.readline().rstrip('\n').split('\t')

    has_header = first_line[0].lower() in [
        'chr', 'chrom', 'chromosome', '#chr', '#chrom', 'gene_id', 'geneid']

    if has_header:
        genes = pd.read_csv(genes_file, sep='\t')
    else:
        ncol = len(first_line)
        base = ['chr', 'start', 'end', 'name', 'score', 'strand',
                'gene_id', 'gene_type']
        if ncol <= len(base):
            col_names = base[:ncol]
        else:
            col_names = base + [f'extra_{i}' for i in range(len(base), ncol)]
        genes = pd.read_csv(genes_file, sep='\t', header=None, names=col_names)

    # Ensure 'name' column exists for gene identification.
    if 'name' not in genes.columns and 'gene_name' not in genes.columns:
        if 'symbol' in genes.columns:
            genes['name'] = genes['symbol']
        elif len(genes.columns) >= 4:
            genes['name'] = genes.iloc[:, 3]

    # The name field may encode 'GENE_SYMBOL;ENSEMBL_ID'; derive a clean
    # gene symbol and keep the Ensembl id where available.
    if 'name' in genes.columns:
        name_str = genes['name'].astype(str)
        genes['gene_name'] = name_str.str.split(';').str[0]
        if 'gene_id' not in genes.columns:
            genes['gene_id'] = name_str.str.split(';').str[-1]

    # Strand-aware TSS (5'-most position of the gene).
    if 'TSS' not in genes.columns and {'start', 'end', 'strand'}.issubset(genes.columns):
        genes['TSS'] = genes.apply(
            lambda r: int(r['start']) if str(r['strand']) != '-' else int(r['end']),
            axis=1)
    
    chrom_sizes = {}
    with open(chrom_sizes_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                chrom_sizes[parts[0]] = int(parts[1])
    
    # Initialize analyzer
    analyzer = NeighborhoodAnalyzer(regions, genes, chrom_sizes)
    
    # Quantify accessibility
    acc_col = f'{accessibility_type}_signal'
    if accessibility_file.endswith('.bw') or accessibility_file.endswith('.bigwig'):
        analyzer.quantify_signal_from_bigwig(accessibility_file, acc_col)
    else:
        analyzer.quantify_signal_from_tagalign(accessibility_file, acc_col)
    
    # Quantify histone marks
    histone_cols = []
    if histone_files:
        for mark, filepath in histone_files.items():
            if filepath and os.path.exists(filepath):
                col_name = f'{mark}_signal'
                if filepath.endswith('.bw') or filepath.endswith('.bigwig'):
                    analyzer.quantify_signal_from_bigwig(filepath, col_name)
                else:
                    analyzer.quantify_signal_from_bam(filepath, col_name)
                histone_cols.append(col_name)
    
    # Quantify inhibitory signals: repressive histone marks + methylation.
    inhibitory_cols = []
    if inhibitory_histone_files:
        for mark, filepath in inhibitory_histone_files.items():
            if filepath and os.path.exists(filepath):
                col_name = f'{mark}_signal'
                if filepath.endswith('.bw') or filepath.endswith('.bigwig'):
                    analyzer.quantify_signal_from_bigwig(filepath, col_name)
                else:
                    analyzer.quantify_signal_from_bam(filepath, col_name)
                inhibitory_cols.append(col_name)

    # Quantify methylation
    if methylation_file and os.path.exists(methylation_file):
        analyzer.quantify_methylation(methylation_file, 'methylation_signal')
        inhibitory_cols.append('methylation_signal')

    # Calculate activity
    analyzer.calculate_activity(
        accessibility_column=acc_col,
        histone_columns=histone_cols if histone_cols else None,
        histone_weights=histone_weights,
        inhibitory_columns=inhibitory_cols if inhibitory_cols else None,
        inhibitory_weights=inhibitory_weights,
        method=activity_method,
        accessibility_weight=accessibility_weight,
    )
    
    # Write output files
    enhancer_file = os.path.join(output_dir, 'EnhancerList.txt')
    gene_file = os.path.join(output_dir, 'GeneList.txt')
    
    analyzer.create_enhancer_list(enhancer_file)
    analyzer.create_gene_list(gene_file, expression_file)
    
    return enhancer_file, gene_file
