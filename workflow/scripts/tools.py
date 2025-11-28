#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PACE Utility Functions

Core utility functions for the PACE pipeline.
Based on ABC-Enhancer-Gene-Prediction with enhancements for livestock genomics.

Author: Linyong Shen @ Northwest A&F University
"""

import os
import sys
import gzip
import numpy as np
import pandas as pd
from typing import Optional, List, Dict, Union, Tuple
import subprocess
import tempfile


def read_bed(filename: str, 
             names: Optional[List[str]] = None,
             usecols: Optional[List[int]] = None,
             **kwargs) -> pd.DataFrame:
    """
    Read BED file into DataFrame.
    
    Args:
        filename: Path to BED file (can be gzipped)
        names: Column names
        usecols: Columns to read
        **kwargs: Additional pandas read_csv arguments
    
    Returns:
        DataFrame with BED data
    """
    if names is None:
        names = ['chr', 'start', 'end', 'name', 'score', 'strand']
    
    open_func = gzip.open if filename.endswith('.gz') else open
    
    return pd.read_csv(
        filename,
        sep='\t',
        header=None,
        names=names[:len(usecols)] if usecols else names,
        usecols=usecols,
        **kwargs
    )


def write_bed(df: pd.DataFrame, 
              filename: str,
              columns: Optional[List[str]] = None) -> None:
    """
    Write DataFrame to BED file.
    
    Args:
        df: DataFrame to write
        filename: Output filename
        columns: Columns to write (default: all)
    """
    if columns:
        df = df[columns]
    
    if filename.endswith('.gz'):
        df.to_csv(filename, sep='\t', header=False, index=False, compression='gzip')
    else:
        df.to_csv(filename, sep='\t', header=False, index=False)


def run_command(cmd: Union[str, List[str]], 
                shell: bool = True,
                check: bool = True,
                capture_output: bool = False) -> subprocess.CompletedProcess:
    """
    Run shell command.
    
    Args:
        cmd: Command string or list
        shell: Use shell execution
        check: Raise on non-zero exit
        capture_output: Capture stdout/stderr
    
    Returns:
        CompletedProcess instance
    """
    return subprocess.run(
        cmd,
        shell=shell,
        check=check,
        capture_output=capture_output,
        text=True
    )


def bedtools_intersect(a_file: str, 
                       b_file: str,
                       output: Optional[str] = None,
                       options: str = "-wa -wb") -> pd.DataFrame:
    """
    Run bedtools intersect.
    
    Args:
        a_file: First BED file
        b_file: Second BED file
        output: Output file (if None, return DataFrame)
        options: bedtools intersect options
    
    Returns:
        DataFrame if output is None
    """
    if output:
        cmd = f"bedtools intersect -a {a_file} -b {b_file} {options} > {output}"
        run_command(cmd)
        return None
    else:
        with tempfile.NamedTemporaryFile(suffix='.bed', delete=False) as tmp:
            cmd = f"bedtools intersect -a {a_file} -b {b_file} {options} > {tmp.name}"
            run_command(cmd)
            df = pd.read_csv(tmp.name, sep='\t', header=None)
            os.unlink(tmp.name)
            return df


def bedtools_coverage(a_file: str,
                      b_file: str,
                      output: Optional[str] = None) -> pd.DataFrame:
    """
    Run bedtools coverage.
    
    Args:
        a_file: Regions file
        b_file: Signal file
        output: Output file
    
    Returns:
        DataFrame with coverage
    """
    if output:
        cmd = f"bedtools coverage -a {a_file} -b {b_file} > {output}"
        run_command(cmd)
        return None
    else:
        with tempfile.NamedTemporaryFile(suffix='.bed', delete=False) as tmp:
            cmd = f"bedtools coverage -a {a_file} -b {b_file} > {tmp.name}"
            run_command(cmd)
            df = pd.read_csv(tmp.name, sep='\t', header=None)
            os.unlink(tmp.name)
            return df


def normalize_signal(signal: np.ndarray, 
                     method: str = 'quantile',
                     reference: Optional[np.ndarray] = None) -> np.ndarray:
    """
    Normalize signal values.
    
    Args:
        signal: Signal array
        method: Normalization method (quantile, zscore, minmax, rpm)
        reference: Reference for quantile normalization
    
    Returns:
        Normalized signal
    """
    if method == 'quantile' and reference is not None:
        # Quantile normalization
        sorted_idx = np.argsort(signal)
        sorted_ref = np.sort(reference)
        normalized = np.zeros_like(signal)
        normalized[sorted_idx] = sorted_ref[:len(signal)]
        return normalized
    
    elif method == 'zscore':
        mean = np.mean(signal)
        std = np.std(signal)
        return (signal - mean) / std if std > 0 else signal - mean
    
    elif method == 'minmax':
        min_val = np.min(signal)
        max_val = np.max(signal)
        return (signal - min_val) / (max_val - min_val) if max_val > min_val else signal
    
    elif method == 'rpm':
        total = np.sum(signal)
        return signal * 1e6 / total if total > 0 else signal
    
    else:
        return signal


def calculate_distance(enhancer_mid: int, tss: int) -> int:
    """
    Calculate distance from enhancer midpoint to TSS.
    
    Args:
        enhancer_mid: Enhancer midpoint
        tss: TSS position
    
    Returns:
        Distance (absolute value)
    """
    return abs(enhancer_mid - tss)


def power_law_contact(distance: Union[int, np.ndarray],
                      hic_gamma: float = 1.024238616787792,
                      hic_scale: float = 5.9594510043736655,
                      hic_pseudocount: int = 5000) -> Union[float, np.ndarray]:
    """
    Estimate contact frequency using power-law.
    
    Args:
        distance: Distance(s) in bp
        hic_gamma: Power-law exponent
        hic_scale: Scale factor
        hic_pseudocount: Pseudocount for distance
    
    Returns:
        Estimated contact frequency
    """
    return hic_scale / np.power(distance + hic_pseudocount, hic_gamma)


def geometric_mean(values: List[float], weights: Optional[List[float]] = None) -> float:
    """
    Calculate geometric mean.
    
    Args:
        values: List of values
        weights: Optional weights for weighted geometric mean
    
    Returns:
        Geometric mean
    """
    values = np.array(values)
    values = values[values > 0]  # Remove zeros
    
    if len(values) == 0:
        return 0.0
    
    if weights is None:
        return np.exp(np.mean(np.log(values)))
    else:
        weights = np.array(weights[:len(values)])
        weights = weights / np.sum(weights)
        return np.exp(np.sum(weights * np.log(values)))


def weighted_geometric_mean(values: List[float], weights: List[float]) -> float:
    """
    Calculate weighted geometric mean: prod(v_i ^ w_i).
    
    Args:
        values: List of values
        weights: Weights for each value
    
    Returns:
        Weighted geometric mean
    """
    result = 1.0
    for v, w in zip(values, weights):
        if v > 0:
            result *= np.power(v, w)
    return result


def load_chromosome_sizes(chrom_sizes_file: str) -> Dict[str, int]:
    """
    Load chromosome sizes from file.
    
    Args:
        chrom_sizes_file: Path to chrom.sizes file
    
    Returns:
        Dictionary of chromosome -> size
    """
    sizes = {}
    with open(chrom_sizes_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                sizes[parts[0]] = int(parts[1])
    return sizes


def filter_chromosomes(df: pd.DataFrame,
                       chrom_sizes: Dict[str, int],
                       chrom_column: str = 'chr') -> pd.DataFrame:
    """
    Filter DataFrame to valid chromosomes.
    
    Args:
        df: DataFrame with chromosome column
        chrom_sizes: Valid chromosome sizes
        chrom_column: Name of chromosome column
    
    Returns:
        Filtered DataFrame
    """
    valid_chroms = set(chrom_sizes.keys())
    return df[df[chrom_column].isin(valid_chroms)].copy()


def merge_overlapping_regions(df: pd.DataFrame,
                              distance: int = 0) -> pd.DataFrame:
    """
    Merge overlapping or nearby regions.
    
    Args:
        df: DataFrame with chr, start, end columns
        distance: Maximum distance to merge
    
    Returns:
        DataFrame with merged regions
    """
    # Write to temp file
    with tempfile.NamedTemporaryFile(suffix='.bed', delete=False, mode='w') as tmp_in:
        df[['chr', 'start', 'end']].to_csv(tmp_in, sep='\t', header=False, index=False)
        tmp_in_name = tmp_in.name
    
    # Run bedtools merge
    with tempfile.NamedTemporaryFile(suffix='.bed', delete=False) as tmp_out:
        cmd = f"sort -k1,1 -k2,2n {tmp_in_name} | bedtools merge -d {distance} > {tmp_out.name}"
        run_command(cmd)
        merged = pd.read_csv(tmp_out.name, sep='\t', header=None, names=['chr', 'start', 'end'])
        os.unlink(tmp_out.name)
    
    os.unlink(tmp_in_name)
    return merged


def create_enhancer_gene_pairs(enhancers: pd.DataFrame,
                               genes: pd.DataFrame,
                               max_distance: int = 5000000) -> pd.DataFrame:
    """
    Create all enhancer-gene pairs within distance threshold.
    
    Args:
        enhancers: DataFrame with enhancer regions
        genes: DataFrame with gene info (must have TSS column)
        max_distance: Maximum E-G distance
    
    Returns:
        DataFrame with E-G pairs
    """
    pairs = []
    
    for _, enh in enhancers.iterrows():
        enh_mid = (enh['start'] + enh['end']) // 2
        enh_chr = enh['chr']
        
        # Find genes on same chromosome within distance
        chr_genes = genes[genes['chr'] == enh_chr]
        
        for _, gene in chr_genes.iterrows():
            tss = gene.get('TSS', gene.get('start', 0))
            dist = abs(enh_mid - tss)
            
            if dist <= max_distance:
                pairs.append({
                    'chr': enh_chr,
                    'start': enh['start'],
                    'end': enh['end'],
                    'name': enh.get('name', f"{enh_chr}:{enh['start']}-{enh['end']}"),
                    'TargetGene': gene.get('gene_name', gene.get('name', '')),
                    'TargetGeneTSS': tss,
                    'distance': dist
                })
    
    return pd.DataFrame(pairs)


def assign_enhancer_class(distance: int, 
                          is_promoter_overlap: bool = False) -> str:
    """
    Assign enhancer class based on distance to TSS.
    
    Args:
        distance: Distance to TSS
        is_promoter_overlap: Whether enhancer overlaps promoter
    
    Returns:
        Class label
    """
    if is_promoter_overlap or distance < 500:
        return 'promoter'
    elif distance < 2000:
        return 'proximal'
    else:
        return 'distal'


def safe_divide(numerator: Union[float, np.ndarray],
                denominator: Union[float, np.ndarray],
                default: float = 0.0) -> Union[float, np.ndarray]:
    """
    Safe division avoiding divide by zero.
    
    Args:
        numerator: Numerator value(s)
        denominator: Denominator value(s)
        default: Default value when denominator is zero
    
    Returns:
        Result of division
    """
    if isinstance(denominator, np.ndarray):
        result = np.where(denominator != 0, numerator / denominator, default)
    else:
        result = numerator / denominator if denominator != 0 else default
    return result


def log_transform(values: np.ndarray, 
                  pseudocount: float = 1.0,
                  base: int = 2) -> np.ndarray:
    """
    Log transform values.
    
    Args:
        values: Values to transform
        pseudocount: Pseudocount to add
        base: Log base (2 or 10)
    
    Returns:
        Log-transformed values
    """
    if base == 2:
        return np.log2(values + pseudocount)
    elif base == 10:
        return np.log10(values + pseudocount)
    else:
        return np.log(values + pseudocount)


def ensure_dir(path: str) -> None:
    """Create directory if it doesn't exist."""
    os.makedirs(path, exist_ok=True)


def get_file_handle(filename: str, mode: str = 'r'):
    """Get file handle, handling gzip automatically."""
    if filename.endswith('.gz'):
        return gzip.open(filename, mode + 't')
    else:
        return open(filename, mode)


class PACELogger:
    """Simple logger for PACE pipeline."""
    
    def __init__(self, verbose: bool = True):
        self.verbose = verbose
    
    def info(self, msg: str) -> None:
        if self.verbose:
            print(f"[PACE INFO] {msg}", file=sys.stderr)
    
    def warning(self, msg: str) -> None:
        print(f"[PACE WARNING] {msg}", file=sys.stderr)
    
    def error(self, msg: str) -> None:
        print(f"[PACE ERROR] {msg}", file=sys.stderr)


# Global logger instance
logger = PACELogger()
