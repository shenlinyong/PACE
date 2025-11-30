#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PACE Peak Processing Module

Functions for processing ATAC-seq/DNase-seq peaks and creating candidate regions.

Author: Linyong Shen @ Northwest A&F University
"""

import os
import numpy as np
import pandas as pd
from typing import Optional, List, Dict, Tuple
import tempfile

from tools import (
    read_bed, write_bed, run_command, bedtools_intersect,
    filter_chromosomes, merge_overlapping_regions, logger
)


def load_narrowpeak(filename: str) -> pd.DataFrame:
    """
    Load MACS2 narrowPeak file.
    
    Args:
        filename: Path to narrowPeak file
    
    Returns:
        DataFrame with peak data
    """
    columns = ['chr', 'start', 'end', 'name', 'score', 'strand', 
               'signalValue', 'pValue', 'qValue', 'peak']
    
    df = pd.read_csv(filename, sep='\t', header=None, names=columns)
    
    # Calculate summit position
    df['summit'] = df['start'] + df['peak']
    
    return df


def select_strongest_peaks(peaks: pd.DataFrame,
                          n_peaks: int = 150000,
                          score_column: str = 'signalValue') -> pd.DataFrame:
    """
    Select N strongest peaks by signal value.
    
    Args:
        peaks: DataFrame with peak data
        n_peaks: Number of peaks to retain
        score_column: Column to use for ranking
    
    Returns:
        DataFrame with top peaks
    """
    if len(peaks) <= n_peaks:
        return peaks
    
    return peaks.nlargest(n_peaks, score_column).copy()


def extend_peaks_from_summit(peaks: pd.DataFrame,
                            extend_bp: int = 250,
                            chrom_sizes: Optional[Dict[str, int]] = None) -> pd.DataFrame:
    """
    Create candidate regions by extending from peak summit.
    
    Args:
        peaks: DataFrame with summit column
        extend_bp: Base pairs to extend in each direction
        chrom_sizes: Chromosome sizes for boundary checking
    
    Returns:
        DataFrame with extended regions
    """
    regions = peaks.copy()
    
    # Create new start/end centered on summit
    regions['start'] = regions['summit'] - extend_bp
    regions['end'] = regions['summit'] + extend_bp
    
    # Ensure non-negative start
    regions['start'] = regions['start'].clip(lower=0)
    
    # Clip to chromosome boundaries if sizes provided
    if chrom_sizes:
        for chrom, size in chrom_sizes.items():
            mask = regions['chr'] == chrom
            regions.loc[mask, 'end'] = regions.loc[mask, 'end'].clip(upper=size)
    
    return regions


def remove_blacklist_regions(regions: pd.DataFrame,
                            blacklist_file: str) -> pd.DataFrame:
    """
    Remove regions overlapping blacklist.
    
    Args:
        regions: DataFrame with regions
        blacklist_file: Path to blacklist BED file
    
    Returns:
        Filtered DataFrame
    """
    if not os.path.exists(blacklist_file) or os.path.getsize(blacklist_file) == 0:
        logger.info("No blacklist file or empty, skipping blacklist filtering")
        return regions
    
    # Write regions to temp file
    with tempfile.NamedTemporaryFile(suffix='.bed', delete=False, mode='w') as tmp:
        regions[['chr', 'start', 'end', 'name']].to_csv(
            tmp, sep='\t', header=False, index=False
        )
        tmp_name = tmp.name
    
    # Use bedtools to find non-overlapping regions
    with tempfile.NamedTemporaryFile(suffix='.bed', delete=False) as tmp_out:
        cmd = f"bedtools intersect -a {tmp_name} -b {blacklist_file} -v > {tmp_out.name}"
        run_command(cmd)
        
        if os.path.getsize(tmp_out.name) > 0:
            filtered = pd.read_csv(
                tmp_out.name, sep='\t', header=None,
                names=['chr', 'start', 'end', 'name']
            )
        else:
            filtered = pd.DataFrame(columns=['chr', 'start', 'end', 'name'])
        
        os.unlink(tmp_out.name)
    
    os.unlink(tmp_name)
    
    # Merge back other columns
    if len(filtered) > 0:
        result = regions[regions['name'].isin(filtered['name'])].copy()
    else:
        result = pd.DataFrame(columns=regions.columns)
    
    logger.info(f"Removed {len(regions) - len(result)} blacklisted regions")
    
    return result


def remove_tss_regions(regions: pd.DataFrame,
                       tss_file: str,
                       keep_proximal: bool = True) -> pd.DataFrame:
    """
    Remove or mark regions overlapping TSS.
    
    Args:
        regions: DataFrame with regions
        tss_file: Path to TSS BED file
        keep_proximal: If True, mark instead of remove
    
    Returns:
        Filtered/annotated DataFrame
    """
    # Write regions to temp file
    with tempfile.NamedTemporaryFile(suffix='.bed', delete=False, mode='w') as tmp:
        regions[['chr', 'start', 'end', 'name']].to_csv(
            tmp, sep='\t', header=False, index=False
        )
        tmp_name = tmp.name
    
    # Find overlapping TSS
    with tempfile.NamedTemporaryFile(suffix='.bed', delete=False) as tmp_out:
        cmd = f"bedtools intersect -a {tmp_name} -b {tss_file} -wa -u > {tmp_out.name}"
        run_command(cmd)
        
        if os.path.getsize(tmp_out.name) > 0:
            tss_overlap = pd.read_csv(
                tmp_out.name, sep='\t', header=None, names=['chr', 'start', 'end', 'name']
            )
            tss_overlap_names = set(tss_overlap['name'])
        else:
            tss_overlap_names = set()
        
        os.unlink(tmp_out.name)
    
    os.unlink(tmp_name)
    
    if keep_proximal:
        # Mark TSS-overlapping regions
        regions = regions.copy()
        regions['isTSS'] = regions['name'].isin(tss_overlap_names)
        return regions
    else:
        # Remove TSS-overlapping regions
        return regions[~regions['name'].isin(tss_overlap_names)].copy()


def quantile_normalize_signal(signal: pd.Series,
                             reference: pd.Series) -> pd.Series:
    """
    Quantile normalize signal to reference distribution.
    
    Args:
        signal: Signal values to normalize
        reference: Reference distribution
    
    Returns:
        Normalized signal
    """
    # Rank the signal values
    ranked = signal.rank(method='average')
    
    # Get reference quantiles
    ref_sorted = reference.sort_values()
    
    # Map ranks to reference values
    n_signal = len(signal)
    n_ref = len(reference)
    
    # Interpolate reference values
    ref_indices = (ranked - 1) * (n_ref - 1) / (n_signal - 1)
    normalized = np.interp(ref_indices, np.arange(n_ref), ref_sorted.values)
    
    return pd.Series(normalized, index=signal.index)


def count_reads_in_regions(regions: pd.DataFrame,
                          bam_or_tagalign: str,
                          output_column: str = 'readCount') -> pd.DataFrame:
    """
    Count reads overlapping each region.
    
    Args:
        regions: DataFrame with regions
        bam_or_tagalign: Path to BAM or tagAlign file
        output_column: Name for count column
    
    Returns:
        DataFrame with read counts
    """
    # Write regions to temp file
    with tempfile.NamedTemporaryFile(suffix='.bed', delete=False, mode='w') as tmp:
        regions[['chr', 'start', 'end', 'name']].to_csv(
            tmp, sep='\t', header=False, index=False
        )
        tmp_name = tmp.name
    
    # Run bedtools coverage
    with tempfile.NamedTemporaryFile(suffix='.bed', delete=False) as tmp_out:
        cmd = f"bedtools coverage -a {tmp_name} -b {bam_or_tagalign} -counts > {tmp_out.name}"
        run_command(cmd)
        
        counts = pd.read_csv(
            tmp_out.name, sep='\t', header=None,
            names=['chr', 'start', 'end', 'name', output_column]
        )
        
        os.unlink(tmp_out.name)
    
    os.unlink(tmp_name)
    
    # Merge counts back
    result = regions.copy()
    result = result.merge(
        counts[['name', output_column]], 
        on='name', 
        how='left'
    )
    result[output_column] = result[output_column].fillna(0)
    
    return result


def calculate_rpm(counts: pd.Series, total_reads: int) -> pd.Series:
    """
    Calculate Reads Per Million.
    
    Args:
        counts: Read counts
        total_reads: Total mapped reads
    
    Returns:
        RPM values
    """
    return counts * 1e6 / total_reads


def create_candidate_regions(narrowpeak_file: str,
                            chrom_sizes_file: str,
                            output_file: str,
                            blacklist_file: Optional[str] = None,
                            tss_file: Optional[str] = None,
                            n_peaks: int = 150000,
                            extend_bp: int = 250,
                            min_region_size: int = 500) -> pd.DataFrame:
    """
    Create candidate enhancer regions from MACS2 peaks.
    
    This is the main function for generating candidate regions.
    
    Args:
        narrowpeak_file: MACS2 narrowPeak file
        chrom_sizes_file: Chromosome sizes file
        output_file: Output BED file
        blacklist_file: Blacklist regions file
        tss_file: TSS regions file
        n_peaks: Number of strongest peaks to retain
        extend_bp: Extension from summit
        min_region_size: Minimum region size
    
    Returns:
        DataFrame with candidate regions
    """
    logger.info(f"Loading peaks from {narrowpeak_file}")
    peaks = load_narrowpeak(narrowpeak_file)
    logger.info(f"Loaded {len(peaks)} peaks")
    
    # Load chromosome sizes
    chrom_sizes = {}
    with open(chrom_sizes_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                chrom_sizes[parts[0]] = int(parts[1])
    
    # Filter to valid chromosomes
    peaks = filter_chromosomes(peaks, chrom_sizes)
    logger.info(f"After chromosome filter: {len(peaks)} peaks")
    
    # Select strongest peaks
    peaks = select_strongest_peaks(peaks, n_peaks)
    logger.info(f"Selected top {len(peaks)} peaks")
    
    # Extend from summit
    regions = extend_peaks_from_summit(peaks, extend_bp, chrom_sizes)
    
    # Ensure minimum size
    region_sizes = regions['end'] - regions['start']
    small_regions = region_sizes < min_region_size
    if small_regions.any():
        # Extend small regions
        midpoints = (regions.loc[small_regions, 'start'] + regions.loc[small_regions, 'end']) // 2
        regions.loc[small_regions, 'start'] = midpoints - min_region_size // 2
        regions.loc[small_regions, 'end'] = midpoints + min_region_size // 2
        regions.loc[small_regions, 'start'] = regions.loc[small_regions, 'start'].clip(lower=0)
    
    # Remove blacklist regions
    if blacklist_file:
        regions = remove_blacklist_regions(regions, blacklist_file)
        logger.info(f"After blacklist filter: {len(regions)} regions")
    
    # Mark TSS-overlapping regions
    if tss_file:
        regions = remove_tss_regions(regions, tss_file, keep_proximal=True)
        logger.info(f"Marked {regions['isTSS'].sum()} TSS-overlapping regions")
    
    # Create unique region names
    regions['name'] = [f"region_{i}" for i in range(len(regions))]
    
    # Sort by position
    regions = regions.sort_values(['chr', 'start']).reset_index(drop=True)
    
    # Write output
    output_cols = ['chr', 'start', 'end', 'name']
    if 'isTSS' in regions.columns:
        output_cols.append('isTSS')
    
    regions[output_cols].to_csv(output_file, sep='\t', header=False, index=False)
    logger.info(f"Wrote {len(regions)} candidate regions to {output_file}")
    
    return regions


def merge_multiple_peak_files(peak_files: List[str],
                             output_file: str,
                             merge_distance: int = 0) -> pd.DataFrame:
    """
    Merge peaks from multiple samples.
    
    Args:
        peak_files: List of narrowPeak files
        output_file: Output merged BED file
        merge_distance: Distance for merging nearby peaks
    
    Returns:
        DataFrame with merged peaks
    """
    all_peaks = []
    
    for pf in peak_files:
        peaks = load_narrowpeak(pf)
        all_peaks.append(peaks[['chr', 'start', 'end']])
    
    combined = pd.concat(all_peaks, ignore_index=True)
    
    # Merge overlapping regions
    merged = merge_overlapping_regions(combined, merge_distance)
    
    # Add names
    merged['name'] = [f"merged_peak_{i}" for i in range(len(merged))]
    
    merged.to_csv(output_file, sep='\t', header=False, index=False)
    
    return merged
