#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PACE Candidate Regions Generator

Generate candidate enhancer regions from MACS2 peaks.

Usage:
    python pace_candidate_regions.py \
        --narrowPeak peaks.narrowPeak \
        --chrom_sizes genome.chrom.sizes \
        --output candidates.bed \
        [--blacklist blacklist.bed] \
        [--tss_regions tss.bed] \
        [--nStrongestPeaks 150000] \
        [--peakExtendFromSummit 250]

Author: Linyong Shen @ Northwest A&F University
"""

import argparse
import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from peaks import create_candidate_regions
from tools import logger


def main():
    parser = argparse.ArgumentParser(
        description='Generate candidate enhancer regions from MACS2 peaks'
    )
    
    # Required arguments
    parser.add_argument('--narrowPeak', required=True,
                       help='MACS2 narrowPeak file')
    parser.add_argument('--chrom_sizes', required=True,
                       help='Chromosome sizes file')
    parser.add_argument('--output', '-o', required=True,
                       help='Output candidate regions file')
    
    # Optional arguments
    parser.add_argument('--blacklist', default=None,
                       help='Blacklist regions to exclude')
    parser.add_argument('--tss_regions', default=None,
                       help='TSS regions file for annotation')
    parser.add_argument('--nStrongestPeaks', type=int, default=150000,
                       help='Number of strongest peaks to retain (default: 150000)')
    parser.add_argument('--peakExtendFromSummit', type=int, default=250,
                       help='BP to extend from peak summit (default: 250)')
    parser.add_argument('--minPeakWidth', type=int, default=500,
                       help='Minimum peak width (default: 500)')
    
    args = parser.parse_args()
    
    # Validate inputs
    if not os.path.exists(args.narrowPeak):
        logger.error(f"narrowPeak file not found: {args.narrowPeak}")
        sys.exit(1)
    
    if not os.path.exists(args.chrom_sizes):
        logger.error(f"Chromosome sizes file not found: {args.chrom_sizes}")
        sys.exit(1)
    
    # Create output directory
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    
    # Run candidate region generation
    logger.info("Starting candidate region generation")
    
    regions = create_candidate_regions(
        narrowpeak_file=args.narrowPeak,
        chrom_sizes_file=args.chrom_sizes,
        output_file=args.output,
        blacklist_file=args.blacklist,
        tss_file=args.tss_regions,
        n_peaks=args.nStrongestPeaks,
        extend_bp=args.peakExtendFromSummit,
        min_region_size=args.minPeakWidth
    )
    
    logger.info(f"Generated {len(regions)} candidate regions")
    logger.info("Done!")


if __name__ == '__main__':
    main()
