#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PACE Neighborhood Analysis

Quantify chromatin features in candidate enhancer regions.

Usage:
    python pace_neighborhoods.py \
        --candidate_regions candidates.bed \
        --genes genes.bed \
        --chrom_sizes genome.chrom.sizes \
        --output_dir output/ \
        --accessibility_file atac.tagAlign.gz \
        [--H3K27ac h3k27ac.bam] \
        [--H3K4me1 h3k4me1.bam] \
        [--methylation methylation.bw] \
        [--expression expression.tsv] \
        [--activity_method geometric_mean]

Author: Linyong Shen @ Northwest A&F University
"""

import argparse
import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from neighborhoods import run_neighborhoods
from tools import logger


def main():
    parser = argparse.ArgumentParser(
        description='Run PACE neighborhood analysis'
    )
    
    # Required arguments
    parser.add_argument('--candidate_regions', required=True,
                       help='Candidate enhancer regions file')
    parser.add_argument('--genes', required=True,
                       help='Gene annotation file')
    parser.add_argument('--chrom_sizes', required=True,
                       help='Chromosome sizes file')
    parser.add_argument('--output_dir', '-o', required=True,
                       help='Output directory')
    
    # Accessibility (required)
    parser.add_argument('--accessibility_file', required=True,
                       help='ATAC-seq or DNase-seq file')
    parser.add_argument('--accessibility_type', default='ATAC',
                       choices=['ATAC', 'DHS'],
                       help='Type of accessibility data')
    
    # Histone modifications (optional)
    parser.add_argument('--H3K27ac', default=None,
                       help='H3K27ac ChIP-seq file')
    parser.add_argument('--H3K4me1', default=None,
                       help='H3K4me1 ChIP-seq file')
    parser.add_argument('--H3K4me3', default=None,
                       help='H3K4me3 ChIP-seq file')
    parser.add_argument('--H3K36me3', default=None,
                       help='H3K36me3 ChIP-seq file')
    parser.add_argument('--H3K9ac', default=None,
                       help='H3K9ac ChIP-seq file')
    parser.add_argument('--H3K27me3', default=None,
                       help='H3K27me3 ChIP-seq file (inhibitory)')
    parser.add_argument('--H3K9me3', default=None,
                       help='H3K9me3 ChIP-seq file (inhibitory)')
    
    # Methylation (optional)
    parser.add_argument('--methylation', default=None,
                       help='DNA methylation file (bigWig or bedGraph)')
    
    # Expression (optional)
    parser.add_argument('--expression', default=None,
                       help='RNA-seq expression file')
    
    # Activity calculation options
    parser.add_argument('--activity_method', default='geometric_mean',
                       choices=['geometric_mean', 'weighted_geometric', 'weighted_sum', 'arithmetic_mean'],
                       help='Activity calculation method')
    
    # Weights
    parser.add_argument('--accessibility_weight', type=float, default=1.5,
                       help='Weight for accessibility signal')
    parser.add_argument('--H3K27ac_weight', type=float, default=1.0,
                       help='Weight for H3K27ac signal')
    parser.add_argument('--H3K4me1_weight', type=float, default=0.8,
                       help='Weight for H3K4me1 signal')
    parser.add_argument('--H3K4me3_weight', type=float, default=0.5,
                       help='Weight for H3K4me3 signal')
    parser.add_argument('--methylation_weight', type=float, default=0.5,
                       help='Weight for methylation (inhibitory)')
    
    args = parser.parse_args()
    
    # Validate inputs
    if not os.path.exists(args.candidate_regions):
        logger.error(f"Candidate regions file not found: {args.candidate_regions}")
        sys.exit(1)
    
    if not os.path.exists(args.genes):
        logger.error(f"Genes file not found: {args.genes}")
        sys.exit(1)
    
    # Collect histone files
    histone_files = {}
    histone_weights = {}
    
    if args.H3K27ac and os.path.exists(args.H3K27ac):
        histone_files['H3K27ac'] = args.H3K27ac
        histone_weights['H3K27ac_signal'] = args.H3K27ac_weight
    
    if args.H3K4me1 and os.path.exists(args.H3K4me1):
        histone_files['H3K4me1'] = args.H3K4me1
        histone_weights['H3K4me1_signal'] = args.H3K4me1_weight
    
    if args.H3K4me3 and os.path.exists(args.H3K4me3):
        histone_files['H3K4me3'] = args.H3K4me3
        histone_weights['H3K4me3_signal'] = args.H3K4me3_weight
    
    if args.H3K36me3 and os.path.exists(args.H3K36me3):
        histone_files['H3K36me3'] = args.H3K36me3
        histone_weights['H3K36me3_signal'] = 0.3
    
    if args.H3K9ac and os.path.exists(args.H3K9ac):
        histone_files['H3K9ac'] = args.H3K9ac
        histone_weights['H3K9ac_signal'] = 0.6
    
    # Run neighborhood analysis
    logger.info("Starting neighborhood analysis")
    
    enhancer_file, gene_file = run_neighborhoods(
        candidate_regions_file=args.candidate_regions,
        genes_file=args.genes,
        chrom_sizes_file=args.chrom_sizes,
        output_dir=args.output_dir,
        accessibility_file=args.accessibility_file,
        accessibility_type=args.accessibility_type,
        histone_files=histone_files if histone_files else None,
        methylation_file=args.methylation,
        expression_file=args.expression,
        activity_method=args.activity_method,
        histone_weights=histone_weights if histone_weights else None
    )
    
    logger.info(f"Created enhancer list: {enhancer_file}")
    logger.info(f"Created gene list: {gene_file}")
    logger.info("Done!")


if __name__ == '__main__':
    main()
