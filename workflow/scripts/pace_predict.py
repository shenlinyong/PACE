#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PACE Prediction

Predict enhancer-gene regulatory interactions.

Usage:
    python pace_predict.py \
        --enhancers EnhancerList.txt \
        --genes GeneList.txt \
        --output predictions.tsv.gz \
        [--hic_file sample.hic] \
        [--hic_type hic] \
        [--expression expression.tsv] \
        [--use_expression_weight] \
        [--threshold 0.02]

Author: Linyong Shen @ Northwest A&F University
"""

import argparse
import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from predictor import run_predictions
from tools import logger


def main():
    parser = argparse.ArgumentParser(
        description='Run PACE enhancer-gene prediction'
    )
    
    # Required arguments
    parser.add_argument('--enhancers', required=True,
                       help='Enhancer list file from neighborhood analysis')
    parser.add_argument('--genes', required=True,
                       help='Gene list file from neighborhood analysis')
    parser.add_argument('--output', '-o', required=True,
                       help='Output predictions file')
    
    # Hi-C options
    parser.add_argument('--hic_file', default=None,
                       help='Hi-C file for 3D contact (optional)')
    parser.add_argument('--hic_type', default='hic',
                       choices=['hic', 'bedpe', 'cool', 'avg'],
                       help='Type of Hi-C file')
    parser.add_argument('--hic_resolution', type=int, default=5000,
                       help='Hi-C resolution in bp (default: 5000)')
    
    # Power-law parameters
    parser.add_argument('--hic_gamma', type=float, default=1.024238616787792,
                       help='Power-law exponent for distance decay')
    parser.add_argument('--hic_scale', type=float, default=5.9594510043736655,
                       help='Power-law scale factor')
    parser.add_argument('--scale_hic_using_powerlaw', action='store_true',
                       help='Scale Hi-C values using power-law')
    
    # Expression options
    parser.add_argument('--expression', default=None,
                       help='Expression file for filtering/weighting')
    parser.add_argument('--use_expression_weight', action='store_true',
                       help='Weight predictions by gene expression')
    parser.add_argument('--weight_method', default='log',
                       choices=['binary', 'linear', 'log'],
                       help='Expression weight method')
    parser.add_argument('--min_expression', type=float, default=1.0,
                       help='Minimum expression for filtering (TPM)')
    
    # Prediction options
    parser.add_argument('--max_distance', type=int, default=5000000,
                       help='Maximum enhancer-gene distance (default: 5Mb)')
    parser.add_argument('--threshold', type=float, default=0.02,
                       help='Score threshold for filtering (default: 0.02)')
    parser.add_argument('--include_self_promoter', action='store_true', default=True,
                       help='Include self-promoter interactions')
    parser.add_argument('--no_self_promoter', action='store_true',
                       help='Exclude self-promoter interactions')
    
    # Output options
    parser.add_argument('--score_column', default='ABC.Score',
                       help='Name for score column')
    
    args = parser.parse_args()
    
    # Validate inputs
    if not os.path.exists(args.enhancers):
        logger.error(f"Enhancer file not found: {args.enhancers}")
        sys.exit(1)
    
    if not os.path.exists(args.genes):
        logger.error(f"Gene file not found: {args.genes}")
        sys.exit(1)
    
    # Create output directory
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    
    # Handle self-promoter option
    include_self = args.include_self_promoter and not args.no_self_promoter
    
    # Run predictions
    logger.info("Starting PACE prediction")
    
    output_file = run_predictions(
        enhancer_file=args.enhancers,
        gene_file=args.genes,
        output_file=args.output,
        hic_file=args.hic_file,
        hic_type=args.hic_type,
        hic_resolution=args.hic_resolution,
        expression_file=args.expression,
        use_expression_weight=args.use_expression_weight,
        weight_method=args.weight_method,
        min_expression=args.min_expression,
        max_distance=args.max_distance,
        hic_gamma=args.hic_gamma,
        hic_scale=args.hic_scale,
        score_threshold=args.threshold,
        include_self_promoter=include_self
    )
    
    logger.info(f"Predictions written to: {output_file}")
    logger.info("Done!")


if __name__ == '__main__':
    main()
