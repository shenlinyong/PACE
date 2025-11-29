#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PACE Metrics Generator

Generate quality control metrics and plots for predictions.

Usage:
    python pace_metrics.py \
        --predictions predictions.tsv.gz \
        --output_dir output/Metrics \
        [--sample_name sample1] \
        [--score_column ABC.Score]

Author: Linyong Shen @ Northwest A&F University
"""

import argparse
import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from metrics import generate_metrics, compare_samples
from tools import logger


def main():
    parser = argparse.ArgumentParser(
        description='Generate PACE quality control metrics'
    )
    
    # Required arguments
    parser.add_argument('--predictions', required=True,
                       help='Predictions file')
    parser.add_argument('--output_dir', '-o', required=True,
                       help='Output directory for metrics')
    
    # Optional arguments
    parser.add_argument('--sample_name', default='sample',
                       help='Sample name for output files')
    parser.add_argument('--score_column', default='ABC.Score',
                       help='Score column name')
    
    # Comparison mode
    parser.add_argument('--compare', nargs='+',
                       help='Additional prediction files for comparison')
    parser.add_argument('--compare_names', nargs='+',
                       help='Names for comparison samples')
    
    args = parser.parse_args()
    
    # Validate inputs
    if not os.path.exists(args.predictions):
        logger.error(f"Predictions file not found: {args.predictions}")
        sys.exit(1)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Generate metrics
    logger.info("Generating quality control metrics")
    
    metrics_file, plot_file = generate_metrics(
        predictions_file=args.predictions,
        output_dir=args.output_dir,
        score_column=args.score_column,
        sample_name=args.sample_name
    )
    
    logger.info(f"Metrics saved to: {metrics_file}")
    logger.info(f"Plots saved to: {plot_file}")
    
    # Run comparison if requested
    if args.compare:
        logger.info("Running sample comparison")
        
        pred_files = {args.sample_name: args.predictions}
        
        for i, pred_file in enumerate(args.compare):
            if os.path.exists(pred_file):
                if args.compare_names and i < len(args.compare_names):
                    name = args.compare_names[i]
                else:
                    name = f"sample_{i+2}"
                pred_files[name] = pred_file
        
        comparison_file = os.path.join(args.output_dir, 'QC_Comparison.tsv')
        compare_samples(pred_files, comparison_file, args.score_column)
        logger.info(f"Comparison saved to: {comparison_file}")
    
    logger.info("Done!")


if __name__ == '__main__':
    main()
