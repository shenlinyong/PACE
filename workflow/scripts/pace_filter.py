#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PACE Filter Predictions

Filter predictions by score threshold and other criteria.

Usage:
    python pace_filter.py \
        --predictions predictions.tsv.gz \
        --output filtered_predictions.tsv \
        [--threshold 0.02] \
        [--only_expressed] \
        [--output_full] \
        [--output_slim]

Author: Linyong Shen @ Northwest A&F University
"""

import argparse
import sys
import os
import pandas as pd

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from tools import logger


def filter_predictions(predictions_file: str,
                       output_file: str,
                       threshold: float = 0.02,
                       score_column: str = 'ABC.Score',
                       only_expressed: bool = False,
                       output_full: bool = True,
                       output_slim: bool = True) -> pd.DataFrame:
    """
    Filter predictions by threshold and criteria.
    
    Args:
        predictions_file: Input predictions file
        output_file: Output filtered file
        threshold: Score threshold
        score_column: Score column name
        only_expressed: Filter to expressed genes only
        output_full: Write full output with all columns
        output_slim: Write slim output with key columns only
    
    Returns:
        Filtered DataFrame
    """
    # Load predictions
    logger.info(f"Loading predictions from {predictions_file}")
    
    if predictions_file.endswith('.gz'):
        df = pd.read_csv(predictions_file, sep='\t', compression='gzip')
    else:
        df = pd.read_csv(predictions_file, sep='\t')
    
    logger.info(f"Loaded {len(df)} predictions")
    
    # Filter by threshold
    filtered = df[df[score_column] >= threshold].copy()
    logger.info(f"After score filter (>= {threshold}): {len(filtered)} predictions")
    
    # Filter to expressed genes
    if only_expressed and 'isExpressed' in filtered.columns:
        filtered = filtered[filtered['isExpressed']]
        logger.info(f"After expression filter: {len(filtered)} predictions")
    
    # Sort by score
    filtered = filtered.sort_values(score_column, ascending=False)
    
    # Write full output
    if output_full:
        full_file = output_file.replace('.tsv', '_Full.tsv')
        filtered.to_csv(full_file, sep='\t', index=False)
        logger.info(f"Wrote full output to {full_file}")
    
    # Write slim output
    if output_slim:
        slim_columns = [
            'chr', 'start', 'end', 'name',
            'TargetGene', 'TargetGeneTSS', 'TargetGeneStrand',
            score_column, 'distance', 'class'
        ]
        
        # Only keep columns that exist
        slim_columns = [c for c in slim_columns if c in filtered.columns]
        
        slim_file = output_file.replace('.tsv', '.tsv')
        filtered[slim_columns].to_csv(slim_file, sep='\t', index=False)
        logger.info(f"Wrote slim output to {slim_file}")
    
    return filtered


def main():
    parser = argparse.ArgumentParser(
        description='Filter PACE predictions'
    )
    
    # Required arguments
    parser.add_argument('--predictions', required=True,
                       help='Input predictions file')
    parser.add_argument('--output', '-o', required=True,
                       help='Output filtered predictions file')
    
    # Filter options
    parser.add_argument('--threshold', type=float, default=0.02,
                       help='Score threshold (default: 0.02)')
    parser.add_argument('--score_column', default='ABC.Score',
                       help='Score column name')
    parser.add_argument('--only_expressed', action='store_true',
                       help='Filter to expressed genes only')
    
    # Output options
    parser.add_argument('--output_full', action='store_true', default=True,
                       help='Write full output with all columns')
    parser.add_argument('--output_slim', action='store_true', default=True,
                       help='Write slim output with key columns')
    parser.add_argument('--no_full', action='store_true',
                       help='Do not write full output')
    parser.add_argument('--no_slim', action='store_true',
                       help='Do not write slim output')
    
    args = parser.parse_args()
    
    # Validate inputs
    if not os.path.exists(args.predictions):
        logger.error(f"Predictions file not found: {args.predictions}")
        sys.exit(1)
    
    # Create output directory
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    
    # Run filtering
    filter_predictions(
        predictions_file=args.predictions,
        output_file=args.output,
        threshold=args.threshold,
        score_column=args.score_column,
        only_expressed=args.only_expressed,
        output_full=args.output_full and not args.no_full,
        output_slim=args.output_slim and not args.no_slim
    )
    
    logger.info("Done!")


if __name__ == '__main__':
    main()
