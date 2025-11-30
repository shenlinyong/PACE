#!/usr/bin/env python3
"""
PACE Machine Learning Integration Script

Train ML models to enhance enhancer-gene predictions using validated data.

Usage:
    # Train model with eQTL validation data
    python scripts/pace_ml.py train \
        --predictions results/Sample/Predictions/EnhancerPredictionsAllPutative.tsv.gz \
        --validation eqtl_validated_pairs.tsv \
        --output models/sample_model.pkl
    
    # Apply trained model to new predictions
    python scripts/pace_ml.py predict \
        --predictions results/NewSample/Predictions/EnhancerPredictionsAllPutative.tsv.gz \
        --model models/sample_model.pkl \
        --output results/NewSample/Predictions/EnhancerPredictions_ML.tsv

Author: Linyong Shen @ Northwest A&F University
"""

import argparse
import sys
import os
import pickle
import numpy as np
import pandas as pd
import logging

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from ml_integration import (
    PACEMLPredictor, MLConfig, FeatureEngineering,
    train_pace_ml_model, SKLEARN_AVAILABLE
)

logging.basicConfig(level=logging.INFO, format='[PACE ML] %(message)s')
logger = logging.getLogger(__name__)


def load_predictions(filepath: str) -> pd.DataFrame:
    """Load predictions file (supports .gz)"""
    if filepath.endswith('.gz'):
        return pd.read_csv(filepath, sep='\t', compression='gzip')
    return pd.read_csv(filepath, sep='\t')


def load_validation(filepath: str, 
                   enhancer_col: str = 'enhancer',
                   gene_col: str = 'gene',
                   label_col: str = None) -> pd.DataFrame:
    """
    Load validation data (e.g., eQTL pairs).
    
    Expected format:
        - enhancer: enhancer region (chr:start-end or chr_start_end)
        - gene: gene name
        - validated: 1 for positive, 0 for negative (optional)
    """
    df = pd.read_csv(filepath, sep='\t')
    
    # Standardize column names
    if enhancer_col in df.columns:
        df['enhancer_id'] = df[enhancer_col]
    if gene_col in df.columns:
        df['gene_id'] = df[gene_col]
    
    # Set validation label
    if label_col and label_col in df.columns:
        df['validated'] = df[label_col]
    elif 'validated' not in df.columns:
        # Assume all provided pairs are validated positives
        df['validated'] = 1
    
    return df


def merge_with_validation(predictions: pd.DataFrame, 
                         validation: pd.DataFrame) -> pd.DataFrame:
    """
    Merge predictions with validation data to create training set.
    """
    # Create enhancer ID for matching
    predictions['enhancer_id'] = (
        predictions['chr'].astype(str) + ':' + 
        predictions['start'].astype(str) + '-' + 
        predictions['end'].astype(str)
    )
    
    # Alternative format
    predictions['enhancer_id_alt'] = (
        predictions['chr'].astype(str) + '_' + 
        predictions['start'].astype(str) + '_' + 
        predictions['end'].astype(str)
    )
    
    # Get gene column
    gene_col = 'TargetGene' if 'TargetGene' in predictions.columns else 'gene_name'
    
    # Create validation set lookup
    val_set = set()
    for _, row in validation.iterrows():
        if row.get('validated', 1) == 1:
            val_set.add((row['enhancer_id'], row['gene_id']))
    
    # Label predictions
    def is_validated(row):
        for enh_id in [row['enhancer_id'], row['enhancer_id_alt']]:
            if (enh_id, row[gene_col]) in val_set:
                return 1
        return 0
    
    predictions['validated'] = predictions.apply(is_validated, axis=1)
    
    logger.info(f"Matched {predictions['validated'].sum()} validated E-G pairs")
    
    return predictions


def get_feature_columns(df: pd.DataFrame) -> list:
    """Auto-detect feature columns"""
    exclude = {
        'chr', 'start', 'end', 'name', 'TargetGene', 'TargetGeneEnsemblID',
        'TargetGeneTSS', 'TargetGeneStrand', 'class', 'validated',
        'enhancer_id', 'enhancer_id_alt', 'gene_id',
        'ABC.Score', 'PACE.Score', 'ML.Score'
    }
    
    numeric_cols = df.select_dtypes(include=[np.number]).columns
    feature_cols = [c for c in numeric_cols if c not in exclude]
    
    return feature_cols


def train_model(args):
    """Train ML model from predictions and validation data"""
    if not SKLEARN_AVAILABLE:
        logger.error("scikit-learn not installed. Run: pip install scikit-learn")
        sys.exit(1)
    
    logger.info("Loading predictions...")
    predictions = load_predictions(args.predictions)
    logger.info(f"Loaded {len(predictions)} predictions")
    
    logger.info("Loading validation data...")
    validation = load_validation(
        args.validation,
        enhancer_col=args.enhancer_column,
        gene_col=args.gene_column
    )
    logger.info(f"Loaded {len(validation)} validation pairs")
    
    # Merge
    logger.info("Merging predictions with validation...")
    training_data = merge_with_validation(predictions, validation)
    
    n_positive = training_data['validated'].sum()
    n_negative = len(training_data) - n_positive
    
    if n_positive < 10:
        logger.error(f"Too few validated pairs ({n_positive}). Need at least 10.")
        sys.exit(1)
    
    logger.info(f"Training set: {n_positive} positive, {n_negative} negative")
    
    # Balance classes if needed
    if args.balance_classes and n_negative > n_positive * 10:
        logger.info("Balancing classes by downsampling negatives...")
        pos_data = training_data[training_data['validated'] == 1]
        neg_data = training_data[training_data['validated'] == 0].sample(
            n=min(n_positive * 5, n_negative), random_state=42
        )
        training_data = pd.concat([pos_data, neg_data])
        logger.info(f"Balanced set: {len(pos_data)} positive, {len(neg_data)} negative")
    
    # Feature engineering
    logger.info("Engineering features...")
    training_data = FeatureEngineering.create_features(training_data)
    
    # Get features
    feature_cols = args.features.split(',') if args.features else get_feature_columns(training_data)
    logger.info(f"Using {len(feature_cols)} features: {feature_cols[:5]}...")
    
    # Configure model
    config = MLConfig(
        model_type=args.model_type,
        n_estimators=args.n_estimators,
        max_depth=args.max_depth,
        learning_rate=args.learning_rate,
        n_cv_folds=args.cv_folds
    )
    
    # Train
    model, metrics = train_pace_ml_model(
        training_data,
        label_column='validated',
        feature_columns=feature_cols,
        config=config
    )
    
    # Save model
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    
    with open(args.output, 'wb') as f:
        pickle.dump({
            'model': model,
            'feature_columns': feature_cols,
            'metrics': metrics,
            'config': config
        }, f)
    
    logger.info(f"Model saved to: {args.output}")
    
    # Print results
    print("\n" + "=" * 50)
    print("Training Results")
    print("=" * 50)
    print(f"CV AUC: {metrics['cv_auc_mean']:.3f} ± {metrics['cv_auc_std']:.3f}")
    print(f"Samples: {metrics['n_samples']}")
    print(f"Features: {metrics['n_features']}")
    print("\nTop Feature Importance:")
    importance = model.get_feature_importance()
    for _, row in importance.head(10).iterrows():
        print(f"  {row['feature']}: {row['importance']:.4f}")
    print("=" * 50)


def predict_with_model(args):
    """Apply trained model to predictions"""
    if not SKLEARN_AVAILABLE:
        logger.error("scikit-learn not installed. Run: pip install scikit-learn")
        sys.exit(1)
    
    # Load model
    logger.info(f"Loading model from {args.model}...")
    with open(args.model, 'rb') as f:
        saved = pickle.load(f)
    
    model = saved['model']
    feature_cols = saved['feature_columns']
    
    logger.info(f"Model uses {len(feature_cols)} features")
    
    # Load predictions
    logger.info("Loading predictions...")
    predictions = load_predictions(args.predictions)
    logger.info(f"Loaded {len(predictions)} predictions")
    
    # Feature engineering
    predictions = FeatureEngineering.create_features(predictions)
    
    # Check features
    missing_features = [f for f in feature_cols if f not in predictions.columns]
    if missing_features:
        logger.warning(f"Missing features (will use 0): {missing_features}")
        for f in missing_features:
            predictions[f] = 0
    
    # Prepare features
    X, _ = model.prepare_features(predictions, feature_cols)
    
    # Predict
    logger.info("Generating ML predictions...")
    probabilities, _ = model.predict(X)
    
    # Add ML score
    predictions['ML.Score'] = probabilities
    
    # Combine with ABC score if available
    if 'ABC.Score' in predictions.columns:
        # Weighted combination
        alpha = args.abc_weight
        predictions['Combined.Score'] = (
            alpha * predictions['ABC.Score'] + 
            (1 - alpha) * predictions['ML.Score']
        )
        score_col = 'Combined.Score'
    else:
        score_col = 'ML.Score'
    
    # Sort by score
    predictions = predictions.sort_values(score_col, ascending=False)
    
    # Save
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    predictions.to_csv(args.output, sep='\t', index=False)
    
    logger.info(f"Predictions saved to: {args.output}")
    
    # Print summary
    print("\n" + "=" * 50)
    print("Prediction Results")
    print("=" * 50)
    print(f"Total predictions: {len(predictions)}")
    print(f"ML Score range: {predictions['ML.Score'].min():.4f} - {predictions['ML.Score'].max():.4f}")
    if 'Combined.Score' in predictions.columns:
        print(f"Combined Score range: {predictions['Combined.Score'].min():.4f} - {predictions['Combined.Score'].max():.4f}")
    print("=" * 50)


def main():
    parser = argparse.ArgumentParser(
        description='PACE Machine Learning Integration',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Train model with eQTL validation
  python scripts/pace_ml.py train \\
      --predictions results/Sample/Predictions/EnhancerPredictionsAllPutative.tsv.gz \\
      --validation data/eqtl_pairs.tsv \\
      --output models/eqtl_model.pkl

  # Apply model to new predictions
  python scripts/pace_ml.py predict \\
      --predictions results/NewSample/Predictions/EnhancerPredictionsAllPutative.tsv.gz \\
      --model models/eqtl_model.pkl \\
      --output results/NewSample/Predictions/EnhancerPredictions_ML.tsv
        """
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Commands')
    
    # Train command
    train_parser = subparsers.add_parser('train', help='Train ML model')
    train_parser.add_argument('--predictions', required=True,
                             help='PACE predictions file')
    train_parser.add_argument('--validation', required=True,
                             help='Validation data (e.g., eQTL pairs)')
    train_parser.add_argument('--output', '-o', required=True,
                             help='Output model file (.pkl)')
    train_parser.add_argument('--enhancer_column', default='enhancer',
                             help='Enhancer column in validation file')
    train_parser.add_argument('--gene_column', default='gene',
                             help='Gene column in validation file')
    train_parser.add_argument('--features', default=None,
                             help='Comma-separated feature columns (auto-detect if not specified)')
    train_parser.add_argument('--model_type', default='gradient_boosting',
                             choices=['gradient_boosting', 'random_forest'],
                             help='ML model type')
    train_parser.add_argument('--n_estimators', type=int, default=100,
                             help='Number of trees')
    train_parser.add_argument('--max_depth', type=int, default=5,
                             help='Maximum tree depth')
    train_parser.add_argument('--learning_rate', type=float, default=0.1,
                             help='Learning rate (for gradient boosting)')
    train_parser.add_argument('--cv_folds', type=int, default=5,
                             help='Number of cross-validation folds')
    train_parser.add_argument('--balance_classes', action='store_true',
                             help='Balance classes by downsampling negatives')
    
    # Predict command
    predict_parser = subparsers.add_parser('predict', help='Apply ML model')
    predict_parser.add_argument('--predictions', required=True,
                               help='PACE predictions file')
    predict_parser.add_argument('--model', required=True,
                               help='Trained model file (.pkl)')
    predict_parser.add_argument('--output', '-o', required=True,
                               help='Output predictions file')
    predict_parser.add_argument('--abc_weight', type=float, default=0.5,
                               help='Weight for ABC score in combined score (0-1)')
    
    args = parser.parse_args()
    
    if args.command == 'train':
        train_model(args)
    elif args.command == 'predict':
        predict_with_model(args)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == '__main__':
    main()
