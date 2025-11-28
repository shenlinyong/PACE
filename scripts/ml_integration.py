#!/usr/bin/env python3
"""
PACE: Machine Learning Integration Module

This module provides ML-based integration of multi-omics features
for enhanced enhancer-gene prediction.

Features:
1. Gradient Boosting for feature integration
2. Feature importance analysis
3. Cross-validation and model evaluation
4. Transfer learning across species

Author: Linyong Shen @ Northwest A&F University
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass
import logging
import warnings

# Suppress sklearn warnings
warnings.filterwarnings('ignore')

try:
    from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
    from sklearn.model_selection import cross_val_score, StratifiedKFold
    from sklearn.preprocessing import StandardScaler
    from sklearn.metrics import precision_recall_curve, roc_auc_score, average_precision_score
    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False
    logging.warning("scikit-learn not available. ML features disabled.")

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class MLConfig:
    """Configuration for ML model"""
    model_type: str = "gradient_boosting"
    n_estimators: int = 100
    max_depth: int = 5
    learning_rate: float = 0.1
    min_samples_split: int = 10
    n_cv_folds: int = 5
    random_state: int = 42


class PACEMLPredictor:
    """
    Machine learning-based enhancer-gene prediction.
    
    This class integrates multiple epigenomic features using ML models
    to predict enhancer-gene regulatory relationships.
    
    Innovation: Instead of fixed formulas, learn optimal feature combinations
    from data with validation.
    """
    
    def __init__(self, config: MLConfig = None):
        """
        Initialize ML predictor.
        
        Args:
            config: ML configuration
        """
        if not SKLEARN_AVAILABLE:
            raise ImportError("scikit-learn is required for ML features")
            
        self.config = config or MLConfig()
        self.model = None
        self.scaler = StandardScaler()
        self.feature_names = []
        self.feature_importance = {}
        
    def _create_model(self):
        """Create ML model based on configuration"""
        if self.config.model_type == "gradient_boosting":
            return GradientBoostingClassifier(
                n_estimators=self.config.n_estimators,
                max_depth=self.config.max_depth,
                learning_rate=self.config.learning_rate,
                min_samples_split=self.config.min_samples_split,
                random_state=self.config.random_state
            )
        elif self.config.model_type == "random_forest":
            return RandomForestClassifier(
                n_estimators=self.config.n_estimators,
                max_depth=self.config.max_depth,
                min_samples_split=self.config.min_samples_split,
                random_state=self.config.random_state
            )
        else:
            raise ValueError(f"Unknown model type: {self.config.model_type}")
    
    def prepare_features(self, 
                         eg_pairs: pd.DataFrame,
                         feature_columns: List[str]) -> Tuple[np.ndarray, List[str]]:
        """
        Prepare feature matrix from E-G pair data.
        
        Args:
            eg_pairs: DataFrame with E-G pairs and features
            feature_columns: List of feature column names
            
        Returns:
            Feature matrix and feature names
        """
        # Filter to available features
        available_features = [f for f in feature_columns if f in eg_pairs.columns]
        
        if len(available_features) < len(feature_columns):
            missing = set(feature_columns) - set(available_features)
            logger.warning(f"Missing features: {missing}")
        
        X = eg_pairs[available_features].values
        
        # Handle missing values
        X = np.nan_to_num(X, nan=0.0)
        
        self.feature_names = available_features
        return X, available_features
    
    def train(self,
              X: np.ndarray,
              y: np.ndarray,
              feature_names: List[str] = None) -> Dict:
        """
        Train ML model.
        
        Args:
            X: Feature matrix
            y: Labels (1 = validated E-G pair, 0 = non-validated)
            feature_names: Names of features
            
        Returns:
            Training metrics
        """
        logger.info(f"Training {self.config.model_type} model...")
        logger.info(f"  Samples: {len(X)}, Features: {X.shape[1]}")
        logger.info(f"  Positive: {y.sum()}, Negative: {(1-y).sum()}")
        
        if feature_names:
            self.feature_names = feature_names
        
        # Scale features
        X_scaled = self.scaler.fit_transform(X)
        
        # Create and train model
        self.model = self._create_model()
        
        # Cross-validation
        cv = StratifiedKFold(n_splits=self.config.n_cv_folds, shuffle=True, 
                            random_state=self.config.random_state)
        
        cv_scores = cross_val_score(self.model, X_scaled, y, cv=cv, scoring='roc_auc')
        
        # Train final model on all data
        self.model.fit(X_scaled, y)
        
        # Get feature importance
        if hasattr(self.model, 'feature_importances_'):
            importance = self.model.feature_importances_
            self.feature_importance = dict(zip(self.feature_names, importance))
        
        metrics = {
            'cv_auc_mean': cv_scores.mean(),
            'cv_auc_std': cv_scores.std(),
            'n_samples': len(X),
            'n_features': X.shape[1],
            'n_positive': int(y.sum()),
            'n_negative': int((1-y).sum())
        }
        
        logger.info(f"  CV AUC: {metrics['cv_auc_mean']:.3f} ± {metrics['cv_auc_std']:.3f}")
        
        return metrics
    
    def predict(self, X: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        Predict probabilities for E-G pairs.
        
        Args:
            X: Feature matrix
            
        Returns:
            Predicted probabilities and binary predictions
        """
        if self.model is None:
            raise ValueError("Model not trained. Call train() first.")
        
        X_scaled = self.scaler.transform(X)
        
        probabilities = self.model.predict_proba(X_scaled)[:, 1]
        predictions = self.model.predict(X_scaled)
        
        return probabilities, predictions
    
    def get_feature_importance(self) -> pd.DataFrame:
        """
        Get feature importance ranking.
        
        Returns:
            DataFrame with feature importance
        """
        if not self.feature_importance:
            raise ValueError("No feature importance available. Train model first.")
        
        importance_df = pd.DataFrame([
            {'feature': k, 'importance': v}
            for k, v in self.feature_importance.items()
        ]).sort_values('importance', ascending=False)
        
        return importance_df
    
    def evaluate(self, X: np.ndarray, y: np.ndarray) -> Dict:
        """
        Evaluate model performance.
        
        Args:
            X: Feature matrix
            y: True labels
            
        Returns:
            Evaluation metrics
        """
        probabilities, predictions = self.predict(X)
        
        metrics = {
            'roc_auc': roc_auc_score(y, probabilities),
            'avg_precision': average_precision_score(y, probabilities),
            'accuracy': (predictions == y).mean()
        }
        
        # Precision-recall at different thresholds
        precision, recall, thresholds = precision_recall_curve(y, probabilities)
        
        # Find threshold for 70% recall (like ABC default)
        idx_70_recall = np.argmin(np.abs(recall - 0.7))
        metrics['precision_at_70_recall'] = precision[idx_70_recall]
        metrics['threshold_at_70_recall'] = thresholds[idx_70_recall] if idx_70_recall < len(thresholds) else 0.5
        
        return metrics
    
    def save_model(self, filepath: str):
        """Save trained model to file"""
        import pickle
        
        model_data = {
            'model': self.model,
            'scaler': self.scaler,
            'feature_names': self.feature_names,
            'feature_importance': self.feature_importance,
            'config': self.config
        }
        
        with open(filepath, 'wb') as f:
            pickle.dump(model_data, f)
        
        logger.info(f"Model saved to {filepath}")
    
    def load_model(self, filepath: str):
        """Load trained model from file"""
        import pickle
        
        with open(filepath, 'rb') as f:
            model_data = pickle.load(f)
        
        self.model = model_data['model']
        self.scaler = model_data['scaler']
        self.feature_names = model_data['feature_names']
        self.feature_importance = model_data['feature_importance']
        self.config = model_data['config']
        
        logger.info(f"Model loaded from {filepath}")


class FeatureEngineering:
    """
    Feature engineering for PACE ML predictions.
    
    Creates derived features from raw signals to improve prediction.
    """
    
    @staticmethod
    def create_features(eg_pairs: pd.DataFrame) -> pd.DataFrame:
        """
        Create engineered features from raw signals.
        
        Args:
            eg_pairs: DataFrame with E-G pairs and raw signals
            
        Returns:
            DataFrame with additional engineered features
        """
        df = eg_pairs.copy()
        
        # Distance-based features
        if 'distance' in df.columns:
            df['log_distance'] = np.log10(df['distance'] + 1)
            df['inverse_distance'] = 1 / (df['distance'] + 1000)
        
        # Activity features
        activity_cols = [c for c in df.columns if 'activity' in c.lower() or 
                        c in ['ATAC', 'DHS', 'H3K27ac', 'H3K4me1', 'H3K4me3']]
        
        if len(activity_cols) >= 2:
            # Geometric mean of activity signals
            activity_product = np.ones(len(df))
            for col in activity_cols:
                if col in df.columns:
                    activity_product *= np.maximum(df[col], 1e-10)
            df['activity_geometric_mean'] = np.power(activity_product, 1/len(activity_cols))
            
            # Sum of activity signals
            df['activity_sum'] = df[activity_cols].sum(axis=1)
        
        # Ratio features
        if 'H3K27ac' in df.columns and 'H3K4me1' in df.columns:
            df['H3K27ac_H3K4me1_ratio'] = df['H3K27ac'] / (df['H3K4me1'] + 0.01)
        
        if 'H3K4me3' in df.columns and 'H3K4me1' in df.columns:
            # Promoter vs enhancer ratio
            df['promoter_enhancer_ratio'] = df['H3K4me3'] / (df['H3K4me1'] + 0.01)
        
        # Contact × Activity interaction
        if 'contact' in df.columns and 'activity' in df.columns:
            df['activity_contact_product'] = df['activity'] * df['contact']
        
        # Expression interaction
        if 'expression' in df.columns:
            df['log_expression'] = np.log2(df['expression'] + 1)
            if 'activity' in df.columns:
                df['activity_expression_product'] = df['activity'] * df['log_expression']
        
        # Binary features
        if 'distance' in df.columns:
            df['is_proximal'] = (df['distance'] < 10000).astype(int)
            df['is_distal'] = (df['distance'] > 100000).astype(int)
        
        if 'methylation' in df.columns:
            df['is_methylated'] = (df['methylation'] > 0.5).astype(int)
        
        return df


class TransferLearning:
    """
    Transfer learning across species.
    
    Trains on well-annotated species (e.g., human, mouse) and
    transfers to livestock species.
    """
    
    def __init__(self, base_model: PACEMLPredictor):
        """
        Initialize transfer learning.
        
        Args:
            base_model: Pre-trained model on source species
        """
        self.base_model = base_model
        self.fine_tuned_model = None
        
    def fine_tune(self,
                  target_X: np.ndarray,
                  target_y: np.ndarray,
                  n_epochs: int = 10) -> Dict:
        """
        Fine-tune model on target species data.
        
        Uses base model predictions as additional features.
        
        Args:
            target_X: Target species feature matrix
            target_y: Target species labels
            n_epochs: Number of fine-tuning epochs
            
        Returns:
            Fine-tuning metrics
        """
        logger.info("Fine-tuning model on target species...")
        
        # Get base model predictions as features
        base_probs, _ = self.base_model.predict(target_X)
        
        # Augment features with base model predictions
        X_augmented = np.column_stack([target_X, base_probs])
        
        # Create and train fine-tuned model
        self.fine_tuned_model = PACEMLPredictor(self.base_model.config)
        
        feature_names = self.base_model.feature_names + ['base_model_prediction']
        
        metrics = self.fine_tuned_model.train(
            X_augmented, target_y, feature_names
        )
        
        return metrics
    
    def predict(self, X: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Predict using fine-tuned model"""
        if self.fine_tuned_model is None:
            logger.warning("No fine-tuned model. Using base model.")
            return self.base_model.predict(X)
        
        # Get base predictions
        base_probs, _ = self.base_model.predict(X)
        
        # Augment and predict
        X_augmented = np.column_stack([X, base_probs])
        
        return self.fine_tuned_model.predict(X_augmented)


def train_pace_ml_model(
    training_data: pd.DataFrame,
    label_column: str = 'validated',
    feature_columns: List[str] = None,
    config: MLConfig = None,
    output_path: str = None
) -> Tuple[PACEMLPredictor, Dict]:
    """
    Train PACE ML model from E-G pair data.
    
    Args:
        training_data: DataFrame with E-G pairs and labels
        label_column: Column containing labels
        feature_columns: List of feature columns (None = auto-detect)
        config: ML configuration
        output_path: Path to save trained model
        
    Returns:
        Trained model and metrics
    """
    # Auto-detect features if not specified
    if feature_columns is None:
        exclude_cols = {'chr', 'start', 'end', 'gene_id', 'gene_name', 
                       label_column, 'PACE.Score', 'ABC.Score'}
        feature_columns = [c for c in training_data.columns 
                         if c not in exclude_cols and training_data[c].dtype in ['float64', 'int64']]
    
    logger.info(f"Using features: {feature_columns}")
    
    # Feature engineering
    training_data = FeatureEngineering.create_features(training_data)
    
    # Prepare data
    predictor = PACEMLPredictor(config)
    X, features = predictor.prepare_features(training_data, feature_columns)
    y = training_data[label_column].values
    
    # Train
    metrics = predictor.train(X, y, features)
    
    # Save if requested
    if output_path:
        predictor.save_model(output_path)
    
    # Print feature importance
    importance = predictor.get_feature_importance()
    logger.info("\nFeature Importance:")
    for _, row in importance.head(10).iterrows():
        logger.info(f"  {row['feature']}: {row['importance']:.4f}")
    
    return predictor, metrics


if __name__ == "__main__":
    print("PACE Machine Learning Module")
    print("=" * 50)
    
    if not SKLEARN_AVAILABLE:
        print("scikit-learn not available. Install with:")
        print("  pip install scikit-learn")
        exit(1)
    
    # Example usage with synthetic data
    np.random.seed(42)
    n_samples = 1000
    
    # Create synthetic E-G pair data
    data = pd.DataFrame({
        'distance': np.random.exponential(100000, n_samples),
        'ATAC': np.random.rand(n_samples),
        'H3K27ac': np.random.rand(n_samples),
        'H3K4me1': np.random.rand(n_samples),
        'contact': np.random.rand(n_samples),
        'expression': np.random.exponential(10, n_samples),
        'methylation': np.random.rand(n_samples),
    })
    
    # Create labels (correlated with features)
    prob = 1 / (1 + np.exp(-(
        data['ATAC'] * 2 + 
        data['H3K27ac'] * 3 + 
        data['contact'] * 2 - 
        data['methylation'] * 1.5 - 
        np.log10(data['distance'] + 1) * 0.5
    )))
    data['validated'] = (np.random.rand(n_samples) < prob).astype(int)
    
    print(f"\nSynthetic data: {n_samples} E-G pairs")
    print(f"Validated pairs: {data['validated'].sum()}")
    
    # Train model
    feature_cols = ['distance', 'ATAC', 'H3K27ac', 'H3K4me1', 
                   'contact', 'expression', 'methylation']
    
    model, metrics = train_pace_ml_model(
        data, 
        label_column='validated',
        feature_columns=feature_cols
    )
    
    print(f"\nTraining Metrics:")
    for k, v in metrics.items():
        print(f"  {k}: {v}")
    
    print("\nDone!")
