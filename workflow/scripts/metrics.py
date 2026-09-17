#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PACE Metrics Module

Quality control metrics and plotting functions.

Author: Linyong Shen @ Northwest A&F University
"""

import os
import numpy as np
import pandas as pd
from typing import Optional, Dict, List, Tuple
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from tools import logger


def gene_keys(df):
    """Use stable identifiers, with sample and chromosome scope when present."""
    identifier = 'TargetGeneEnsemblID' if 'TargetGeneEnsemblID' in df else 'TargetGene'
    return [c for c in ['sample_id', 'chr', identifier] if c in df]


class PACEMetrics:
    """
    Calculator for PACE quality control metrics.
    """
    
    def __init__(self, predictions: pd.DataFrame, score_column: str = 'PACE.Score'):
        """
        Initialize metrics calculator.
        
        Args:
            predictions: DataFrame with predictions
            score_column: Name of score column
        """
        self.predictions = predictions
        self.score_column = score_column
        self.metrics = {}
    
    def calculate_basic_metrics(self) -> Dict:
        """
        Calculate basic prediction metrics.
        
        Returns:
            Dictionary of metrics
        """
        df = self.predictions
        
        metrics = {
            'total_predictions': len(df),
            'unique_enhancers': len(df[[c for c in ['sample_id', 'chr', 'start', 'end'] if c in df]].drop_duplicates())
                if {'chr', 'start', 'end'}.issubset(df) else 0,
            'unique_genes': len(df[gene_keys(df)].drop_duplicates())
                if {'TargetGeneEnsemblID', 'TargetGene'}.intersection(df) else 0,
        }
        
        if self.score_column in df.columns:
            metrics['mean_score'] = df[self.score_column].mean()
            metrics['median_score'] = df[self.score_column].median()
            metrics['max_score'] = df[self.score_column].max()
            metrics['score_std'] = df[self.score_column].std()
            
            # Count by threshold
            for thresh in [0.01, 0.02, 0.03, 0.05]:
                metrics[f'n_above_{thresh}'] = (df[self.score_column] >= thresh).sum()
        
        if 'distance' in df.columns:
            metrics['mean_distance'] = df['distance'].mean()
            metrics['median_distance'] = df['distance'].median()
        
        if 'class' in df.columns:
            class_counts = df['class'].value_counts().to_dict()
            for cls, count in class_counts.items():
                metrics[f'n_{cls}'] = count
        
        self.metrics.update(metrics)
        return metrics
    
    def calculate_activity_metrics(self) -> Dict:
        """
        Calculate activity-related metrics.
        
        Returns:
            Dictionary of metrics
        """
        df = self.predictions
        metrics = {}
        
        if 'activity' in df.columns:
            metrics['mean_activity'] = df['activity'].mean()
            metrics['median_activity'] = df['activity'].median()
            metrics['activity_std'] = df['activity'].std()
        
        if 'contact' in df.columns:
            metrics['mean_contact'] = df['contact'].mean()
            metrics['median_contact'] = df['contact'].median()
        
        self.metrics.update(metrics)
        return metrics
    
    def calculate_gene_metrics(self) -> Dict:
        """
        Calculate per-gene metrics.
        
        Returns:
            Dictionary of metrics
        """
        df = self.predictions
        metrics = {}
        
        if {'TargetGeneEnsemblID', 'TargetGene'}.intersection(df):
            enhancers_per_gene = df.groupby(gene_keys(df)).size()
            metrics['mean_enhancers_per_gene'] = enhancers_per_gene.mean()
            metrics['median_enhancers_per_gene'] = enhancers_per_gene.median()
            metrics['max_enhancers_per_gene'] = enhancers_per_gene.max()
            
            # A descriptive score cutoff is not a confidence calibration.
            if self.score_column in df.columns:
                selected = df[df[self.score_column] >= 0.02]
                metrics['genes_with_score_ge_0_02'] = len(selected[gene_keys(df)].drop_duplicates())
        
        self.metrics.update(metrics)
        return metrics
    
    def calculate_all_metrics(self) -> Dict:
        """
        Calculate all metrics.
        
        Returns:
            Dictionary of all metrics
        """
        self.calculate_basic_metrics()
        self.calculate_activity_metrics()
        self.calculate_gene_metrics()
        return self.metrics
    
    def to_dataframe(self) -> pd.DataFrame:
        """
        Convert metrics to DataFrame.
        
        Returns:
            DataFrame with metrics
        """
        return pd.DataFrame([self.metrics])
    
    def save(self, output_file: str) -> None:
        """
        Save metrics to file.
        
        Args:
            output_file: Output file path
        """
        df = self.to_dataframe()
        df.to_csv(output_file, sep='\t', index=False)
        logger.info(f"Saved metrics to {output_file}")


class PACEPlotter:
    """
    Plotter for PACE quality control visualizations.
    """
    
    def __init__(self, predictions: pd.DataFrame, score_column: str = 'PACE.Score'):
        """
        Initialize plotter.
        
        Args:
            predictions: DataFrame with predictions
            score_column: Name of score column
        """
        self.predictions = predictions
        self.score_column = score_column
        
        # Set style
        plt.style.use('default')
        plt.rcParams['figure.dpi'] = 150
    
    def plot_score_distribution(self, ax: Optional[plt.Axes] = None) -> plt.Axes:
        """
        Plot score distribution.
        
        Args:
            ax: Optional matplotlib axes
        
        Returns:
            Matplotlib axes
        """
        if ax is None:
            fig, ax = plt.subplots(figsize=(6, 4))
        
        scores = self.predictions[self.score_column].dropna()
        
        ax.hist(scores, bins=50, edgecolor='black', alpha=0.7)
        ax.axvline(0.02, color='red', linestyle='--', label='Threshold (0.02)')
        ax.set_xlabel('Score')
        ax.set_ylabel('Count')
        ax.set_title('Score Distribution')
        ax.legend()
        
        return ax
    
    def plot_distance_distribution(self, ax: Optional[plt.Axes] = None) -> plt.Axes:
        """
        Plot distance distribution.
        
        Args:
            ax: Optional matplotlib axes
        
        Returns:
            Matplotlib axes
        """
        if ax is None:
            fig, ax = plt.subplots(figsize=(6, 4))
        
        distances = self.predictions['distance'] / 1000  # Convert to kb
        
        ax.hist(distances, bins=50, edgecolor='black', alpha=0.7)
        ax.set_xlabel('Distance to TSS (kb)')
        ax.set_ylabel('Count')
        ax.set_title('Distance Distribution')
        
        return ax
    
    def plot_score_vs_distance(self, ax: Optional[plt.Axes] = None) -> plt.Axes:
        """
        Plot score vs distance.
        
        Args:
            ax: Optional matplotlib axes
        
        Returns:
            Matplotlib axes
        """
        if ax is None:
            fig, ax = plt.subplots(figsize=(6, 4))
        
        df = self.predictions.sample(min(5000, len(self.predictions)), random_state=42)
        
        ax.scatter(
            df['distance'] / 1000,
            df[self.score_column],
            alpha=0.3,
            s=5
        )
        ax.set_xlabel('Distance to TSS (kb)')
        ax.set_ylabel('Score')
        ax.set_title('Score vs Distance')
        ax.axhline(0.02, color='red', linestyle='--', alpha=0.5)
        
        return ax
    
    def plot_enhancer_class(self, ax: Optional[plt.Axes] = None) -> plt.Axes:
        """
        Plot enhancer class distribution.
        
        Args:
            ax: Optional matplotlib axes
        
        Returns:
            Matplotlib axes
        """
        if ax is None:
            fig, ax = plt.subplots(figsize=(6, 4))
        
        if 'class' in self.predictions.columns:
            class_counts = self.predictions['class'].value_counts()
            
            ax.bar(class_counts.index, class_counts.values, edgecolor='black', alpha=0.7)
            ax.set_xlabel('Enhancer Class')
            ax.set_ylabel('Count')
            ax.set_title('Enhancer Class Distribution')
        
        return ax
    
    def plot_activity_distribution(self, ax: Optional[plt.Axes] = None) -> plt.Axes:
        """
        Plot activity distribution.
        
        Args:
            ax: Optional matplotlib axes
        
        Returns:
            Matplotlib axes
        """
        if ax is None:
            fig, ax = plt.subplots(figsize=(6, 4))
        
        if 'activity' in self.predictions.columns:
            activities = self.predictions['activity']
            activities = activities[activities > 0]
            
            ax.hist(np.log10(activities + 1e-10), bins=50, edgecolor='black', alpha=0.7)
            ax.set_xlabel('log10(Activity)')
            ax.set_ylabel('Count')
            ax.set_title('Activity Distribution')
        
        return ax
    
    def plot_enhancers_per_gene(self, ax: Optional[plt.Axes] = None) -> plt.Axes:
        """
        Plot enhancers per gene distribution.
        
        Args:
            ax: Optional matplotlib axes
        
        Returns:
            Matplotlib axes
        """
        if ax is None:
            fig, ax = plt.subplots(figsize=(6, 4))
        
        enhancers_per_gene = self.predictions.groupby(gene_keys(self.predictions)).size()
        
        ax.hist(enhancers_per_gene, bins=50, edgecolor='black', alpha=0.7)
        ax.set_xlabel('Number of Enhancers')
        ax.set_ylabel('Number of Genes')
        ax.set_title('Enhancers per Gene')
        
        return ax
    
    def create_qc_report(self, output_file: str) -> None:
        """
        Create comprehensive QC report with multiple plots.
        
        Args:
            output_file: Output PDF file
        """
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        
        self.plot_score_distribution(axes[0, 0])
        self.plot_distance_distribution(axes[0, 1])
        self.plot_score_vs_distance(axes[0, 2])
        self.plot_enhancer_class(axes[1, 0])
        self.plot_activity_distribution(axes[1, 1])
        self.plot_enhancers_per_gene(axes[1, 2])
        
        plt.tight_layout()
        plt.savefig(output_file, bbox_inches='tight')
        plt.close()
        
        logger.info(f"Saved QC report to {output_file}")


def generate_metrics(predictions_file: str,
                    output_dir: str,
                    score_column: str = 'PACE.Score',
                    sample_name: str = 'sample') -> Tuple[str, str]:
    """
    Generate metrics and plots for predictions.
    
    Args:
        predictions_file: Path to predictions file
        output_dir: Output directory
        score_column: Score column name
        sample_name: Sample name for output files
    
    Returns:
        Tuple of (metrics_file, plot_file)
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Load predictions
    if predictions_file.endswith('.gz'):
        df = pd.read_csv(predictions_file, sep='\t', compression='gzip')
    else:
        df = pd.read_csv(predictions_file, sep='\t')
    
    # Calculate metrics
    metrics = PACEMetrics(df, score_column)
    metrics.calculate_all_metrics()
    
    metrics_file = os.path.join(output_dir, f'QCSummary_{sample_name}.tsv')
    metrics.save(metrics_file)
    
    # Create plots
    plotter = PACEPlotter(df, score_column)
    plot_file = os.path.join(output_dir, f'QCPlots_{sample_name}.pdf')
    plotter.create_qc_report(plot_file)
    
    return metrics_file, plot_file


def compare_samples(predictions_files: Dict[str, str],
                   output_file: str,
                   score_column: str = 'PACE.Score') -> pd.DataFrame:
    """
    Compare metrics across multiple samples.
    
    Args:
        predictions_files: Dictionary of sample_name -> predictions_file
        output_file: Output comparison file
        score_column: Score column name
    
    Returns:
        DataFrame with comparison
    """
    all_metrics = []
    
    for sample_name, pred_file in predictions_files.items():
        if pred_file.endswith('.gz'):
            df = pd.read_csv(pred_file, sep='\t', compression='gzip')
        else:
            df = pd.read_csv(pred_file, sep='\t')
        
        metrics = PACEMetrics(df, score_column)
        metrics.calculate_all_metrics()
        metrics.metrics['sample'] = sample_name
        all_metrics.append(metrics.metrics)
    
    comparison = pd.DataFrame(all_metrics)
    
    # Reorder columns
    cols = ['sample'] + [c for c in comparison.columns if c != 'sample']
    comparison = comparison[cols]
    
    comparison.to_csv(output_file, sep='\t', index=False)
    logger.info(f"Saved comparison to {output_file}")
    
    return comparison
