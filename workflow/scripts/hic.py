#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PACE Hi-C Processing Module

Functions for processing Hi-C data and estimating 3D chromatin contacts.

Author: Linyong Shen @ Northwest A&F University
"""

import os
import numpy as np
import pandas as pd
from typing import Optional, Dict, Tuple, Union
import tempfile

from tools import power_law_contact, logger


class HiCProcessor:
    """
    Processor for Hi-C data and contact estimation.
    """
    
    def __init__(self,
                 resolution: int = 5000,
                 hic_gamma: float = 1.024238616787792,
                 hic_scale: float = 5.9594510043736655,
                 hic_pseudocount: int = 5000):
        """
        Initialize Hi-C processor.
        
        Args:
            resolution: Hi-C resolution in bp
            hic_gamma: Power-law exponent for distance decay
            hic_scale: Scale factor for power-law
            hic_pseudocount: Pseudocount for distance
        """
        self.resolution = resolution
        self.hic_gamma = hic_gamma
        self.hic_scale = hic_scale
        self.hic_pseudocount = hic_pseudocount
        
        self.contact_matrix = {}  # chromosome -> sparse contact matrix
        self.hic_loaded = False
    
    def load_hic_file(self, hic_file: str, hic_type: str = 'hic') -> None:
        """
        Load Hi-C data from file.
        
        Args:
            hic_file: Path to Hi-C file
            hic_type: Type of file (hic, bedpe, cool)
        """
        if hic_type == 'hic':
            self._load_hic_format(hic_file)
        elif hic_type == 'bedpe':
            self._load_bedpe_format(hic_file)
        elif hic_type == 'cool':
            self._load_cool_format(hic_file)
        else:
            logger.warning(f"Unknown Hi-C format: {hic_type}, using power-law")
            return
        
        self.hic_loaded = True
        logger.info(f"Loaded Hi-C data from {hic_file}")
    
    def _load_hic_format(self, hic_file: str) -> None:
        """Load .hic format file using hic-straw."""
        try:
            import hicstraw
            
            hic = hicstraw.HiCFile(hic_file)
            
            for chrom in hic.getChromosomes():
                if chrom.name in ['All', 'ALL']:
                    continue
                
                try:
                    matrix = hic.getMatrixZoomData(
                        chrom.name, chrom.name,
                        'observed', 'KR', 'BP', self.resolution
                    )
                    
                    # Get contact records
                    records = matrix.getRecords(0, chrom.length, 0, chrom.length)
                    
                    # Store as sparse dict
                    contacts = {}
                    for record in records:
                        bin1 = record.binX // self.resolution
                        bin2 = record.binY // self.resolution
                        contacts[(bin1, bin2)] = record.counts
                        contacts[(bin2, bin1)] = record.counts
                    
                    self.contact_matrix[chrom.name] = contacts
                    
                except Exception as e:
                    logger.warning(f"Could not load {chrom.name}: {e}")
        
        except ImportError:
            logger.warning("hicstraw not installed, using power-law estimation")
    
    def _load_bedpe_format(self, bedpe_file: str) -> None:
        """Load BEDPE format Hi-C file."""
        df = pd.read_csv(bedpe_file, sep='\t', header=None,
                        names=['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'score'])
        
        # Only intra-chromosomal contacts
        df = df[df['chr1'] == df['chr2']]
        
        for chrom in df['chr1'].unique():
            chrom_df = df[df['chr1'] == chrom]
            contacts = {}
            
            for _, row in chrom_df.iterrows():
                bin1 = (row['start1'] + row['end1']) // 2 // self.resolution
                bin2 = (row['start2'] + row['end2']) // 2 // self.resolution
                contacts[(bin1, bin2)] = row['score']
                contacts[(bin2, bin1)] = row['score']
            
            self.contact_matrix[chrom] = contacts
    
    def _load_cool_format(self, cool_file: str) -> None:
        """Load cooler format Hi-C file."""
        try:
            import cooler
            
            clr = cooler.Cooler(cool_file)
            
            for chrom in clr.chromnames:
                try:
                    matrix = clr.matrix(balance=True).fetch(chrom)
                    
                    # Convert to sparse dict
                    contacts = {}
                    rows, cols = np.where(~np.isnan(matrix) & (matrix > 0))
                    for r, c in zip(rows, cols):
                        contacts[(r, c)] = matrix[r, c]
                    
                    self.contact_matrix[chrom] = contacts
                    
                except Exception as e:
                    logger.warning(f"Could not load {chrom}: {e}")
        
        except ImportError:
            logger.warning("cooler not installed, using power-law estimation")
    
    def get_contact(self,
                   chrom: str,
                   pos1: int,
                   pos2: int,
                   use_powerlaw_fallback: bool = True) -> float:
        """
        Get contact frequency between two positions.
        
        Args:
            chrom: Chromosome name
            pos1: First position
            pos2: Second position
            use_powerlaw_fallback: Use power-law if no Hi-C data
        
        Returns:
            Contact frequency
        """
        if self.hic_loaded and chrom in self.contact_matrix:
            bin1 = pos1 // self.resolution
            bin2 = pos2 // self.resolution
            
            contacts = self.contact_matrix[chrom]
            contact = contacts.get((bin1, bin2), 0)
            
            if contact > 0:
                return contact
        
        # Fall back to power-law
        if use_powerlaw_fallback:
            distance = abs(pos1 - pos2)
            return self.estimate_contact_powerlaw(distance)
        
        return 0.0
    
    def estimate_contact_powerlaw(self, distance: Union[int, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Estimate contact frequency using power-law.
        
        Args:
            distance: Distance(s) in bp
        
        Returns:
            Estimated contact frequency
        """
        return power_law_contact(
            distance,
            hic_gamma=self.hic_gamma,
            hic_scale=self.hic_scale,
            hic_pseudocount=self.hic_pseudocount
        )
    
    def scale_hic_by_powerlaw(self,
                             predictions: pd.DataFrame,
                             distance_column: str = 'distance',
                             contact_column: str = 'hic_contact') -> pd.DataFrame:
        """
        Scale Hi-C contacts to match power-law at long distances.
        
        Args:
            predictions: DataFrame with predictions
            distance_column: Column with distance values
            contact_column: Column with Hi-C contact values
        
        Returns:
            DataFrame with scaled contacts
        """
        df = predictions.copy()
        
        # Get expected contacts from power-law
        df['powerlaw_contact'] = self.estimate_contact_powerlaw(df[distance_column].values)
        
        # Calculate scaling factor for long-range contacts (>100kb)
        long_range = df[df[distance_column] > 100000]
        
        if len(long_range) > 0:
            scale_factor = (
                long_range['powerlaw_contact'].mean() / 
                (long_range[contact_column].mean() + 1e-10)
            )
            
            df[contact_column] = df[contact_column] * scale_factor
        
        df.drop('powerlaw_contact', axis=1, inplace=True)
        
        return df
    
    def get_contact_matrix(self, chrom: str) -> Optional[Dict]:
        """
        Get contact matrix for a chromosome.
        
        Args:
            chrom: Chromosome name
        
        Returns:
            Sparse contact dictionary or None
        """
        return self.contact_matrix.get(chrom)


def load_hic_avg(avg_file: str) -> Dict[int, float]:
    """
    Load average Hi-C contact by distance file.
    
    Args:
        avg_file: Path to average Hi-C file
    
    Returns:
        Dictionary of distance -> contact
    """
    avg_contacts = {}
    
    with open(avg_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                try:
                    distance = int(parts[0])
                    contact = float(parts[1])
                    avg_contacts[distance] = contact
                except ValueError:
                    continue
    
    return avg_contacts


def interpolate_contact(distance: int,
                       avg_contacts: Dict[int, float],
                       resolution: int = 5000) -> float:
    """
    Interpolate contact frequency from average Hi-C.
    
    Args:
        distance: Distance in bp
        avg_contacts: Average contacts by distance
        resolution: Hi-C resolution
    
    Returns:
        Interpolated contact
    """
    # Round to nearest resolution
    dist_bin = (distance // resolution) * resolution
    
    if dist_bin in avg_contacts:
        return avg_contacts[dist_bin]
    
    # Find nearest distances
    distances = sorted(avg_contacts.keys())
    
    if not distances:
        return 0.0
    
    if dist_bin <= distances[0]:
        return avg_contacts[distances[0]]
    
    if dist_bin >= distances[-1]:
        return avg_contacts[distances[-1]]
    
    # Linear interpolation
    for i, d in enumerate(distances[:-1]):
        if d <= dist_bin < distances[i + 1]:
            d1, d2 = d, distances[i + 1]
            c1, c2 = avg_contacts[d1], avg_contacts[d2]
            return c1 + (c2 - c1) * (dist_bin - d1) / (d2 - d1)
    
    return 0.0


class ContactEstimator:
    """
    Unified contact estimator supporting multiple methods.
    """
    
    def __init__(self,
                 method: str = 'powerlaw',
                 hic_file: Optional[str] = None,
                 hic_type: str = 'hic',
                 resolution: int = 5000,
                 hic_gamma: float = 1.024238616787792,
                 hic_scale: float = 5.9594510043736655):
        """
        Initialize contact estimator.
        
        Args:
            method: Estimation method (powerlaw, hic, avg)
            hic_file: Path to Hi-C file
            hic_type: Type of Hi-C file
            resolution: Hi-C resolution
            hic_gamma: Power-law exponent
            hic_scale: Power-law scale
        """
        self.method = method
        self.resolution = resolution
        self.hic_gamma = hic_gamma
        self.hic_scale = hic_scale
        
        self.hic_processor = None
        self.avg_contacts = None
        
        if method == 'hic' and hic_file:
            self.hic_processor = HiCProcessor(
                resolution=resolution,
                hic_gamma=hic_gamma,
                hic_scale=hic_scale
            )
            self.hic_processor.load_hic_file(hic_file, hic_type)
        
        elif method == 'avg' and hic_file:
            self.avg_contacts = load_hic_avg(hic_file)
    
    def estimate(self,
                chrom: str,
                pos1: int,
                pos2: int) -> float:
        """
        Estimate contact between two positions.
        
        Args:
            chrom: Chromosome
            pos1: First position
            pos2: Second position
        
        Returns:
            Estimated contact
        """
        distance = abs(pos1 - pos2)
        
        if self.method == 'powerlaw':
            return power_law_contact(
                distance,
                hic_gamma=self.hic_gamma,
                hic_scale=self.hic_scale
            )
        
        elif self.method == 'hic' and self.hic_processor:
            return self.hic_processor.get_contact(chrom, pos1, pos2)
        
        elif self.method == 'avg' and self.avg_contacts:
            return interpolate_contact(distance, self.avg_contacts, self.resolution)
        
        else:
            # Default to power-law
            return power_law_contact(distance, hic_gamma=self.hic_gamma, hic_scale=self.hic_scale)
    
    def estimate_batch(self,
                      df: pd.DataFrame,
                      chrom_col: str = 'chr',
                      pos1_col: str = 'enhancer_mid',
                      pos2_col: str = 'TargetGeneTSS') -> np.ndarray:
        """
        Estimate contacts for batch of pairs.
        
        Args:
            df: DataFrame with enhancer-gene pairs
            chrom_col: Chromosome column
            pos1_col: Enhancer position column
            pos2_col: Gene TSS column
        
        Returns:
            Array of contact estimates
        """
        contacts = []
        
        for _, row in df.iterrows():
            contact = self.estimate(
                row[chrom_col],
                row[pos1_col],
                row[pos2_col]
            )
            contacts.append(contact)
        
        return np.array(contacts)
