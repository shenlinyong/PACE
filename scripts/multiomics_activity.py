#!/usr/bin/env python3
"""
PACE: Multi-omics Activity Calculator

This module provides flexible activity score calculation by integrating
multiple epigenomic signals with configurable weights and aggregation
methods. It is the reference implementation of the PACE multi-layer
enhancer activity model.

All signal quantification is performed directly from the input files
(bigWig via pyBigWig, or BED/tagAlign/BAM via bedtools coverage). No part
of this module fabricates or simulates data: given the same inputs and the
same configuration it always returns the same activity scores.

Supported data types
--------------------
- Chromatin accessibility: ATAC-seq, DNase-seq
- Histone modifications: H3K27ac, H3K4me1, H3K4me3, H3K9ac, H3K36me3,
  H3K27me3, H3K9me3
- Transcription factors: TF ChIP-seq (activator OR repressor; see
  ``is_inhibitory`` flag)
- DNA methylation: WGBS / RRBS (inhibitory)
- Gene expression: RNA-seq

Author: Linyong Shen @ Northwest A&F University
"""

import os
import gzip
import tempfile
import subprocess
import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Union
from dataclasses import dataclass, field
from enum import Enum
import logging

try:
    import pyBigWig
    HAS_PYBIGWIG = True
except ImportError:  # pragma: no cover - optional dependency
    pyBigWig = None
    HAS_PYBIGWIG = False

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class AggregationMethod(Enum):
    """Methods for aggregating multiple signals."""
    GEOMETRIC_MEAN = "geometric_mean"        # Original ABC method
    ARITHMETIC_MEAN = "arithmetic_mean"
    WEIGHTED_SUM = "weighted_sum"
    WEIGHTED_GEOMETRIC = "weighted_geometric"
    MAX = "max"
    PRODUCT = "product"


class SignalType(Enum):
    """Types of epigenomic signals."""
    ACCESSIBILITY = "accessibility"        # ATAC, DNase
    ACTIVE_ENHANCER = "active_enhancer"    # H3K27ac
    ENHANCER_MARK = "enhancer_mark"        # H3K4me1
    PROMOTER_MARK = "promoter_mark"        # H3K4me3
    TRANSCRIPTION = "transcription"        # H3K36me3
    REPRESSIVE = "repressive"              # H3K27me3, H3K9me3
    TF_BINDING = "tf_binding"              # TF ChIP-seq
    METHYLATION = "methylation"            # DNA methylation (inhibitory)
    EXPRESSION = "expression"              # RNA-seq, eRNA


# Signals that are inhibitory by default (down-weight enhancer activity).
DEFAULT_INHIBITORY_TYPES = {SignalType.REPRESSIVE, SignalType.METHYLATION}


@dataclass
class SignalConfig:
    """Configuration for a single signal type."""
    name: str
    signal_type: SignalType
    file_path: str
    weight: float = 1.0
    is_inhibitory: bool = False   # True for repressive marks, methylation, TF repressors
    normalize: str = "rpkm"       # rpkm | quantile | minmax | none
    log_transform: bool = False
    pseudocount: float = 1.0
    stat: str = "mean"            # bigWig summary statistic


@dataclass
class ActivityConfig:
    """Configuration for activity calculation."""
    signals: List[SignalConfig] = field(default_factory=list)
    aggregation_method: AggregationMethod = AggregationMethod.GEOMETRIC_MEAN
    min_signals_required: int = 1
    accessibility_required: bool = True
    output_individual_scores: bool = True


# --------------------------------------------------------------------------- #
# Real signal quantification (no simulated data)
# --------------------------------------------------------------------------- #

def _run(cmd: str) -> None:
    """Run a shell command, raising on failure."""
    subprocess.run(cmd, shell=True, check=True)


def _is_bigwig(path: str) -> bool:
    return path.lower().endswith((".bw", ".bigwig"))


def quantify_signal_over_regions(regions: pd.DataFrame,
                                 file_path: str,
                                 stat: str = "mean",
                                 region_size_norm: bool = True) -> np.ndarray:
    """Quantify a genomic signal over a set of regions.

    Reads coverage directly from the data file. bigWig files are read with
    pyBigWig; BED/tagAlign/BAM files are quantified with ``bedtools coverage``.
    The result is deterministic for a given input.

    Args:
        regions: DataFrame with ``chr``, ``start``, ``end`` (and optionally
            ``name``) columns.
        file_path: Path to the signal file.
        stat: Statistic for bigWig quantification (mean/max/sum).
        region_size_norm: For read-count based files, divide by region length
            (kb) so the result is read density (reads/kb), making regions of
            different sizes comparable.

    Returns:
        np.ndarray of per-region signal values aligned with ``regions``.

    Raises:
        FileNotFoundError: if the signal file does not exist.
        ImportError: if a bigWig is supplied but pyBigWig is unavailable.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Signal file not found: {file_path}")

    n = len(regions)
    if n == 0:
        return np.array([], dtype=float)

    if _is_bigwig(file_path):
        if not HAS_PYBIGWIG:
            raise ImportError(
                "pyBigWig is required to read bigWig files. "
                "Install with: pip install pyBigWig"
            )
        bw = pyBigWig.open(file_path)
        values = np.zeros(n, dtype=float)
        for i, (_, row) in enumerate(regions.iterrows()):
            try:
                v = bw.stats(str(row["chr"]), int(row["start"]),
                             int(row["end"]), type=stat)[0]
                values[i] = v if v is not None else 0.0
            except (RuntimeError, ValueError):
                values[i] = 0.0
        bw.close()
        return values

    # BED / tagAlign / BAM -> bedtools coverage -counts.
    # Use a unique internal id (not the user 'name') for the round-trip so
    # regions with duplicate names are not collapsed when mapping counts back.
    work = regions.copy().reset_index(drop=True)
    work["_pace_id"] = [f"_r{i}" for i in range(len(work))]

    with tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False) as tmp_in:
        work[["chr", "start", "end", "_pace_id"]].to_csv(
            tmp_in, sep="\t", header=False, index=False)
        in_name = tmp_in.name
    out_name = in_name + ".cov"
    try:
        _run(f"bedtools coverage -a {in_name} -b {file_path} -counts "
             f"> {out_name}")
        cov = pd.read_csv(out_name, sep="\t", header=None,
                          names=["chr", "start", "end", "_pace_id", "count"])
    finally:
        for p in (in_name, out_name):
            if os.path.exists(p):
                os.unlink(p)

    count_map = dict(zip(cov["_pace_id"], cov["count"]))
    counts = work["_pace_id"].map(count_map).fillna(0).to_numpy(dtype=float)

    if region_size_norm:
        widths_kb = (work["end"] - work["start"]).to_numpy(dtype=float) / 1000.0
        widths_kb[widths_kb <= 0] = 1e-6
        counts = counts / widths_kb
    return counts


def _normalize(values: np.ndarray, method: str) -> np.ndarray:
    """Normalize a signal vector with the requested method."""
    values = np.asarray(values, dtype=float)
    if method == "none" or len(values) == 0:
        return values
    if method == "rpkm":
        total = values.sum()
        return values * 1e6 / total if total > 0 else values
    if method == "minmax":
        lo, hi = values.min(), values.max()
        return (values - lo) / (hi - lo) if hi > lo else np.zeros_like(values)
    if method == "quantile":
        # Map to uniform quantiles in (0, 1] – robust, scale-free.
        order = values.argsort()
        ranks = np.empty_like(order, dtype=float)
        ranks[order] = np.arange(1, len(values) + 1)
        return ranks / len(values)
    return values


class MultiOmicsActivityCalculator:
    """Calculate enhancer activity from multiple epigenomic signals.

    This is the core of PACE, allowing flexible integration of diverse data
    types with configurable weights and aggregation methods.

    Example
    -------
    >>> calc = MultiOmicsActivityCalculator(
    ...     method=AggregationMethod.WEIGHTED_GEOMETRIC)
    >>> calc.add_signal("ATAC", SignalType.ACCESSIBILITY, "atac.bw", weight=1.5)
    >>> calc.add_signal("H3K27ac", SignalType.ACTIVE_ENHANCER,
    ...                 "h3k27ac.bw", weight=1.0)
    >>> activity = calc.calculate_activity(regions)
    """

    def __init__(self,
                 method: AggregationMethod = AggregationMethod.GEOMETRIC_MEAN,
                 accessibility_required: bool = True):
        self.signals: List[SignalConfig] = []
        self.method = method
        self.accessibility_required = accessibility_required
        self._signal_data: Dict[str, np.ndarray] = {}

    def add_signal(self,
                   name: str,
                   signal_type: SignalType,
                   file_path: str,
                   weight: float = 1.0,
                   inhibitory: Optional[bool] = None,
                   normalize: str = "rpkm",
                   log_transform: bool = False,
                   stat: str = "mean") -> None:
        """Add a signal to the activity calculation.

        Args:
            name: Signal name (e.g. ``"H3K27ac"``).
            signal_type: Type of signal.
            file_path: Path to signal file (bigWig / BED / tagAlign / BAM).
            weight: Weight for this signal.
            inhibitory: Whether this is an inhibitory signal. If ``None`` the
                default for the signal type is used (repressive marks and
                methylation are inhibitory by default). Transcription factors
                may be flagged inhibitory to model repressors.
            normalize: Normalization method (rpkm/quantile/minmax/none).
            log_transform: Whether to log2-transform the signal.
            stat: bigWig summary statistic.
        """
        if inhibitory is None:
            inhibitory = signal_type in DEFAULT_INHIBITORY_TYPES
        config = SignalConfig(
            name=name, signal_type=signal_type, file_path=file_path,
            weight=weight, is_inhibitory=inhibitory, normalize=normalize,
            log_transform=log_transform, stat=stat,
        )
        self.signals.append(config)
        logger.info("Added signal: %s (type=%s, weight=%.2f, inhibitory=%s)",
                    name, signal_type.value, weight, inhibitory)

    def _load_signal(self, config: SignalConfig,
                     regions: pd.DataFrame) -> np.ndarray:
        """Load and process a real signal for the given regions."""
        logger.info("Quantifying signal %s from %s",
                    config.name, config.file_path)
        values = quantify_signal_over_regions(
            regions, config.file_path, stat=config.stat)

        if config.log_transform:
            values = np.log2(values + config.pseudocount)

        values = _normalize(values, config.normalize)
        return values

    def _aggregate_signals(self,
                           activating: Dict[str, np.ndarray],
                           inhibitory: Dict[str, np.ndarray],
                           weights_act: Dict[str, float],
                           weights_inh: Dict[str, float]) -> np.ndarray:
        """Aggregate signals into a single activity score.

        GEOMETRIC_MEAN (ABC):      Activity = (prod s_i) ** (1/n)
        WEIGHTED_GEOMETRIC (PACE): Activity = (prod s_i ** w_i) ** (1/sum w_i)
        WEIGHTED_SUM:              Activity = sum (w_i * s_i)

        Inhibitory signals: Final = Activating * (1 - clip(sum w_k s_k, 0, 1))
        """
        n_regions = len(next(iter(activating.values())))

        if self.method == AggregationMethod.GEOMETRIC_MEAN:
            log_sum = np.zeros(n_regions)
            for values in activating.values():
                log_sum += np.log(np.maximum(values, 1e-10))
            activating_score = np.exp(log_sum / len(activating))

        elif self.method == AggregationMethod.WEIGHTED_GEOMETRIC:
            # Normalized weighted geometric mean so the activity scale is
            # independent of the absolute weight magnitudes (consistent with
            # neighborhoods.py).
            log_sum = np.zeros(n_regions)
            weight_sum = 0.0
            for name, values in activating.items():
                w = weights_act[name]
                log_sum += w * np.log(np.maximum(values, 1e-10))
                weight_sum += w
            weight_sum = weight_sum if weight_sum > 0 else 1.0
            activating_score = np.exp(log_sum / weight_sum)

        elif self.method == AggregationMethod.WEIGHTED_SUM:
            activating_score = np.zeros(n_regions)
            for name, values in activating.items():
                activating_score += weights_act[name] * values

        elif self.method == AggregationMethod.ARITHMETIC_MEAN:
            activating_score = np.mean(list(activating.values()), axis=0)

        elif self.method == AggregationMethod.MAX:
            activating_score = np.max(list(activating.values()), axis=0)

        elif self.method == AggregationMethod.PRODUCT:
            activating_score = np.ones(n_regions)
            for values in activating.values():
                activating_score *= values
        else:
            raise ValueError(f"Unknown aggregation method: {self.method}")

        if inhibitory:
            inhibitory_score = np.zeros(n_regions)
            for name, values in inhibitory.items():
                # Inhibitory signals are scaled to [0, 1] so the weight is the
                # maximum fractional reduction of activity.
                v = _normalize(values, "minmax")
                inhibitory_score += weights_inh[name] * v
            inhibitory_score = np.clip(inhibitory_score, 0, 1)
            final_activity = activating_score * (1 - inhibitory_score)
        else:
            final_activity = activating_score

        return final_activity

    def calculate_activity(self,
                           regions: pd.DataFrame,
                           return_components: bool = False
                           ) -> Union[np.ndarray, Dict]:
        """Calculate activity scores for the given regions."""
        if self.accessibility_required:
            has_acc = any(s.signal_type == SignalType.ACCESSIBILITY
                          for s in self.signals)
            if not has_acc:
                raise ValueError(
                    "Accessibility signal is required but not provided")

        if not self.signals:
            raise ValueError("No signals added to the calculator")

        activating: Dict[str, np.ndarray] = {}
        inhibitory: Dict[str, np.ndarray] = {}
        weights_act: Dict[str, float] = {}
        weights_inh: Dict[str, float] = {}

        for config in self.signals:
            values = self._load_signal(config, regions)
            self._signal_data[config.name] = values
            if config.is_inhibitory:
                inhibitory[config.name] = values
                weights_inh[config.name] = config.weight
            else:
                activating[config.name] = values
                weights_act[config.name] = config.weight

        if not activating:
            raise ValueError("At least one activating signal is required")

        activity = self._aggregate_signals(
            activating, inhibitory, weights_act, weights_inh)

        if return_components:
            return {
                "activity": activity,
                "activating_signals": activating,
                "inhibitory_signals": inhibitory,
                "weights": {"activating": weights_act,
                            "inhibitory": weights_inh},
            }
        return activity


class ExpressionFilter:
    """Filter / weight enhancer-gene predictions by gene expression."""

    def __init__(self, expression_file: str, min_expression: float = 1.0,
                 expression_column: str = "TPM",
                 gene_id_column: Optional[str] = None):
        self.expression_file = expression_file
        self.min_expression = min_expression
        self.expression_column = expression_column
        self.gene_id_column = gene_id_column
        self._expression_data: Optional[pd.DataFrame] = None

    def load_expression(self) -> pd.DataFrame:
        logger.info("Loading expression data from %s", self.expression_file)
        self._expression_data = pd.read_csv(self.expression_file, sep="\t")
        return self._expression_data

    def _id_col(self) -> str:
        if self.gene_id_column:
            return self.gene_id_column
        for c in ("gene_id", "gene_name", "GeneID", "Gene", "TargetGene"):
            if c in self._expression_data.columns:
                return c
        return self._expression_data.columns[0]

    def get_expressed_genes(self) -> List[str]:
        if self._expression_data is None:
            self.load_expression()
        expressed = self._expression_data[
            self._expression_data[self.expression_column] >= self.min_expression
        ]
        return expressed[self._id_col()].tolist()

    def filter_predictions(self, predictions: pd.DataFrame) -> pd.DataFrame:
        expressed = set(self.get_expressed_genes())
        gene_col = "TargetGene" if "TargetGene" in predictions.columns else "gene_name"
        filtered = predictions[predictions[gene_col].isin(expressed)]
        logger.info("Filtered %d -> %d predictions (%d expressed genes)",
                    len(predictions), len(filtered), len(expressed))
        return filtered

    def add_expression_weight(self, predictions: pd.DataFrame,
                              weight_method: str = "log") -> pd.DataFrame:
        if self._expression_data is None:
            self.load_expression()
        expr_dict = dict(zip(self._expression_data[self._id_col()],
                             self._expression_data[self.expression_column]))
        predictions = predictions.copy()
        gene_col = "TargetGene" if "TargetGene" in predictions.columns else "gene_name"
        predictions["expression"] = predictions[gene_col].map(expr_dict).fillna(0)

        if weight_method == "log":
            log_expr = np.log2(predictions["expression"] + 1)
            denom = log_expr.max()
            predictions["expression_weight"] = log_expr / denom if denom > 0 else 0.0
        elif weight_method == "linear":
            max_expr = predictions["expression"].max()
            predictions["expression_weight"] = (
                predictions["expression"] / max_expr if max_expr > 0 else 0.0)
        elif weight_method == "binary":
            predictions["expression_weight"] = (
                predictions["expression"] >= self.min_expression).astype(float)
        else:
            raise ValueError(f"Unknown weight method: {weight_method}")
        return predictions


class eQTLValidator:
    """Validate enhancer-gene predictions using eQTL data.

    Computes a real overlap between predicted enhancer regions and significant
    cis-eQTL variants and quantifies enrichment with Fisher's exact test,
    returning the full 2x2 contingency table, odds ratio and p-value.
    For genome-wide enrichment with size-matched controls and score
    stratification, see ``scripts/eqtl_enrichment.py``.
    """

    def __init__(self, eqtl_file: str,
                 chr_col: str = "chr", pos_col: str = "pos",
                 gene_col: str = "gene", pval_col: str = "pvalue",
                 beta_col: str = "beta"):
        self.eqtl_file = eqtl_file
        self.chr_col = chr_col
        self.pos_col = pos_col
        self.gene_col = gene_col
        self.pval_col = pval_col
        self.beta_col = beta_col
        self._eqtl_data: Optional[pd.DataFrame] = None

    def load_eqtl(self) -> pd.DataFrame:
        logger.info("Loading eQTL data from %s", self.eqtl_file)
        self._eqtl_data = pd.read_csv(self.eqtl_file, sep="\t")
        return self._eqtl_data

    def validate_predictions(self, predictions: pd.DataFrame,
                             pvalue_threshold: float = 1e-5,
                             window: int = 1000) -> pd.DataFrame:
        """Annotate each prediction with overlapping significant eQTL evidence."""
        if self._eqtl_data is None:
            self.load_eqtl()
        eqtl = self._eqtl_data
        sig = eqtl[eqtl[self.pval_col] <= pvalue_threshold].copy()

        predictions = predictions.copy()
        predictions["eQTL_support"] = False
        predictions["eQTL_pvalue"] = np.nan
        predictions["eQTL_beta"] = np.nan

        if len(sig) == 0:
            logger.warning("No significant eQTLs at p <= %g", pvalue_threshold)
            return predictions

        gene_col = "TargetGene" if "TargetGene" in predictions.columns else "gene_name"
        # Index eQTLs by gene for matching to predicted target gene.
        by_gene = {g: d for g, d in sig.groupby(self.gene_col)}

        for idx, row in predictions.iterrows():
            g = row[gene_col]
            if g not in by_gene:
                continue
            cand = by_gene[g]
            same_chr = cand[cand[self.chr_col].astype(str) == str(row["chr"])]
            in_window = same_chr[
                (same_chr[self.pos_col] >= row["start"] - window) &
                (same_chr[self.pos_col] <= row["end"] + window)
            ]
            if len(in_window) > 0:
                best = in_window.loc[in_window[self.pval_col].idxmin()]
                predictions.at[idx, "eQTL_support"] = True
                predictions.at[idx, "eQTL_pvalue"] = best[self.pval_col]
                predictions.at[idx, "eQTL_beta"] = best.get(self.beta_col, np.nan)

        logger.info("eQTL validation: %d/%d predictions supported",
                    int(predictions["eQTL_support"].sum()), len(predictions))
        return predictions

    def calculate_enrichment(self, predictions: pd.DataFrame,
                             score_column: str = "PACE.Score",
                             threshold: float = 0.02) -> Dict:
        """Fisher's exact enrichment of eQTL support in high-score predictions.

        Returns the contingency table, odds ratio and p-value.
        """
        from scipy.stats import fisher_exact
        if "eQTL_support" not in predictions.columns:
            predictions = self.validate_predictions(predictions)

        high = predictions[score_column] >= threshold
        supp = predictions["eQTL_support"].astype(bool)

        a = int((high & supp).sum())      # high-score, eQTL-supported
        b = int((high & ~supp).sum())     # high-score, not supported
        c = int((~high & supp).sum())     # low-score, supported
        d = int((~high & ~supp).sum())    # low-score, not supported

        table = [[a, b], [c, d]]
        odds_ratio, pvalue = fisher_exact(table, alternative="greater")
        return {
            "n_predictions": len(predictions),
            "contingency_table": table,
            "high_supported": a, "high_unsupported": b,
            "low_supported": c, "low_unsupported": d,
            "odds_ratio": odds_ratio,
            "pvalue": pvalue,
            "fraction_supported": float(supp.mean()),
        }


class MethylationIntegrator:
    """Integrate DNA methylation data as a regulatory inhibitor.

    High DNA methylation at an enhancer typically indicates a repressed /
    inactive state. Methylation levels are read directly from the data file
    (bigWig or bedGraph of per-region or per-CpG fractional methylation).
    """

    def __init__(self, methylation_file: str, context: str = "CG"):
        self.methylation_file = methylation_file
        self.context = context

    def get_methylation_scores(self, regions: pd.DataFrame) -> np.ndarray:
        """Return mean methylation fraction (0-1) per region from real data."""
        if not os.path.exists(self.methylation_file):
            raise FileNotFoundError(
                f"Methylation file not found: {self.methylation_file}")

        if len(regions) == 0:
            return np.array([], dtype=float)

        if _is_bigwig(self.methylation_file):
            vals = quantify_signal_over_regions(
                regions, self.methylation_file, stat="mean",
                region_size_norm=False)
        else:
            # bedGraph: chr start end value(0-1 or 0-100).
            # Sort both inputs to plain temp files (no bash process
            # substitution, which is unavailable under /bin/sh=dash).
            work = regions.copy().reset_index(drop=True)
            if "name" not in work.columns:
                work["name"] = [f"region_{i}" for i in range(len(work))]
            regions_bed = sorted_meth = out = None
            try:
                with tempfile.NamedTemporaryFile("w", suffix=".bed",
                                                 delete=False) as t:
                    work[["chr", "start", "end", "name"]].sort_values(
                        ["chr", "start"]).to_csv(
                        t, sep="\t", header=False, index=False)
                    regions_bed = t.name
                sorted_meth = regions_bed + ".meth.sorted"
                out = regions_bed + ".map"
                _run(f"sort -k1,1 -k2,2n {self.methylation_file} > {sorted_meth}")
                _run(f"bedtools map -a {regions_bed} -b {sorted_meth} "
                     f"-c 4 -o mean > {out}")
                res = pd.read_csv(out, sep="\t", header=None,
                                  names=["chr", "start", "end", "name", "m"])
                res["m"] = pd.to_numeric(res["m"], errors="coerce").fillna(0)
                m_map = dict(zip(res["name"], res["m"]))
                vals = work["name"].map(m_map).fillna(0).to_numpy(float)
            finally:
                for p in (regions_bed, sorted_meth, out):
                    if p and os.path.exists(p):
                        os.unlink(p)

        # Normalize percentage to fraction if necessary.
        if np.nanmax(vals) > 1.5:
            vals = vals / 100.0
        return np.clip(vals, 0, 1)

    def calculate_inhibitory_factor(self, methylation: np.ndarray,
                                    method: str = "linear") -> np.ndarray:
        if method == "linear":
            return methylation
        if method == "sigmoid":
            return 1 / (1 + np.exp(-10 * (methylation - 0.5)))
        raise ValueError(f"Unknown method: {method}")


def calculate_pace_activity(
    regions: pd.DataFrame,
    atac_file: str,
    h3k27ac_file: Optional[str] = None,
    h3k4me1_file: Optional[str] = None,
    h3k4me3_file: Optional[str] = None,
    tf_files: Optional[Dict[str, str]] = None,
    tf_inhibitory: Optional[Dict[str, bool]] = None,
    methylation_file: Optional[str] = None,
    method: AggregationMethod = AggregationMethod.WEIGHTED_GEOMETRIC,
    weights: Optional[Dict[str, float]] = None,
) -> np.ndarray:
    """Convenience function to calculate the PACE activity score.

    Args:
        regions: DataFrame with enhancer regions (chr/start/end).
        atac_file: ATAC-seq / DNase-seq signal file (required).
        h3k27ac_file/h3k4me1_file/h3k4me3_file: optional histone signal files.
        tf_files: mapping of TF name -> signal file.
        tf_inhibitory: mapping of TF name -> whether the TF is a repressor.
            Defaults to activator.
        methylation_file: optional methylation file (inhibitory).
        method: aggregation method.
        weights: custom per-signal weights.

    Returns:
        Activity scores per region.
    """
    weights = weights or {}
    tf_inhibitory = tf_inhibitory or {}

    calc = MultiOmicsActivityCalculator(method=method)
    calc.add_signal("ATAC", SignalType.ACCESSIBILITY, atac_file,
                    weight=weights.get("ATAC", 1.5))
    if h3k27ac_file:
        calc.add_signal("H3K27ac", SignalType.ACTIVE_ENHANCER, h3k27ac_file,
                        weight=weights.get("H3K27ac", 1.0))
    if h3k4me1_file:
        calc.add_signal("H3K4me1", SignalType.ENHANCER_MARK, h3k4me1_file,
                        weight=weights.get("H3K4me1", 0.8))
    if h3k4me3_file:
        calc.add_signal("H3K4me3", SignalType.PROMOTER_MARK, h3k4me3_file,
                        weight=weights.get("H3K4me3", 0.5))
    if tf_files:
        for tf_name, tf_file in tf_files.items():
            calc.add_signal(tf_name, SignalType.TF_BINDING, tf_file,
                            weight=weights.get(tf_name, 0.3),
                            inhibitory=tf_inhibitory.get(tf_name, False))
    if methylation_file:
        calc.add_signal("methylation", SignalType.METHYLATION, methylation_file,
                        weight=weights.get("methylation", 0.5), inhibitory=True)
    return calc.calculate_activity(regions)


if __name__ == "__main__":
    print("PACE Multi-omics Activity Calculator")
    print("=" * 50)
    print("This module quantifies real epigenomic signals over regions.")
    print("Provide bigWig/BED/tagAlign/BAM files to compute activity scores.")
    print("See scripts/calculate_pace_score.py for the full pipeline.")
