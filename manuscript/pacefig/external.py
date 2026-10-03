"""Comparison methods: ENCODE ABC predictions, distance baselines and a supervised GBM."""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier

from . import bench


def map_intervals(pairs: pd.DataFrame, pred: pd.DataFrame, cols=("ABC.Score",)):
    """Max prediction over predicted elements overlapping each tested element of the same gene."""
    p = pairs[["chrom", "chromStart", "chromEnd", "gene_id"]].reset_index().rename(columns={"index": "pid"})
    m = p.merge(pred, on=["chrom", "gene_id"], how="inner")
    m = m[(m.start < m.chromEnd) & (m.end > m.chromStart)]
    out = {}
    for c in cols:
        s = m.groupby("pid")[c].max()
        a = np.full(len(pairs), np.nan)
        a[s.index.to_numpy()] = s.to_numpy()
        out[c] = a
    return out


def load_abc():
    return pd.read_csv(bench.WORK + "abc_encode_k562.tsv.gz", sep="\t")


FEATURES = ["DHS.RPM", "H3K27ac.RPM", "H3K27ac.RPM.expandedRegion", "H3K4me1.RPM", "CTCF.RPM", "H3K27me3.RPM",
            "H3K27ac_peak_overlap", "CTCF_peak_overlap", "H3K27me3_peak_overlap", "H3K4me1_peak_overlap"]


def gbm_features(pairs: pd.DataFrame, extra: dict):
    X = pairs[FEATURES].copy()
    X["log_distance"] = np.log10(pairs.distanceToTSS.clip(lower=1))
    for k, v in extra.items():
        X[k] = v
    return X


def gbm_scores(pairs: pd.DataFrame, X: pd.DataFrame, seed=0):
    """Leave-one-chromosome-out predictions for the combined benchmark; the held-out set is scored by a model
    trained on the whole combined benchmark."""
    y = pairs.y.to_numpy()
    main = (pairs.set == "main").to_numpy()
    held = ~main
    s = np.full(len(pairs), np.nan)
    chroms = pairs.chrom.to_numpy()

    def model():
        return HistGradientBoostingClassifier(max_depth=4, learning_rate=0.05, max_iter=200, min_samples_leaf=30,
                                              l2_regularization=1.0, random_state=seed)

    for ch in np.unique(chroms[main]):
        te = main & (chroms == ch)
        tr = main & (chroms != ch)
        s[te] = model().fit(X[tr], y[tr]).predict_proba(X[te])[:, 1]
    s[held] = model().fit(X[main], y[main]).predict_proba(X[held])[:, 1]
    return s
