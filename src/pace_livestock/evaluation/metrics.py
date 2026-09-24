"""Tie-aware average precision and AUROC, with coverage kept separate from scores."""

import math

import numpy as np

from ..errors import PaceError


def binary_metrics(labels, scores, *, threshold=None):
    # Validate before integer conversion: fractional labels must never become
    # apparently valid negatives/positives through silent truncation.
    try:
        y, s = np.asarray(labels, dtype=float), np.asarray(scores, dtype=float)
    except (TypeError, ValueError) as exc:
        raise PaceError("Metrics require binary labels and numeric scores") from exc
    if y.shape != s.shape or y.ndim != 1 or not np.all(np.isin(y, [0, 1])):
        raise PaceError("Metrics require paired binary labels and one-dimensional scores")
    y = y.astype(int)
    valid = np.isfinite(s)
    missed_positive = int(np.sum((y == 1) & ~valid))
    all_positive = int(y.sum())
    yv, sv = y[valid], s[valid]
    result = {
        "n_tested": len(y),
        "n_scored": int(valid.sum()),
        "n_positive": all_positive,
        "coverage": float(valid.mean()) if len(y) else math.nan,
        "missed_positive": missed_positive,
        "average_precision": math.nan,
        "auroc": math.nan,
        "precision": math.nan,
        "recall_candidate": math.nan,
        "recall_end_to_end": math.nan,
        "threshold": threshold,
        "reason": "resolved",
    }
    if len(yv) == 0 or yv.sum() == 0 or yv.sum() == len(yv):
        result["reason"] = "no_scores_or_single_class"
    else:
        order = np.argsort(-sv, kind="stable")
        ys, ss = yv[order], sv[order]
        ends = np.r_[np.flatnonzero(np.diff(ss)), len(ss) - 1]
        tp = np.cumsum(ys)[ends]
        fp = (ends + 1) - tp
        recall = tp / yv.sum()
        precision = tp / (ends + 1)
        result["average_precision"] = float(np.sum(np.diff(np.r_[0, recall]) * precision))
        # Trapezoidal ROC (unlike AP); ties enter together, giving half credit.
        fpr, tpr = np.r_[0, fp / (len(yv) - yv.sum())], np.r_[0, recall]
        result["auroc"] = float(np.sum(np.diff(fpr) * (tpr[:-1] + tpr[1:]) / 2))
    if threshold is not None:
        if not math.isfinite(threshold):
            raise PaceError("Frozen threshold must be finite")
        called = sv >= threshold
        tp = int(yv[called].sum())
        result["precision"] = tp / int(called.sum()) if called.any() else math.nan
        result["recall_candidate"] = tp / int(yv.sum()) if yv.sum() else math.nan
        result["recall_end_to_end"] = tp / all_positive if all_positive else math.nan
    return result
