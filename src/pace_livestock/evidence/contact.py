"""Distance priors and explicit linear shrinkage; no implicit human-species defaults."""

import math

import numpy as np

from ..core.scoring import nonnegative
from ..errors import PaceError


def distance_prior(distance: float, prior: dict) -> float:
    for k in ("a", "gamma", "d_min", "d_ref"):
        if not math.isfinite(prior[k]) or prior[k] <= 0:
            raise PaceError(f"Contact prior {k} must be positive and finite")
    return math.exp(
        math.log(prior["a"])
        - prior["gamma"] * (math.log(max(distance, prior["d_min"])) - math.log(prior["d_ref"]))
    )


def shrink(observed: float, prior: float, reliability: float) -> float:
    if not math.isfinite(reliability) or not 0 <= reliability <= 1:
        raise PaceError("Contact reliability must be in [0,1]")
    if reliability == 0:
        return float(nonnegative([prior])[0])
    if reliability == 1:
        return float(nonnegative([observed])[0])
    c = nonnegative([observed, prior])
    return float(c @ np.array([reliability, 1 - reliability]))


def fit_distance_prior(distances, contacts, *, bin_edges, d_ref: float, d_min: float) -> dict:
    d, c = nonnegative(distances), nonnegative(contacts)
    boundaries = np.asarray(bin_edges, dtype=float)
    if d.shape != c.shape or np.any(~np.isfinite(d)) or np.any(~np.isfinite(c)):
        raise PaceError(
            "Prior fitting requires finite paired distances and contacts, including zeros"
        )
    if (
        d_ref <= 0
        or d_min <= 0
        or np.any(~np.isfinite(boundaries))
        or len(boundaries) < 3
        or np.any(np.diff(boundaries) <= 0)
    ):
        raise PaceError("Prior fit needs increasing bin edges and positive d_ref/d_min")
    bins = []
    for low, high in zip(boundaries[:-1], boundaries[1:], strict=True):
        mask = (d >= low) & (d < high)
        n = int(mask.sum())
        mean = float(np.mean(c[mask])) if n else None
        center = float(np.exp(np.mean(np.log(np.maximum(d[mask], d_min))))) if n else None
        bins.append(
            {
                "low": float(low),
                "high": float(high),
                "n_pairs": n,
                "mean_contact": mean,
                "geometric_distance": center,
                "fit_status": "used" if n and mean > 0 else "zero_mean_excluded" if n else "empty",
            }
        )
    used = [b for b in bins if b["fit_status"] == "used"]
    if len(used) < 2:
        raise PaceError("Prior fit needs at least two nonzero mean distance bins")
    x = np.log([b["geometric_distance"] / d_ref for b in used])
    y = np.log([b["mean_contact"] for b in used])
    w = np.sqrt([b["n_pairs"] for b in used])
    design = np.column_stack([np.ones(len(x)), x]) * w[:, None]
    if np.linalg.matrix_rank(design) < 2:
        raise PaceError("Contact prior slope is unidentifiable")
    beta = np.linalg.lstsq(design, y * w, rcond=None)[0]
    a, gamma = math.exp(beta[0]), -float(beta[1])
    if gamma <= 0 or not math.isfinite(a):
        raise PaceError("Fitted data do not follow a positive decreasing power-law prior")
    return {
        "a": a,
        "gamma": gamma,
        "d_ref": d_ref,
        "d_min": d_min,
        "bins": bins,
        "fit_weight": "valid_pair_count",
        "fit_range": [float(boundaries[0]), float(boundaries[-1])],
        "n_pairs_outside_bins": int(np.sum((d < boundaries[0]) | (d >= boundaries[-1]))),
    }
