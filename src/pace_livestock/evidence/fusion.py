"""Closed-form constrained calibration against an independently measured target."""

import math

import numpy as np

from ..core.scoring import fuse, nonnegative
from ..errors import PaceError


def fit_fusion(observed, predicted, target, *, scale: float, minimum_samples: int) -> dict:
    arrays = [nonnegative(x) for x in (observed, predicted, target)]
    if scale <= 0 or not math.isfinite(scale) or minimum_samples < 2:
        raise PaceError(
            "Fusion requires positive scale and a protocol-defined minimum_samples >= 2"
        )
    if len({a.shape for a in arrays}) != 1 or arrays[0].ndim != 1:
        raise PaceError("Fusion vectors must have the same one-dimensional shape")
    valid = np.all(np.isfinite(arrays), axis=0)
    if int(valid.sum()) < minimum_samples:
        raise PaceError("Too few complete independent calibration measurements")
    obs, seq, truth = [
        np.logaddexp(
            np.log(a[valid], where=a[valid] > 0, out=np.full(valid.sum(), -np.inf)), math.log(scale)
        )
        - math.log(scale)
        for a in arrays
    ]
    delta = obs - seq
    denominator = float(delta @ delta)
    if denominator == 0:
        return {
            "identifiable": False,
            "weight": None,
            "scale": scale,
            "n": int(valid.sum()),
            "reason": "identical_sources",
            "fallback": "observed_then_prediction",
        }
    weight = float(np.clip(delta @ (truth - seq) / denominator, 0, 1))
    return {
        "identifiable": True,
        "weight": weight,
        "scale": scale,
        "n": int(valid.sum()),
        "calibration_mse": {
            "observed": float(np.mean((obs - truth) ** 2)),
            "sequence": float(np.mean((seq - truth) ** 2)),
            "fused": float(np.mean((weight * obs + (1 - weight) * seq - truth) ** 2)),
        },
    }


def resolve_signal(
    observed: float, predicted: float, *, regime: str, calibration=None
) -> tuple[float, str, float]:
    if regime == "measured":
        return observed, "observed", 1.0
    if regime == "genome_only":
        return predicted, "sequence_prediction", 0.0
    if calibration and calibration["identifiable"]:
        w = calibration["weight"]
        if w == 0 or w == 1 or (math.isfinite(observed) and math.isfinite(predicted)):
            return (
                fuse(observed, predicted, weight=w, scale=calibration["scale"]),
                "fused" if 0 < w < 1 else "observed" if w == 1 else "sequence_prediction",
                w,
            )
    if math.isfinite(observed):
        return observed, "observed", 1.0
    return predicted, "sequence_prediction", 0.0
