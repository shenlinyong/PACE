"""Quantify nonnegative window means with explicit missing-pixel semantics."""

import math

import numpy as np

from ..errors import PaceError


def quantify_bigwig(path, units, *, missing_is_measured_zero=False, minimum_callable_fraction=0.0):
    if type(missing_is_measured_zero) is not bool:
        raise PaceError("missing_is_measured_zero must be a YAML boolean")
    try:
        import pyBigWig
    except ImportError as exc:
        raise PaceError("bigWig support requires pip install 'pace-livestock[io]'") from exc
    if not 0 <= minimum_callable_fraction <= 1:
        raise PaceError("minimum_callable_fraction must be in [0,1]")
    rows = []
    with pyBigWig.open(str(path)) as bw:
        sizes = bw.chroms()
        for unit in units:
            chrom, start, end = unit["chrom"], int(unit["start"]), int(unit["end"])
            if chrom not in sizes or end > sizes[chrom] or start < 0 or end <= start:
                raise PaceError(f"bigWig/reference mismatch for {chrom}:{start}-{end}")
            values = np.asarray(bw.values(chrom, start, end, numpy=True), dtype=float)
            valid = np.isfinite(values)
            if np.any(values[valid] < 0) or np.any(np.isinf(values)):
                raise PaceError("Negative or infinite bigWig signal cannot be used as activity")
            stored_fraction = float(valid.mean())
            if missing_is_measured_zero:
                values = np.where(valid, values, 0)
                fraction = 1.0
            else:
                values = values[valid]
                fraction = stored_fraction
            status = (
                "observed"
                if len(values) and fraction >= minimum_callable_fraction
                else "low_coverage"
                if len(values)
                else "unmeasured"
            )
            rows.append(
                {
                    "element_id": unit["element_id"],
                    "signal": float(values.mean()) if status == "observed" else math.nan,
                    "callable_fraction": fraction,
                    "stored_fraction": stored_fraction,
                    "measurement_status": status,
                }
            )
    return rows
