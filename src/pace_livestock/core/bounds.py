"""Sensitivity bounds for a fixed candidate universe, not confidence intervals."""

import math

import numpy as np

from ..errors import PaceError


def score_intervals(rows):
    """Bound S_i / sum(S) given exact or explicitly bounded nonnegative support.

    Unknown support has [0, infinity) unless the caller supplies justified bounds.
    Bounds condition on positive total support. Prefix/suffix sums avoid both
    quadratic work and subtraction of nearly equal, very large numbers.
    """
    lower, upper, kinds = [], [], []
    for row in rows:
        if row["scoreable"]:
            lower.append(row["log_support"])
            upper.append(row["log_support"])
            kinds.append("exact_resolved_support")
            continue
        lo = row.get("support_lower", 0.0)
        hi = row.get("support_upper", math.inf)
        hi = math.inf if hi is None or math.isnan(hi) else hi
        if not math.isfinite(lo) or lo < 0 or hi < lo:
            raise PaceError("Missing support bounds require 0 <= lower <= upper")
        lower.append(math.log(lo) if lo else -math.inf)
        upper.append(math.log(hi) if hi else -math.inf)
        kinds.append("provided_missing_support" if row.get("bound_source") else "nonnegative_only")

    def excluding(values):
        n = len(values)
        prefix, suffix = [-math.inf] * (n + 1), [-math.inf] * (n + 1)
        for i in range(n):
            prefix[i + 1] = float(np.logaddexp(prefix[i], values[i]))
            suffix[n - i - 1] = float(np.logaddexp(suffix[n - i], values[n - i - 1]))
        return [float(np.logaddexp(prefix[i], suffix[i + 1])) for i in range(n)], prefix[-1]

    other_lo, _ = excluding(lower)
    other_hi, possible_total = excluding(upper)
    for i, row in enumerate(rows):
        lo, hi = lower[i], upper[i]
        if possible_total == -math.inf:
            a = b = math.nan
        else:
            a = (
                1.0
                if other_hi[i] == -math.inf and hi > -math.inf
                else 0.0
                if lo == -math.inf or other_hi[i] == math.inf
                else math.exp(lo - float(np.logaddexp(lo, other_hi[i])))
            )
            b = (
                0.0
                if hi == -math.inf
                else 1.0
                if hi == math.inf or other_lo[i] == -math.inf
                else math.exp(hi - float(np.logaddexp(hi, other_lo[i])))
            )
        row.update(pace_score_lo=a, pace_score_hi=b, bound_basis=kinds[i])
