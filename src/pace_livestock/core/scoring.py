"""PACE September 2026 equations, evaluated in log space without pseudocounts."""

from __future__ import annotations

import math
from collections import defaultdict

import numpy as np

from ..errors import PaceError
from ..provenance import digest


def nonnegative(values) -> np.ndarray:
    x = np.asarray(values, dtype=np.float64)
    if np.any(np.isinf(x)) or np.any(x < 0):
        raise PaceError("Signals must be finite nonnegative values or NA")
    return x


def activity(values) -> float:
    x = nonnegative(values)
    if x.size == 0 or np.any(np.isnan(x)):
        return math.nan
    if np.any(x == 0):
        return 0.0
    return float(np.exp(np.mean(np.log(x))))


def bulk_mean(values, weights=None) -> np.ndarray:
    """Average marginal assay signals across declared copies BEFORE constructing activity."""
    x = nonnegative(values)
    if x.ndim < 1 or len(x) == 0:
        raise PaceError("bulk_mean requires at least one copy")
    w = np.full(len(x), 1 / len(x)) if weights is None else nonnegative(weights)
    if (
        w.shape != (len(x),)
        or not np.all(np.isfinite(w))
        or not math.isclose(float(w.sum()), 1, abs_tol=1e-12)
    ):
        raise PaceError("Copy weights must be finite, nonnegative and sum to one")
    # Zero-weight sources are not required; scaled summation avoids finite-value overflow.
    used = w > 0
    return np.sum(x[used] * w[used].reshape((-1,) + (1,) * (x.ndim - 1)), axis=0)


def fuse(observed: float, predicted: float, *, weight: float, scale: float) -> float:
    if not math.isfinite(weight) or not 0 <= weight <= 1 or not math.isfinite(scale) or scale <= 0:
        raise PaceError("Fusion requires weight in [0,1] and positive finite scale")
    if weight == 1:
        return float(nonnegative([observed])[0])
    if weight == 0:
        return float(nonnegative([predicted])[0])
    x = nonnegative([observed, predicted])
    # logaddexp avoids overflow in x / scale for otherwise finite input.
    z = np.logaddexp(np.log(x, where=x > 0, out=np.full(2, -np.inf)), math.log(scale)) - math.log(
        scale
    )
    combined = float(z @ np.array([weight, 1 - weight]))
    if combined < 700:
        return scale * math.expm1(combined)
    return safe_exp(math.log(scale) + combined) - scale


def logsum(values) -> float:
    x = np.asarray(values, dtype=float)
    if not len(x):
        return -math.inf
    return float(np.logaddexp.reduce(x))


def safe_exp(value: float) -> float:
    if math.isnan(value) or value > math.log(np.finfo(float).max):
        return math.nan
    return math.exp(value)


def log_normalize(log_support) -> tuple[np.ndarray, float]:
    x = np.asarray(log_support, dtype=float)
    valid = ~np.isnan(x)
    total = logsum(x[valid])
    result = np.full(len(x), np.nan)
    if math.isfinite(total):
        result[valid] = np.exp(x[valid] - total)
    return result, total


def tss_contact(contacts, weights) -> float:
    c, w = nonnegative(contacts), nonnegative(weights)
    if (
        c.shape != w.shape
        or not np.all(np.isfinite(w))
        or not math.isclose(float(w.sum()), 1, abs_tol=1e-10)
    ):
        raise PaceError("TSS weights must be finite and sum to one for the complete gene")
    if np.any(np.isnan(c[w > 0])):
        return math.nan
    return float(np.sum(c[w > 0] * w[w > 0]))


def score(edges: list[dict], *, eta: int = 0) -> tuple[list[dict], list[dict]]:
    """Score sparse E–G records carrying A_used and Cbar; never drop an input edge.

    All planned genes for each element MUST be present, including missing contacts.
    This pure kernel assumes the catalog and promoter mappings were validated upstream.
    """
    if eta not in (0, 1):
        raise PaceError("eta must be 0 or 1")
    rows = [dict(r) for r in sorted(edges, key=lambda r: (r["gene_id"], r["element_id"]))]
    seen, by_element, by_gene = set(), defaultdict(list), defaultdict(list)
    for i, row in enumerate(rows):
        key = row["element_id"], row["gene_id"]
        if key in seen:
            raise PaceError(f"Duplicate candidate edge: {key}")
        seen.add(key)
        nonnegative([row["A_used"], row["Cbar"]])
        by_element[key[0]].append(i)
        by_gene[key[1]].append(i)
    for indices in by_element.values():
        c = np.array([rows[i]["Cbar"] for i in indices])
        complete = np.all(~np.isnan(c))
        log_csum = (
            logsum(np.log(c, where=c > 0, out=np.full(len(c), -np.inf))) if complete else math.nan
        )
        for i in indices:
            row = rows[i]
            a, contact = row["A_used"], row["Cbar"]
            row["B"] = math.nan
            row["log_support"] = math.nan
            reason = row.get("reason", "")
            if math.isnan(a):
                reason = reason or "unresolved_activity"
            elif math.isnan(contact):
                reason = reason or "unresolved_contact"
            elif eta and not complete:
                reason = reason or "fixed_gene_set_contact_missing"
            elif a == 0 or contact == 0:
                row["log_support"] = -math.inf
                if eta and math.isfinite(log_csum):
                    row["B"] = 0.0 if contact == 0 else math.exp(math.log(contact) - log_csum)
                reason = "zero_support"
            else:
                log_b = math.log(contact) - log_csum if eta else 0.0
                row["B"] = math.exp(log_b) if eta else math.nan
                row["log_support"] = math.log(a) + math.log(contact) + log_b
                reason = "resolved"
            row["support"] = safe_exp(row["log_support"])
            row["support_status"] = (
                "overflow"
                if math.isfinite(row["log_support"]) and math.isnan(row["support"])
                else reason
            )
            row["reason"] = reason
            # TSV cannot store -inf; structural/observed zero remains recoverable as support=0.
            row["scoreable"] = not math.isnan(row["log_support"])
    summaries = []
    for gene, indices in sorted(by_gene.items()):
        logs = [rows[i]["log_support"] for i in indices]
        scores, log_total = log_normalize(logs)
        valid = [i for i in indices if rows[i]["scoreable"]]
        state = (
            "empty"
            if not valid
            else "zero_support"
            if log_total == -math.inf
            else "complete"
            if len(valid) == len(indices)
            else "partial"
        )
        universe = digest(sorted(rows[i]["element_id"] for i in valid))
        summary = {
            "gene_id": gene,
            "n_candidates": len(indices),
            "n_scoreable": len(valid),
            "n_positive": sum(math.isfinite(rows[i]["log_support"]) for i in valid),
            "entry_coverage": len(valid) / len(indices),
            "denominator": safe_exp(log_total),
            "log_denominator": log_total,
            "normalization_status": state,
            "normalization_universe_id": universe,
        }
        summaries.append(summary)
        for i, p in zip(indices, scores, strict=True):
            rows[i].update(
                {
                    k: summary[k]
                    for k in [
                        "denominator",
                        "log_denominator",
                        "normalization_status",
                        "normalization_universe_id",
                    ]
                }
            )
            rows[i]["pace_score"] = float(p)
    return rows, summaries
