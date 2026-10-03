"""Vectorised implementation of the PACE score (Equations 1-4) for whole-genome benchmarks.

The reference implementation (``pace run``) resolves every element-TSS pair in pure Python
and needs tens of GB for the ~37 million pairs of the K562 benchmark.  This module evaluates
the same equations with numpy, one chromosome at a time.  ``manuscript/verify_engine.py``
checks that it reproduces the scores of ``pace run`` on a gene subset.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import pandas as pd

from pace_livestock.evidence.contact import fit_binned_prior

HUMAN_ABC_GAMMA = 1.024238616787792


# --------------------------------------------------------------------------- priors
@dataclass
class Prior:
    """Equation 2: P(d) = a (max(d, d_min) / d_ref)^-gamma, evaluated at bin centres."""

    gamma: float
    a: float = 1.0
    d_ref: float = 5000.0
    d_min: float = 10000.0
    res: int = 10000
    centered: bool = True  # False: evaluate at the genomic anchors themselves (the built-in human preset)

    def __call__(self, d):
        d = np.maximum(np.asarray(d, dtype=float), self.d_min)
        return self.a * (d / self.d_ref) ** (-self.gamma)


def human_reference_prior() -> Prior:
    """Case C: the fixed human ABC reference decay (relative scale, 5 kb bin centres)."""
    return Prior(gamma=HUMAN_ABC_GAMMA, a=1.0, d_ref=5000.0, d_min=5000.0, res=5000, centered=False)


# --------------------------------------------------------------------------- Hi-C bands
def ice_balance(band: np.ndarray, ignore_diags: int = 2, mad_max: float = 5.0, min_nnz: int = 10,
                tol: float = 1e-5, max_iters: int = 200, maxd_bins: int | None = None) -> np.ndarray:
    """Iterative correction of a banded cis matrix (cooler defaults).  Returns weights (NaN = masked)."""
    from scipy import sparse
    n, K = band.shape
    if maxd_bins is not None:
        K = min(K, maxd_bins)
    b = band[:, :K].astype(np.float64)
    b[:, :ignore_diags] = 0.0
    U = np.zeros((K, n))
    for k in range(K):
        U[k, k:] = b[: n - k, k]
    U = sparse.dia_array((U, np.arange(K)), shape=(n, n)).tocsr()
    S = (U + U.T).tocsr()  # diagonal ignored (ignore_diags >= 1), so no double counting
    marg = np.asarray(S.sum(axis=1)).ravel()
    nnz = np.diff(S.indptr)
    with np.errstate(divide="ignore"):
        lm = np.log(marg[marg > 0])
    low = -np.inf
    if len(lm):
        med = np.median(lm)
        low = med - mad_max * 1.4826 * np.median(np.abs(lm - med))
    bad = (marg <= 0) | (np.log(np.where(marg > 0, marg, 1)) < low) | (nnz < min_nnz)
    w = np.where(bad, 0.0, 1.0)
    valid = ~bad
    for _ in range(max_iters):
        m = w * (S @ w)
        mean = m[valid].mean() if valid.any() else 1.0
        bias = np.where(valid & (m > 0), m / mean, 1.0)
        w = np.where(valid, w / bias, 0.0)
        if np.max(np.abs(bias[valid] - 1.0)) < tol:
            break
    m = w * (S @ w)
    w = w / math.sqrt(m[valid].mean() if valid.any() else 1.0)
    w[bad] = np.nan
    return w


def distance_bins(bands: dict, weights: dict, res: int, min_distance=None, max_distance=5_000_000,
                  n_bins=30):
    """Per-distance-class sums of balanced contact and valid-bin opportunities (as ``pace fit-prior``)."""
    min_distance = res if min_distance is None else min_distance
    edges = np.geomspace(min_distance, max_distance, n_bins + 1)
    edges[0], edges[-1] = min_distance, max_distance
    n_pairs = np.zeros(n_bins)
    logd = np.zeros(n_bins)
    csum = np.zeros(n_bins)
    for ch, band in bands.items():
        w = weights[ch]
        valid = np.isfinite(w) & (w > 0)
        n = len(valid)
        K = band.shape[1]
        wz = np.where(valid, w, 0.0)
        fft = 1 << max(1, (2 * n - 1).bit_length())
        t = np.fft.rfft(valid.astype(float), n=fft)
        opp = np.maximum(np.rint(np.fft.irfft(t * t.conj(), n=fft)[:n]), 0)
        ks = np.arange(n)
        dist = ks * res
        idx = np.searchsorted(edges, dist, side="right") - 1
        use = (idx >= 0) & (idx < n_bins)
        n_pairs += np.bincount(idx[use], weights=opp[use], minlength=n_bins)
        logd += np.bincount(idx[use], weights=opp[use] * np.log(np.maximum(dist[use], 1)), minlength=n_bins)
        for k in range(K):
            d = k * res
            j = np.searchsorted(edges, d, side="right") - 1
            if j < 0 or j >= n_bins or k >= n:
                continue
            csum[j] += float(np.dot(band[: n - k, k].astype(np.float64), wz[: n - k] * wz[k:]))
    rows = []
    for i in range(n_bins):
        n = n_pairs[i]
        rows.append(dict(low=edges[i], high=edges[i + 1], n_pairs=int(n),
                         contact_sum=csum[i], log_distance_sum=logd[i],
                         mean_contact=csum[i] / n if n else None,
                         geometric_distance=math.exp(logd[i] / n) if n else None,
                         fit_status="used" if n and csum[i] > 0 else ("zero_mean_excluded" if n else "empty")))
    return rows


def fit_prior(bands: dict, weights: dict, res: int, d_ref=5000.0) -> tuple[Prior, list]:
    rows = distance_bins(bands, weights, res)
    fit = fit_binned_prior(rows, d_ref=d_ref, d_min=float(res))
    return Prior(gamma=fit["gamma"], a=fit["a"], d_ref=d_ref, d_min=float(res), res=res), rows


# --------------------------------------------------------------------------- activity
def element_activity(obs: pd.DataFrame, panel: dict, units_index: pd.Index, reps: dict | None = None):
    """Equation 1 activity: replicate means per assay, then the geometric mean over the panel.

    ``obs`` is the long observed-activity table of ``pace activity``; ``panel`` maps assay -> list of
    sample ids to use.  Unmeasured replicates are skipped; an assay without any measurement gives NA.
    A measured zero in any assay gives zero activity.
    """
    cols = []
    for assay, samples in panel.items():
        d = obs[(obs.sample_id.isin(samples)) & (obs.measurement_status == "observed")]
        m = d.groupby("element_id")["signal"].mean().reindex(units_index)
        cols.append(m.to_numpy(float))
    X = np.vstack(cols)
    with np.errstate(divide="ignore", invalid="ignore"):
        g = np.exp(np.mean(np.log(np.where(X > 0, X, 1.0)), axis=0))
    g[np.any(X == 0, axis=0)] = 0.0
    g[np.any(np.isnan(X), axis=0)] = np.nan
    return g


# --------------------------------------------------------------------------- scoring
def _expand(cand_gene, ntss, tss_start):
    """Index arrays that repeat every candidate row once per TSS of its gene."""
    reps = ntss[cand_gene]
    row = np.repeat(np.arange(len(cand_gene)), reps)
    first = np.cumsum(reps) - reps
    within = np.arange(len(row)) - np.repeat(first, reps)
    tss_idx = tss_start[cand_gene][row] + within
    return row, tss_idx, first


def contact_values(anchor, tss, mode, *, prior, band=None, weight=None, res=None, kappa=1.0,
                   d0=5000.0, fallback=True):
    """Resolve C~(E, t) for arrays of element anchors and TSS positions (Equation 3 / Cases A-C)."""
    pres = prior.res
    if prior.centered:
        p_d = prior(np.abs(anchor // pres - tss // pres) * pres)
    else:
        p_d = prior(np.abs(anchor - tss))
    if mode in ("B", "C"):
        return p_d
    assert mode == "A"
    i = anchor // res
    j = tss // res
    lo = np.minimum(i, j)
    k = np.abs(i - j)
    same = k == 0
    inside = k < band.shape[1]
    kk = np.where(inside, k, 0)
    lo = np.where(inside, lo, 0)
    obs = band[lo, kk].astype(np.float64) * weight[lo] * weight[lo + kk]
    obs[~inside] = np.nan
    valid = np.isfinite(obs)
    pc = kappa * np.minimum(p_d, prior(d0))
    out = np.where(valid, obs + pc, p_d if fallback else np.nan)
    out = np.where(same, prior(0.0) * np.ones_like(p_d), out)  # near-diagonal pairs use the prior
    return out


def score_chromosome(units: pd.DataFrame, A: np.ndarray, cand_u: np.ndarray, cand_g: np.ndarray,
                     genes: pd.DataFrame, tss: pd.DataFrame, mode: str, prior: Prior, band=None,
                     weight=None, res=None, kappa=1.0, d0=5000.0, fallback=True):
    """Score the candidate pairs (cand_u, cand_g) of one chromosome.

    units: anchor0 column; A: activity per unit; genes: one row per gene (gene index = row);
    tss: columns gene_idx, tss0, pi (sorted by gene_idx).
    Returns arrays (cbar, score_withhold, score_conditional) aligned to the candidates.
    """
    ntss = np.bincount(tss.gene_idx.to_numpy(), minlength=len(genes))
    tss_start = np.cumsum(ntss) - ntss
    row, ti, first = _expand(cand_g, ntss, tss_start)
    anchor = units.anchor0.to_numpy()[cand_u][row]
    c = contact_values(anchor, tss.tss0.to_numpy()[ti], mode, prior=prior, band=band, weight=weight,
                       res=res, kappa=kappa, d0=d0, fallback=fallback)
    pi = tss.pi.to_numpy()[ti]
    cbar = np.add.reduceat(c * pi, first) if len(first) else np.zeros(0)
    nan_any = np.add.reduceat(np.isnan(c).astype(np.int8), first) > 0 if len(first) else np.zeros(0, bool)
    cbar[nan_any] = np.nan
    a = A[cand_u]
    support = a * cbar
    support = np.where((a == 0) & np.isfinite(cbar), 0.0, support)
    sc_with, sc_cond = _normalise(support, cand_g, len(genes))
    return cbar, sc_with, sc_cond


def _normalise(support, gene, ngenes):
    ok = np.isfinite(support)
    tot_ok = np.bincount(gene[ok], weights=support[ok], minlength=ngenes)
    n_bad = np.bincount(gene[~ok], minlength=ngenes)
    with np.errstate(divide="ignore", invalid="ignore"):
        cond = np.where(ok, support / tot_ok[gene], np.nan)
    cond[~np.isfinite(cond)] = np.nan
    with_ = np.where(n_bad[gene] == 0, cond, np.nan)
    return with_, cond
