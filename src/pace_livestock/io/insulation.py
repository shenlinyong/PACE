"""Domain boundaries from a Hi-C map by the insulation score (Crane et al. 2015)."""

import numpy as np

from ..errors import PaceError


def insulation_scores(band, valid, window):
    """Mean contact in the window x window square that spans each bin, from its band.

    band[k, a] holds the contact between bins a and a+k (k = 1..2*window); missing or
    invalid pixels are NaN. Returns log2(insulation / chromosome mean), NaN where undefined.
    """
    n = band.shape[1]
    total, count = np.zeros(n), np.zeros(n)
    i = np.arange(n)
    for k in range(2, 2 * window + 1):
        values = band[k, : n - k]
        ok = np.isfinite(values) & valid[: n - k] & valid[k:]
        prefix_sum = np.r_[0.0, np.cumsum(np.where(ok, values, 0.0))]
        prefix_n = np.r_[0, np.cumsum(ok)]
        lo = np.clip(np.maximum(i - window, i + 1 - k), 0, n - k)
        hi = np.clip(np.minimum(i - 1, i + window - k), -1, n - k - 1)
        use = hi >= lo
        total[use] += prefix_sum[hi[use] + 1] - prefix_sum[lo[use]]
        count[use] += prefix_n[hi[use] + 1] - prefix_n[lo[use]]
    score = np.full(n, np.nan)
    enough = valid & (count >= 0.5 * window * window) & (total > 0)
    score[enough] = total[enough] / count[enough]
    if np.any(enough):
        score[enough] = np.log2(score[enough] / np.mean(score[enough]))
    return score


def boundary_calls(score, window, *, min_depth=0.1):
    """Local minima of the insulation score whose depth below both flanks is >= min_depth."""
    calls = []
    n = len(score)
    for i in np.flatnonzero(np.isfinite(score)):
        left = score[max(0, i - window) : i]
        right = score[i + 1 : i + 1 + window]
        if not len(left) or not len(right):
            continue
        near = score[max(0, i - window // 2) : i + window // 2 + 1]
        if score[i] > np.nanmin(near):
            continue
        if np.all(np.isnan(left)) or np.all(np.isnan(right)):
            continue
        depth = min(np.nanmax(left), np.nanmax(right)) - score[i]
        if depth < min_depth or i in (0, n - 1):
            continue
        # Equal minima in neighbouring bins are one boundary; keep the first.
        if calls and i - calls[-1][0] <= window // 2:
            if depth > calls[-1][1]:
                calls[-1] = (int(i), float(depth))
            continue
        calls.append((int(i), float(depth)))
    return calls


def insulation_boundaries(uri, *, balanced=True, window_bp=500_000, min_depth=0.1):
    """Boundary points at bin centres; strength is the depth scaled by its 95th percentile."""
    try:
        import cooler
    except ImportError as exc:
        raise PaceError("Hi-C boundaries require the io extra") from exc
    c = cooler.Cooler(str(uri))
    resolution = int(c.binsize)
    window = max(2, round(window_bp / resolution))
    matrix = c.matrix(balance=balanced, sparse=True)
    calls = []
    for chrom in c.chromnames:
        bins = c.bins().fetch(chrom)
        n, size = len(bins), int(c.chromsizes[chrom])
        if n < 2 * window + 2:
            continue
        weights = bins["weight"].to_numpy(float) if balanced else np.ones(n)
        valid = np.isfinite(weights) & (weights > 0)
        band = np.full((2 * window + 1, n), np.nan)
        for start in range(0, n, 2000):
            stop = min(n, start + 2000)
            block = matrix.fetch(
                (chrom, start * resolution, min(size, stop * resolution)),
                (chrom, start * resolution, min(size, (stop + 2 * window) * resolution)),
            ).tocoo()
            a = block.row + start
            k = block.col - block.row
            keep = (k >= 1) & (k <= 2 * window)
            band[k[keep], a[keep]] = block.data[keep]
        # Unstored pixels between valid bins are measured zeros.
        for k in range(1, 2 * window + 1):
            row = band[k, : n - k]
            both = valid[: n - k] & valid[k:]
            row[both & np.isnan(row)] = 0.0
        score = insulation_scores(band, valid, window)
        start0 = int(bins["start"].iloc[0])
        for i, depth in boundary_calls(score, window, min_depth=min_depth):
            calls.append(
                dict(
                    chrom=chrom,
                    position0=start0 + i * resolution + resolution // 2,
                    depth=depth,
                    evidence_type="hic_insulation",
                )
            )
    if calls:
        scale = float(np.quantile([r["depth"] for r in calls], 0.95))
        for row in calls:
            row["strength"] = min(1.0, row["depth"] / scale)
    return calls
