"""Sparse cooler bin-pair queries; invalid balanced bins stay unavailable."""

from collections import defaultdict

import numpy as np

from ..errors import PaceError


def query_contacts(
    uri,
    pairs,
    *,
    resolution: int,
    balanced: bool,
    missing_pixels_are_zero: bool,
    diagonal_window_bp=5000,
):
    if type(balanced) is not bool or type(missing_pixels_are_zero) is not bool:
        raise PaceError("balanced and missing_pixels_are_zero must be YAML booleans")
    try:
        import cooler
    except ImportError as exc:
        raise PaceError("cool/mcool support requires pip install 'pace-livestock[io]'") from exc
    c = cooler.Cooler(str(uri))
    if c.binsize != resolution:
        raise PaceError(f"cooler resolution {c.binsize} does not match requested {resolution}")
    if balanced and "weight" not in c.bins().columns:
        raise PaceError("Balanced contact requested but cooler has no weight column")
    bins = c.bins()[:]
    valid = (
        np.isfinite(bins["weight"].to_numpy()) & (bins["weight"].to_numpy() > 0)
        if balanced
        else np.ones(len(bins), dtype=bool)
    )
    offsets = {chrom: c.offset(chrom) for chrom in c.chromnames}
    requests, grouped = [], defaultdict(set)
    for pair in pairs:
        chrom = pair["chrom"]
        if (
            chrom not in offsets
            or min(pair["anchor0"], pair["tss0"]) < 0
            or max(pair["anchor0"], pair["tss0"]) >= c.chromsizes[chrom]
        ):
            raise PaceError("Contact pair coordinates do not match cooler reference")
        i, j = sorted(
            (
                offsets[chrom] + pair["anchor0"] // resolution,
                offsets[chrom] + pair["tss0"] // resolution,
            )
        )
        requests.append((pair, i, j))
        grouped[i].add(j)
    neighbors = {}
    for pair, i, j in requests:
        if i != j:
            continue
        chrom = pair["chrom"]
        first, last = c.extent(chrom)
        radius = max(1, int(np.ceil(diagonal_window_bp / resolution)))
        near = [
            tuple(sorted((i, k)))
            for k in range(max(first, i - radius), min(last, i + radius + 1))
            if k != i
        ]
        neighbors[i] = near
        for x, y in near:
            grouped[x].add(y)
    lookup = {}
    selector = c.matrix(balance=balanced, sparse=True)
    for i, js in grouped.items():
        low, high = min(js), max(js) + 1
        block = selector[i : i + 1, low:high].tocoo()
        stored = {
            low + int(j): float(value) for j, value in zip(block.col, block.data, strict=True)
        }
        for j in js:
            value = stored.get(j, 0.0 if missing_pixels_are_zero else np.nan)
            if not valid[i] or not valid[j]:
                value = np.nan
            if np.isfinite(value) and value < 0:
                raise PaceError("Negative contact value")
            lookup[i, j] = value
    return [
        {
            **pair,
            "contact_value": lookup[i, j],
            "bin_pair_id": f"{i}:{j}",
            "resolution": resolution,
            "measurement_status": "observed" if np.isfinite(lookup[i, j]) else "unmappable",
            "near_diagonal_value": max(
                [lookup[x, y] for x, y in neighbors.get(i, []) if np.isfinite(lookup[x, y])],
                default=np.nan,
            )
            if i == j and valid[i]
            else np.nan,
            "near_diagonal_method": "neighbor_max" if i == j else None,
            "near_diagonal_source_bins": ";".join(f"{x}:{y}" for x, y in neighbors.get(i, []))
            if i == j
            else None,
        }
        for pair, i, j in requests
    ]
