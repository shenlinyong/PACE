"""Canonical interval catalogs and sparse cis candidate generation."""

from bisect import bisect_left, bisect_right
from collections import defaultdict

from ..errors import PaceError
from ..provenance import digest


def canonical_units(
    regions: list[dict],
    chrom_sizes: dict[str, int],
    promoters=(),
    *,
    width: int = 500,
    offset: int = 0,
    include_promoters: bool = True,
):
    if width <= 0 or not 0 <= offset < width:
        raise PaceError("Grid requires width > 0 and 0 <= offset < width")
    cells, memberships, excluded = {}, set(), []
    intervals = list(regions)
    if include_promoters:
        intervals += [
            {
                "chrom": p["chrom"],
                "start": p["tss0"],
                "end": p["tss0"] + 1,
                "region_id": p["promoter_id"],
                "source_id": "promoters",
                "role": "promoter",
            }
            for p in promoters
        ]
    for region in intervals:
        chrom, start, end = region["chrom"], int(region["start"]), int(region["end"])
        if chrom not in chrom_sizes or start < 0 or end <= start or end > chrom_sizes[chrom]:
            raise PaceError(f"Region outside reference bounds: {region}")
        for index in range((start - offset) // width, (end - 1 - offset) // width + 1):
            left, right = offset + index * width, offset + (index + 1) * width
            if left < 0 or right > chrom_sizes[chrom]:
                excluded.append(
                    {"chrom": chrom, "start": left, "end": right, "reason": "incomplete_unit"}
                )
                continue
            key = chrom, left, right
            cells.setdefault(key, set()).add(region.get("role", "enhancer"))
            memberships.add((region["region_id"], f"{chrom}:{left}-{right}", region["source_id"]))
    catalog_id = digest({"width": width, "offset": offset, "cells": sorted(cells)})
    units = [
        {
            "element_id": f"{chrom}:{left}-{right}",
            "chrom": chrom,
            "start": left,
            "end": right,
            "anchor0": (left + right - 1) // 2,
            "element_roles": ";".join(sorted(cells[(chrom, left, right)])),
            "canonical_catalog_id": catalog_id,
        }
        for chrom, left, right in sorted(cells)
    ]
    membership = [
        {"region_id": r, "element_id": e, "source_id": s, "membership_rule": "overlap_at_least_1bp"}
        for r, e, s in sorted(memberships)
    ]
    return units, membership, excluded


def candidate_edges(units: list[dict], promoters: list[dict], *, radius: int = 5_000_000):
    """Binary-search local windows; memory is O(units + emitted edges), never E × G."""
    if radius < 0:
        raise PaceError("Candidate radius must be nonnegative")
    by_chrom = defaultdict(list)
    for unit in units:
        by_chrom[unit["chrom"]].append((unit["anchor0"], unit["element_id"]))
    for pairs in by_chrom.values():
        pairs.sort()
    anchors = {chrom: [x[0] for x in pairs] for chrom, pairs in by_chrom.items()}
    edges = set()
    for p in promoters:
        positions = anchors.get(p["chrom"], [])
        lo, hi = (
            bisect_left(positions, p["tss0"] - radius),
            bisect_right(positions, p["tss0"] + radius),
        )
        for _, element in by_chrom[p["chrom"]][lo:hi]:
            edges.add((element, p["gene_id"]))
    universe = digest(sorted(edges))
    return [
        {"element_id": e, "gene_id": g, "candidate_universe_id": universe} for e, g in sorted(edges)
    ]


def map_labels(labels, memberships):
    mapping = defaultdict(set)
    for row in memberships:
        mapping[row["region_id"]].add(row["element_id"])
    usable, rejected = [], []
    for row in labels:
        elements = mapping[row["assayed_region_id"]]
        status = row["label_status"]
        reason = None
        if len(elements) != 1:
            reason = "ambiguous_or_unmapped_region"
        elif row["effect_direction"] == "up":
            reason = "potential_repressive_or_complex"
        elif status not in ("enhancing_positive", "powered_negative"):
            reason = "unusable_label_status"
        elif status == "enhancing_positive" and row["effect_direction"] != "down":
            reason = "positive_requires_downregulation"
        if reason:
            rejected.append({**row, "reason": reason})
        else:
            usable.append(
                {
                    **row,
                    "element_id": next(iter(elements)),
                    "label": int(status == "enhancing_positive"),
                }
            )
    return usable, rejected
