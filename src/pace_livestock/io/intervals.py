"""Measured BED overlap annotations, with union coverage and optional motif orientation."""

from bisect import bisect_left
from collections import defaultdict

from ..errors import PaceError


def interval_features(
    regions, peaks, *, evidence_id, feature_prefix, entity_type="element", motif_strands=False
):
    """Annotate supplied regions; callers must declare the BED assay's measurement scope.

    For E–P span annotations, use a frozen span as the region and entity_type='edge'.
    Orientation is descriptive and is emitted only for explicit +/- motif records.
    """
    groups = defaultdict(list)
    for peak in peaks:
        if peak["start"] < 0 or peak["end"] <= peak["start"]:
            raise PaceError("Invalid BED peak interval")
        if motif_strands and peak.get("strand") not in ("+", "-"):
            raise PaceError("Motif orientation requires explicit + or - strand")
        groups[peak["chrom"]].append(peak)
    index = {}
    for chrom, rows in groups.items():
        rows.sort(key=lambda r: (r["start"], r["end"]))
        ends, maximum = [], 0
        for r in rows:
            maximum = max(maximum, r["end"])
            ends.append(maximum)
        index[chrom] = (rows, [r["start"] for r in rows], ends)
    output = []
    for region in regions:
        left, right = int(region["start"]), int(region["end"])
        if left < 0 or right <= left:
            raise PaceError("Invalid annotation region")
        rows, starts, ends = index.get(region["chrom"], ([], [], []))
        low, high = bisect_left(ends, left + 1), bisect_left(starts, right)
        overlaps = [r for r in rows[low:high] if r["end"] > left]
        intervals = sorted({(max(left, r["start"]), min(right, r["end"])) for r in overlaps})
        merged = []
        for a, b in intervals:
            if merged and a <= merged[-1][1]:
                merged[-1] = merged[-1][0], max(merged[-1][1], b)
            else:
                merged.append((a, b))
        values = {
            "peak_count": len({(r["start"], r["end"]) for r in overlaps}),
            "overlap_fraction": sum(b - a for a, b in merged) / (right - left),
        }
        if motif_strands:
            values.update(
                {
                    "plus_motifs": sum(r["strand"] == "+" for r in overlaps),
                    "minus_motifs": sum(r["strand"] == "-" for r in overlaps),
                }
            )
        for feature, value in values.items():
            output.append(
                {
                    "entity_type": entity_type,
                    "entity_id": region["element_id"],
                    "feature_name": f"{feature_prefix}:{feature}",
                    "value": value,
                    "evidence_id": evidence_id,
                    "status": "observed",
                }
            )
    return output
