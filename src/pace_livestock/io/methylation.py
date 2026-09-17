"""CpG dyad summaries; site mean and pooled proportion are different estimands."""

import math
from bisect import bisect_left
from collections import defaultdict

from ..errors import PaceError
from .tables import integer, unique


def summarize_methylation(counts, regions, *, minimum_coverage=1, reference_cpg=None):
    unique(counts, ("chrom", "dyad_start0", "sample_id", "assay"), "CpG dyads (already merged)")
    groups = defaultdict(list)
    for row in counts:
        pos = integer(row["dyad_start0"], "dyad_start0")
        m, n = (
            integer(row["methylated_count"], "methylated_count"),
            integer(row["total_count"], "total_count"),
        )
        if m > n:
            raise PaceError("CpG methylated_count cannot exceed total_count")
        groups[row["sample_id"], row["assay"], row["chrom"]].append((pos, m, n))
    indexed = {}
    for key, values in groups.items():
        values.sort()
        indexed[key] = ([v[0] for v in values], values)
    samples = sorted({(r["sample_id"], r["assay"]) for r in counts})
    output = []
    for region in regions:
        for sample, assay in samples:
            positions, values = indexed.get((sample, assay, region["chrom"]), ([], []))
            sub = values[
                bisect_left(positions, region["start"]) : bisect_left(positions, region["end"])
            ]
            covered = [(m, n) for _, m, n in sub if n >= minimum_coverage and n > 0]
            total = sum(n for m, n in covered)
            n_ref = None if reference_cpg is None else reference_cpg.get(region["element_id"])
            if n_ref is not None and n_ref < len(sub):
                raise PaceError("Observed CpG count exceeds supplied reference CpG count")
            status = (
                "observed"
                if total
                else "no_cpg"
                if n_ref == 0
                else "low_coverage"
                if sub
                else "uncovered"
            )
            output.append(
                {
                    "element_id": region["element_id"],
                    "sample_id": sample,
                    "assay": assay,
                    "M_site": sum(m / n for m, n in covered) / len(covered)
                    if covered
                    else math.nan,
                    "M_pooled": sum(m for m, n in covered) / total if total else math.nan,
                    "covered_cpg": len(covered),
                    "reported_cpg": len(sub),
                    "total_reads": total,
                    "reference_cpg": n_ref,
                    "cpg_coverage_fraction": len(covered) / n_ref if n_ref else math.nan,
                    "status": status,
                }
            )
    return output


def merge_stranded_cpg(rows, reference):
    """Convert cytosine calls with explicit strand and 0-based position exactly once."""
    unique(rows, ("chrom", "pos0", "strand", "sample_id", "assay"), "stranded CpG")
    merged = defaultdict(lambda: [0, 0])
    for row in rows:
        if row["strand"] not in ("+", "-"):
            raise PaceError("CpG conversion requires a + or - strand")
        pos = integer(row["pos0"], "CpG pos0") - (row["strand"] == "-")
        if pos < 0 or reference.fetch(row["chrom"], pos, pos + 2).upper() != "CG":
            raise PaceError("CpG dyad does not match reference CG")
        m, n = (
            integer(row["methylated_count"], "methylated_count"),
            integer(row["total_count"], "total_count"),
        )
        if m > n:
            raise PaceError("CpG methylated count exceeds total count")
        key = row["chrom"], pos, row["sample_id"], row["assay"]
        merged[key][0] += m
        merged[key][1] += n
    return [
        {
            "chrom": chrom,
            "dyad_start0": pos,
            "sample_id": sample,
            "assay": assay,
            "methylated_count": m,
            "total_count": n,
        }
        for (chrom, pos, sample, assay), (m, n) in sorted(merged.items())
    ]
