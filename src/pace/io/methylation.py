"""CpG dyad summaries; site mean and pooled proportion are different estimands."""

import math
from bisect import bisect_left
from collections import defaultdict

from ..errors import PaceError
from .tables import integer, read_table, unique


def summarize_methylation(counts, regions, *, minimum_coverage=1, reference_cpg=None):
    minimum_coverage = integer(minimum_coverage, "methylation.minimum_coverage", minimum=1)
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
        start = integer(region["start"], "methylation region start")
        end = integer(region["end"], "methylation region end", minimum=1)
        if end <= start:
            raise PaceError("Methylation region end must exceed start")
        for sample, assay in samples:
            positions, values = indexed.get((sample, assay, region["chrom"]), ([], []))
            sub = values[bisect_left(positions, start) : bisect_left(positions, end)]
            covered = [(m, n) for _, m, n in sub if n >= minimum_coverage and n > 0]
            total = sum(n for m, n in covered)
            n_ref = None if reference_cpg is None else reference_cpg.get(region["element_id"])
            if n_ref is not None:
                n_ref = integer(n_ref, "reference CpG count")
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


def promoter_methylation_regions(promoters, *, upstream_bp=2000, downstream_bp=500):
    """Strand-aware half-open windows, including the TSS base exactly once.

    The window contains ``upstream_bp + downstream_bp + 1`` bases before
    left-edge clipping. These are annotation windows, not contact bins.
    """
    upstream_bp = integer(upstream_bp, "methylation.promoter_upstream_bp")
    downstream_bp = integer(downstream_bp, "methylation.promoter_downstream_bp")
    regions = {}
    for row in promoters:
        tss = integer(row["tss0"], "promoter tss0")
        if row["strand"] not in ("+", "-"):
            raise PaceError("Promoter methylation requires explicit + or - strand")
        left, right = (
            (upstream_bp, downstream_bp) if row["strand"] == "+" else (downstream_bp, upstream_bp)
        )
        region = {
            "element_id": row["promoter_id"],
            "chrom": row["chrom"],
            "start": max(0, tss - left),
            "end": tss + right + 1,
        }
        previous = regions.setdefault(row["promoter_id"], region)
        if previous != region:
            raise PaceError("Physical promoter_id has inconsistent methylation windows")
    return list(regions.values())


def load_reference_cpg(path):
    """Read explicit denominators for elements and physical promoter windows.

    Preferred columns are entity_type/entity_id/n_cpg. Legacy
    element_id/n_cpg files are restricted to element annotations.
    Missing rows mean an unknown denominator, never zero CpGs.
    """
    if not path:
        return {}
    result = {}
    for row in read_table(path, required=["n_cpg"]):
        kind = row.get("entity_type") or "element"
        identifier = row.get("entity_id") or row.get("element_id")
        if kind not in ("element", "promoter") or not identifier:
            raise PaceError("reference_cpg requires element/promoter entity_type and entity_id")
        key = kind, identifier
        if key in result:
            raise PaceError(f"Duplicate reference CpG denominator: {key}")
        result[key] = integer(row["n_cpg"], "n_cpg")
    return result


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
