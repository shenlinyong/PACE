"""Select one promoter set per gene, shared by all its candidate elements."""

import math
from collections import defaultdict

from ..catalog import canonical_cell_bounds
from ..core import activity


def usage_weights(tables, resolved_a, cfg):
    """Set each TSS weight to its share of the gene's promoter activity.

    A TSS is measured by the activity of the grid cell that contains it. TSSs in the same
    cell are one measurable promoter cluster: the cell's activity is counted once and
    shared equally by its TSSs, so annotating more TSSs inside one cell cannot raise the
    weight of that signal. A cluster whose cell is unmeasured or absent takes the mean
    activity of the gene's measured clusters, so it neither gains nor loses weight while the
    measured clusters keep their information. Genes without measured clusters, or whose
    measured activity is all zero, use equal weights per TSS.
    """
    width, offset = cfg["catalog"]["width_bp"], cfg["catalog"]["offset_bp"]
    values = defaultdict(list)
    for r in resolved_a:
        values[r["element_id"]].append(r)
    by_gene = defaultdict(list)
    for p in tables["promoters"]:
        by_gene[p["gene_id"]].append(p)
    for group in by_gene.values():
        clusters = defaultdict(list)
        for p in group:
            start, end = canonical_cell_bounds(p["tss0"], width=width, offset=offset)
            clusters[f"{p['chrom']}:{start}-{end}"].append(p)
        usage = {}
        for cell in clusters:
            rows = values.get(cell, [])
            usage[cell] = (
                activity(
                    [r["resolved_value"] for r in rows],
                    [cfg["activity"]["pseudocounts"].get(r["assay"], 0.0) for r in rows],
                )
                if rows
                else math.nan
            )
        known = [u for u in usage.values() if math.isfinite(u)]
        mean = math.fsum(known) / len(known) if known else math.nan
        filled = {cell: u if math.isfinite(u) else mean for cell, u in usage.items()}
        total = math.fsum(filled.values())
        measured = math.isfinite(total) and total > 0
        for cell, members in clusters.items():
            for p in members:
                p["pi"] = filled[cell] / total / len(members) if measured else 1 / len(group)
                p["pi_source"] = (
                    "equal_physical_tss"
                    if not measured
                    else "promoter_activity"
                    if math.isfinite(usage[cell])
                    else "gene_mean_promoter_activity"
                )


def resolve_promoter_weights(tables, contacts, cfg):
    """Apply explicit TSS filters without changing weights separately for each element.

    The input annotation remains in the report. Filtering is disabled by default.
    A failed retained-weight threshold leaves the gene unscoreable.
    """
    settings = cfg["promoters"]
    by_gene, elements = defaultdict(list), defaultdict(set)
    for p in tables["promoters"]:
        by_gene[p["gene_id"]].append(p)
    for edge in tables["candidates"]:
        elements[edge["gene_id"]].add(edge["element_id"])
    available = {
        (r["element_id"], r["promoter_id"]) for r in contacts if math.isfinite(r["resolved_value"])
    }
    output, summaries = [], {}
    for gene, promoters in sorted(by_gene.items()):
        reasons, missing = {}, {}
        for p in promoters:
            pid = p["promoter_id"]
            missing[pid] = sum((e, pid) not in available for e in elements[gene])
            reasons[pid] = (
                "zero_weight"
                if p["pi"] == 0
                else "below_minimum_weight"
                if p["pi"] < settings["minimum_weight"]
                else "missing_contact"
                if settings["missing_policy"] == "drop_missing" and missing[pid]
                else "retained"
            )
        retained = math.fsum(p["pi"] for p in promoters if reasons[p["promoter_id"]] == "retained")
        accepted = retained > 0 and retained + 1e-12 >= settings["minimum_retained_weight"]
        changed = any(p["pi"] > 0 and reasons[p["promoter_id"]] != "retained" for p in promoters)
        status = (
            "insufficient_retained_weight"
            if not accepted
            else "filtered"
            if changed
            else "unchanged"
        )
        summaries[gene] = {
            "tss_policy_status": status,
            "tss_retained_weight": retained,
            "tss_dropped_ids": ";".join(
                sorted(
                    p["promoter_id"]
                    for p in promoters
                    if p["pi"] > 0 and reasons[p["promoter_id"]] != "retained"
                )
            ),
            "tss_contact_scope": "selected_tss_set" if changed else "all_positive_weight_tss",
        }
        for p in promoters:
            pid = p["promoter_id"]
            output.append(
                {
                    **p,
                    "pi_original": p["pi"],
                    "pi": p["pi"] / retained
                    if accepted and changed and reasons[pid] == "retained"
                    else 0.0
                    if accepted and changed
                    else p["pi"],
                    "tss_selection_reason": reasons[pid],
                    "n_missing_candidate_contacts": missing[pid],
                    **summaries[gene],
                }
            )
    return output, summaries
