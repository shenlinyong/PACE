"""Select one promoter set per gene, shared by all its candidate elements."""

import math
from collections import defaultdict


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
