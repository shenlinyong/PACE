"""Common-denominator comparison from support, never subtraction after an inner join."""

import math
from collections import defaultdict
from pathlib import Path

import numpy as np

from ..config import operation_config
from ..core.scoring import log_normalize, safe_exp
from ..errors import PaceError
from ..io.tables import number, read_table, write_table
from ..provenance import digest, output_directory, read_json, write_json


def load_run(path):
    path = Path(path)
    rows = read_table(
        path / "scores.tsv.gz",
        required=["element_id", "gene_id", "log_support", "support", "pace_score"],
    )
    for r in rows:
        for field in (
            "support",
            "log_support",
            "pace_score",
            "A_used",
            "Cbar",
            "distance_bp",
            "denominator",
            "log_denominator",
        ):
            if field in r:
                r[field] = number(r[field], field, missing=True)
        # Finite logarithms are authoritative: exp(log_support) may underflow to
        # zero even though the relative support remains perfectly measurable.
        # Genuine zeros lose -inf in TSV and are recovered from an explicit state.
        if (
            math.isnan(r["log_support"])
            and r["support"] == 0
            and (r.get("support_status") == "zero_support" or r.get("reason") == "zero_support")
        ):
            r["log_support"] = -math.inf
    return rows, read_json(path / "run_manifest.json")


def compare_rows(left, right, *, minimum_common=2, full_compatible=True, absolute_compatible=True):
    if minimum_common < 2:
        raise PaceError("Compositional comparison requires at least two common units")
    a = {(r["gene_id"], r["element_id"]): r for r in left}
    b = {(r["gene_id"], r["element_id"]): r for r in right}
    if len(a) != len(left) or len(b) != len(right):
        raise PaceError("Comparison input has duplicate candidate edges")
    by_gene = defaultdict(list)
    for gene, element in sorted(set(a) | set(b)):
        by_gene[gene].append(element)
    result = []
    for gene, elements in sorted(by_gene.items()):
        common = [
            e
            for e in elements
            if (gene, e) in a
            and (gene, e) in b
            and not math.isnan(a[gene, e]["log_support"])
            and not math.isnan(b[gene, e]["log_support"])
        ]
        pa, ta = log_normalize([a[gene, e]["log_support"] for e in common])
        pb, tb = log_normalize([b[gene, e]["log_support"] for e in common])
        informative = len(common) >= minimum_common and math.isfinite(ta) and math.isfinite(tb)
        common_a, common_b = dict(zip(common, pa, strict=True)), dict(zip(common, pb, strict=True))
        total_a = log_normalize([a[gene, e]["log_support"] for e in elements if (gene, e) in a])[1]
        total_b = log_normalize([b[gene, e]["log_support"] for e in elements if (gene, e) in b])[1]
        complete = len(common) == len(elements) and all(
            a[gene, e].get("normalization_status")
            == b[gene, e].get("normalization_status")
            == "complete"
            for e in common
        )
        structural_ok = all(
            not any(
                tag in str(row.get("structural_status", ""))
                for tag in ("unsupported", "affects_relationship")
            )
            for e in common
            for row in (a[gene, e], b[gene, e])
        )
        universe = digest({"gene": gene, "common_elements": common})
        for element in elements:
            ra, rb = a.get((gene, element), {}), b.get((gene, element), {})
            original_a, original_b = ra.get("pace_score", math.nan), rb.get("pace_score", math.nan)
            ca, cb = common_a.get(element, math.nan), common_b.get(element, math.nan)
            conditional = cb - ca if informative else math.nan
            full = (
                original_b - original_a
                if informative and complete and full_compatible and structural_ok
                else math.nan
            )
            reason = (
                "complete"
                if math.isfinite(full)
                else "insufficient_common_units"
                if len(common) < minimum_common
                else "zero_common_denominator"
                if not informative
                else "parameter_difference"
                if not full_compatible
                else "technical_missing_or_structure"
            )
            da = (
                rb.get("A_used", math.nan) - ra.get("A_used", math.nan)
                if absolute_compatible
                else math.nan
            )
            ds = (
                rb.get("support", math.nan) - ra.get("support", math.nan)
                if absolute_compatible
                else math.nan
            )
            result.append(
                {
                    "gene_id": gene,
                    "element_id": element,
                    "original_score_a": original_a,
                    "original_score_b": original_b,
                    "common_score_a": ca if informative else math.nan,
                    "common_score_b": cb if informative else math.nan,
                    "full_delta_pace": full,
                    "conditional_delta_pace": conditional,
                    "delta_A": da,
                    "delta_support": ds,
                    "gene_total_support_a": safe_exp(total_a),
                    "gene_total_support_b": safe_exp(total_b),
                    "common_denominator_a": safe_exp(ta),
                    "common_denominator_b": safe_exp(tb),
                    "n_common_units": len(common),
                    "comparison_universe_id": universe,
                    "reason": reason,
                }
            )
    return result


def compare_runs(
    left_path,
    right_path,
    *,
    allow_eta_difference=False,
    allow_evidence_difference=False,
    minimum_common=2,
):
    a, ma = load_run(left_path)
    b, mb = load_run(right_path)
    ca, cb = ma["comparison_contract"], mb["comparison_contract"]
    if any("contact_measurement_contract" not in contract for contract in (ca, cb)):
        raise PaceError(
            "Comparison requires a recorded contact measurement contract; rerun legacy results"
        )
    different = {k for k in set(ca) | set(cb) if ca.get(k) != cb.get(k)}
    allowed = {"eta"} if allow_eta_difference else set()
    if allow_evidence_difference:
        allowed.add("evidence_policy")
    if different - allowed:
        raise PaceError(f"Incompatible comparison contracts: {sorted(different - allowed)}")
    return (
        compare_rows(
            a,
            b,
            minimum_common=minimum_common,
            full_compatible=not different,
            absolute_compatible=not different,
        ),
        ma,
        mb,
    )


def compare_command(config_path, out):
    cfg = operation_config(
        config_path,
        allowed={
            "left",
            "right",
            "minimum_common_units",
            "allow_eta_difference",
            "allow_evidence_difference",
        },
        required=["left", "right"],
        paths=["left", "right"],
    )
    rows, ma, mb = compare_runs(
        cfg["left"],
        cfg["right"],
        minimum_common=cfg.get("minimum_common_units", 2),
        allow_eta_difference=cfg.get("allow_eta_difference", False),
        allow_evidence_difference=cfg.get("allow_evidence_difference", False),
    )
    with output_directory(out) as dest:
        write_table(dest / "comparison.tsv", rows)
        write_json(
            dest / "comparison_manifest.json",
            {
                "left": ma["run_id"],
                "right": mb["run_id"],
                "left_config_hash": ma["config_hash"],
                "right_config_hash": mb["config_hash"],
                "interpretation": "Conditional deltas refer only to the common measurable subset.",
            },
        )
    return rows


def stability_command(config_path, out):
    cfg = operation_config(
        config_path,
        allowed={"replicates", "minimum_common_units"},
        required=["replicates"],
        paths=["replicates"],
    )
    replicates = read_table(cfg["replicates"], required=["run_path", "donor_id", "replicate_type"])
    if len(replicates) < 2:
        raise PaceError("stability needs at least two runs")
    base = Path(cfg["replicates"]).parent
    output = []
    for i, left in enumerate(replicates):
        for right in replicates[i + 1 :]:
            if left["replicate_type"] not in ("biological", "technical") or right[
                "replicate_type"
            ] not in ("biological", "technical"):
                raise PaceError("replicate_type must be biological or technical")
            rows, _, _ = compare_runs(
                base / left["run_path"],
                base / right["run_path"],
                minimum_common=cfg.get("minimum_common_units", 2),
            )
            x = np.array([r["common_score_a"] for r in rows])
            y = np.array([r["common_score_b"] for r in rows])
            valid = np.isfinite(x) & np.isfinite(y)
            correlation = (
                float(np.corrcoef(x[valid], y[valid])[0, 1])
                if valid.sum() >= 2 and np.std(x[valid]) > 0 and np.std(y[valid]) > 0
                else math.nan
            )
            output.append(
                {
                    "run_a": left["run_path"],
                    "run_b": right["run_path"],
                    "pearson": correlation,
                    "n_common_edges": int(valid.sum()),
                    "conditional": any(r["reason"] != "complete" for r in rows),
                    "independence": "same_donor"
                    if left["donor_id"] == right["donor_id"]
                    else "distinct_donors",
                    "reason": "resolved"
                    if math.isfinite(correlation)
                    else "constant_or_insufficient",
                }
            )
    with output_directory(out) as dest:
        write_table(dest / "stability.tsv", output)
        write_json(
            dest / "report.json",
            {
                "n_independent_donors": len({r["donor_id"] for r in replicates}),
                "confidence_intervals": "not_estimated",
                "technical_replicates_are_not_independent_animals": True,
            },
        )
    return output
