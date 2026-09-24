"""Functional-label benchmarks with frozen baselines and explicit missed positives."""

from __future__ import annotations

import copy
import math
from collections import defaultdict

from ..catalog import map_labels
from ..config import load_config, operation_base, operation_config, strict_keys
from ..core.scoring import log_normalize, safe_exp, score
from ..errors import PaceError
from ..io.tables import number, read_table, write_table
from ..pipeline import compute
from ..provenance import output_directory, write_json
from .compare import load_run
from .metrics import binary_metrics


def predict_ml_command(path, out):
    from ..learning.model import predict_score_rows

    cfg = operation_config(
        path,
        allowed={"run", "model", "features"},
        required=["run", "model"],
        paths=["run", "model", "features"],
    )
    rows, manifest = load_run(cfg["run"])
    features = (
        read_table(cfg["features"], required=["entity_type", "entity_id", "feature_name", "value"])
        if cfg.get("features")
        else []
    )
    for f in features:
        f["value"] = number(f["value"], "feature value", missing=True)
    predictions = predict_score_rows(
        rows,
        features,
        cfg["model"],
        execution_profile=manifest["execution_profile"],
        context=manifest["comparison_contract"]["context"],
        feature_contract=manifest.get("ml_feature_contract"),
    )
    with output_directory(out) as dest:
        write_table(dest / "predictions.tsv", predictions)
    return predictions


def benchmark_command(path, out):
    cfg = operation_config(
        path,
        allowed={
            "run_config",
            "labels",
            "region_membership",
            "thresholds",
            "external_methods",
            "stratify",
        },
        required=["run_config", "labels", "region_membership"],
        paths=["run_config", "labels", "region_membership"],
    )
    run_cfg = load_config(cfg["run_config"])
    if run_cfg["promoters"]["weights"] != "provided":
        raise PaceError("Benchmark baselines require explicitly provided promoter weights")
    labels = read_table(
        cfg["labels"],
        required=[
            "label_id",
            "assayed_region_id",
            "gene_id",
            "context_id",
            "effect_direction",
            "label_status",
        ],
    )
    if any(r["context_id"] != run_cfg["context"]["context_id"] for r in labels):
        raise PaceError("Functional labels and scored context differ")
    membership = read_table(cfg["region_membership"], required=["region_id", "element_id"])
    mapped, rejected = map_labels(labels, membership)
    if not mapped:
        raise PaceError(
            "No evaluable unambiguous functional labels; biological benchmarking unavailable"
        )
    keys = [(r["element_id"], r["gene_id"]) for r in mapped]
    if len(set(keys)) != len(keys):
        raise PaceError(
            "Repeated functional assays for one edge need a predeclared label aggregation"
        )
    methods = {}
    calibration = None
    if run_cfg["allocation"].get("weak_model_path"):
        calibrated = compute(run_cfg)
        calibration = calibrated["eta_calibration"]
        methods["PACE_eqtl_weak"] = calibrated["scores"]
    elif run_cfg["allocation"]["labels_path"] or run_cfg["allocation"]["calibrator_path"]:
        calibrated = compute(run_cfg)
        calibration = calibrated["eta_calibration"]
        if any(
            r["gene_id"] in calibration.get("fit_genes", [])
            or r["element_id"] in calibration.get("fit_element_ids", [])
            or r.get("group_id") in calibration.get("fit_group_ids", [])
            for r in mapped
        ):
            raise PaceError(
                "Benchmark functional labels overlap eta fitting genes, elements or groups"
            )
        methods["PACE_calibrated"] = calibrated["scores"]
    elif run_cfg["allocation"]["eta"] not in ("auto", 0, 1):
        methods["PACE_fixed_eta"] = compute(run_cfg)["scores"]
    for eta in (0, 1):
        if run_cfg["allocation"].get("weak_model_path"):
            methods[f"PACE_eta{eta}"] = score(
                calibrated["scores"], eta=eta, partial_policy=run_cfg["scoring"]["partial_policy"]
            )[0]
            continue
        c = copy.deepcopy(run_cfg)
        c["allocation"].update(
            eta=eta, labels_path=None, calibrator_path=None, weak_model_path=None
        )
        methods[f"PACE_eta{eta}"] = compute(c)["scores"]
    # ABC-style single physical TSS: choose smallest tss0, then promoter_id, before evaluation.
    from ..schemas import load_tables

    t = load_tables(run_cfg)
    selected = {}
    for p in sorted(t["promoters"], key=lambda r: (r["gene_id"], r["tss0"], r["promoter_id"])):
        selected.setdefault(p["gene_id"], p["promoter_id"])
    # Reuse resolved per-TSS contacts; this is an explicit ABC-style baseline, not an external reproduction.
    base_cfg = copy.deepcopy(run_cfg)
    base_cfg["allocation"].update(
        eta=0, labels_path=None, calibrator_path=None, weak_model_path=None
    )
    resolved = calibrated if run_cfg["allocation"].get("weak_model_path") else compute(base_cfg)
    contact = {
        (r["element_id"], r["promoter_id"]): r["resolved_value"]
        for r in resolved["resolved_contacts"]
    }
    abc_edges = [
        {**r, "Cbar": contact[r["element_id"], selected[r["gene_id"]]], "reason": ""}
        for r in resolved["scores"]
    ]
    methods["ABC_style_single_TSS"] = score(abc_edges, eta=0)[0]
    # Distinguish allocation from simply steepening the contact exponent.
    power_rows = [dict(r) for r in resolved["scores"]]
    by_gene = defaultdict(list)
    for row in power_rows:
        a, contact = row["A_used"], row["Cbar"]
        row["log_support"] = (
            math.nan
            if math.isnan(a) or math.isnan(contact)
            else -math.inf
            if a == 0 or contact == 0
            else math.log(a) + 2 * math.log(contact)
        )
        row["support"] = safe_exp(row["log_support"])
        by_gene[row["gene_id"]].append(row)
    for group in by_gene.values():
        values, _ = log_normalize([r["log_support"] for r in group])
        complete = all(not math.isnan(r["log_support"]) for r in group)
        for row, value in zip(group, values, strict=True):
            row["pace_score"] = float(value) if complete else math.nan
    methods["contact_power_2"] = power_rows
    reference_rows = methods["PACE_eta0"]
    lookup = {
        name: {(r["element_id"], r["gene_id"]): r for r in rows} for name, rows in methods.items()
    }
    # Intersect scoreable support sets per gene across formula methods, and renormalize each.
    gene_elements = defaultdict(list)
    for row in reference_rows:
        key = row["element_id"], row["gene_id"]
        if all(
            key in method and not math.isnan(method[key]["log_support"])
            for method in lookup.values()
        ):
            gene_elements[row["gene_id"]].append(row["element_id"])
    common_scores = {name: {} for name in methods}
    for gene, elements in gene_elements.items():
        if len(elements) < 2:
            continue
        normalizations = {
            name: log_normalize([method[e, gene]["log_support"] for e in elements])
            for name, method in lookup.items()
        }
        if not all(math.isfinite(total) for _, total in normalizations.values()):
            continue
        for name, (scores, _) in normalizations.items():
            common_scores[name].update(
                {(e, gene): float(s) for e, s in zip(elements, scores, strict=True)}
            )
    y = [r["label"] for r in mapped]
    scores_by_method = {
        name: {key: r["pace_score"] for key, r in rows.items()} for name, rows in lookup.items()
    }
    scores_by_method["negative_distance"] = {
        (r["element_id"], r["gene_id"]): -r["distance_bp"] for r in reference_rows
    }

    external_status = []
    for external in cfg.get("external_methods", []):
        strict_keys(external, {"name", "path", "version", "configuration"}, "external_method")
        if not external.get("version") or not external.get("configuration"):
            raise PaceError(
                "External benchmark method requires version and configuration provenance"
            )
        if external["name"] in scores_by_method:
            raise PaceError("External method name collides with a built-in baseline")
        ep = (operation_base(path) / external["path"]).resolve()
        if not ep.is_file():
            external_status.append(
                {"method": external["name"], "status": "not_available", "reason": "scores_missing"}
            )
            continue
        external_rows = read_table(ep, required=["element_id", "gene_id", "score"])
        ek = [(r["element_id"], r["gene_id"]) for r in external_rows]
        if len(ek) != len(set(ek)):
            raise PaceError("Duplicate external method edge")
        scores_by_method[external["name"]] = {
            key: number(r["score"], "external score", missing=True)
            for key, r in zip(ek, external_rows, strict=True)
        }
        external_status.append(
            {
                **external,
                "status": "available",
                "common_denominator": "not_assessed_without_support",
            }
        )
    report = []
    strata = cfg.get("stratify", [])
    for field in strata:
        if any(field not in r for r in mapped):
            raise PaceError(f"Unknown label stratification field: {field}")
    for method, predictions in scores_by_method.items():
        threshold = cfg.get("thresholds", {}).get(method)
        if threshold:
            strict_keys(threshold, {"value", "source_split", "source_id"}, "threshold")
            if threshold.get("source_split") not in ("train", "calibration") or not threshold.get(
                "source_id"
            ):
                raise PaceError(
                    "Decision thresholds must have recorded train/calibration provenance"
                )
        value = threshold["value"] if threshold else None
        metric = binary_metrics(y, [predictions.get(k, math.nan) for k in keys], threshold=value)
        report.append({"method": method, "scope": "all_tested_labels", **metric})
        if method in common_scores:
            common = binary_metrics(y, [common_scores[method].get(k, math.nan) for k in keys])
            report.append({"method": method, "scope": "common_normalization", **common})
        for field in strata:
            for level in sorted({r[field] for r in mapped}, key=str):
                indices = [i for i, r in enumerate(mapped) if r[field] == level]
                report.append(
                    {
                        "method": method,
                        "scope": f"{field}={level}",
                        **binary_metrics(
                            [y[i] for i in indices],
                            [predictions.get(keys[i], math.nan) for i in indices],
                            threshold=value,
                        ),
                    }
                )
    with output_directory(out) as dest:
        write_table(dest / "metrics.tsv", report)
        write_table(
            dest / "excluded_labels.tsv",
            rejected,
            fields=None if rejected else ["label_id", "reason"],
        )
        write_json(
            dest / "benchmark_manifest.json",
            {
                "ABC_TSS_rule": "minimum_tss0_then_promoter_id",
                "selected_promoters": selected,
                "average_precision_definition": "sum(recall increments * precision), grouped by score ties; not trapezoidal PR area",
                "external_methods": external_status,
                "thresholds": cfg.get("thresholds", {}),
                "n_labels": len(mapped),
                "n_excluded": len(rejected),
                "input_profile": run_cfg["execution_profile"],
                "eta_calibration": calibration,
                "biological_validation_claim": "none_for_synthetic_data"
                if run_cfg["execution_profile"] == "demonstration"
                else "task_and_dataset_specific_only",
            },
        )
    return report
