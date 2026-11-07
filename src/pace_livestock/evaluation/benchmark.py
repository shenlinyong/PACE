"""Functional-label benchmarks with frozen baselines and explicit missed positives."""

from __future__ import annotations

import copy
import json
import math
from collections import defaultdict
from pathlib import Path

from ..catalog import map_labels
from ..config import load_config, operation_base, operation_config, strict_keys
from ..core.scoring import log_normalize, score
from ..errors import PaceError
from ..evidence.assets import load_asset
from ..io.tables import number, read_table, write_table
from ..labels import classify_label
from ..pipeline import compute
from ..provenance import (
    file_hash,
    output_directory,
    policy_digest,
    read_json,
    setting_files,
    write_json,
)
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


def verify_run_inputs(config_path, run_cfg):
    """Stop when an input of the scored run changed after the run wrote its manifest.

    The baselines are recomputed from the run's input tables, so a modified table would
    give scores that no longer belong to the reported run.
    """
    config_path = Path(config_path).resolve()
    # A result folder is passed as the directory itself; a config file sits inside one.
    folder = config_path if config_path.is_dir() else config_path.parent
    manifest_path = folder / "run_manifest.json"
    if not manifest_path.is_file():
        return
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    stale = (
        "the baselines would not describe the reported run. Rerun the analysis or restore "
        "the run's files and resolved_config.yaml"
    )
    recorded = manifest.get("input_hashes", {})
    used = {name for name, path in run_cfg["inputs"].items() if path}
    if used != set(recorded):
        changes = sorted(f"+{n}" for n in used - set(recorded)) + sorted(
            f"-{n}" for n in set(recorded) - used
        )
        raise PaceError(
            f"Run inputs differ from the ones recorded in run_manifest.json ({', '.join(changes)}); "
            + stale
        )
    for name in sorted(used):
        path = run_cfg["inputs"][name]
        if file_hash(path) != recorded[name]:
            raise PaceError(
                f"Run input {name} ({path}) differs from the hash in run_manifest.json; " + stale
            )
    # The scoring settings are part of the run: a changed partial_policy, eta or contact rule
    # gives other scores from unchanged tables.
    if manifest.get("policy_hash") and policy_digest(run_cfg) != manifest["policy_hash"]:
        raise PaceError(
            "The scoring settings differ from the ones recorded in run_manifest.json; " + stale
        )
    # Files named by settings (prior, calibrators, models): none may be added or removed,
    # and each file must be unchanged.
    if "setting_files" in manifest:
        current = setting_files(run_cfg)
        for name in sorted(set(current) | set(manifest["setting_files"])):
            before, now = manifest["setting_files"].get(name), current.get(name)
            if (before is None) != (now is None):
                change = "removed" if now is None else "added"
                raise PaceError(f"{name} was {change} after the run; " + stale)
            if before != now:
                raise PaceError(
                    f"{name} differs from the file recorded in run_manifest.json; " + stale
                )
    verify_scoring_assets(manifest, run_cfg, stale)


def verify_scoring_assets(manifest, run_cfg, stale):
    """Rebuild the scoring assets the way the run did and compare them with the manifest.

    The contact prior is a fitted asset (prior_path) or the built-in distance prior
    (prior_preset). An eQTL weak calibrator replaces its gamma and beta, so the calibrator
    file is checked first and the derived prior is rebuilt from the source prior.
    """
    from ..evidence.assets import builtin_contact_prior

    prior = None
    if run_cfg["contact"].get("prior_path"):
        prior = load_asset(run_cfg["contact"]["prior_path"], run_cfg, kind="contact_prior")
    elif run_cfg["contact"].get("prior_preset"):
        prior = builtin_contact_prior(run_cfg)
    weak = run_cfg["allocation"].get("weak_model_path")
    recorded_weak = (manifest.get("allocation") or {}).get("artifact_sha256")
    if bool(weak) != bool(recorded_weak):
        raise PaceError(
            "The eQTL weak calibrator was "
            + ("added" if weak else "removed")
            + " after the run; "
            + stale
        )
    if weak:
        if file_hash(weak) != recorded_weak:
            raise PaceError(
                f"The weak calibrator {weak} differs from the one recorded in "
                "run_manifest.json; " + stale
            )
        if prior is not None:
            from ..weak_labels import weak_calibrated_prior, weak_parameters

            prior = weak_calibrated_prior(prior, weak_parameters(read_json(weak)), weak)
    recorded = manifest.get("asset_hashes", {}).get("contact_prior")
    expected = prior["manifest_sha256"] if prior is not None else None
    if expected != recorded:
        what = (
            "was removed"
            if expected is None
            else "was added"
            if recorded is None
            else "differs from the one recorded in run_manifest.json"
        )
        raise PaceError(f"The contact prior {what}; " + stale)


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
    verify_run_inputs(cfg["run_config"], run_cfg)
    labels = read_table(
        cfg["labels"],
        required=[
            "label_id",
            "assayed_region_id",
            "gene_id",
            "species",
            "assembly",
            "context_id",
            "effect_direction",
            "label_status",
        ],
    )
    # Labels of another species or assembly can share tissue names and element IDs.
    for field in ("species", "assembly", "context_id"):
        if any(r[field] != run_cfg["context"][field] for r in labels):
            raise PaceError(f"Functional labels and scored run differ in {field}")
    membership = read_table(cfg["region_membership"], required=["region_id", "element_id"])
    from ..schemas import load_tables

    t = load_tables(run_cfg)
    catalog = {u["element_id"] for u in t["units"]}
    listed = {r["element_id"] for r in membership}
    # Elements outside the catalog are kept (their labels count as missed predictions), but a
    # membership table that shares no element with the run was built for another catalog.
    outside = listed - catalog
    if listed and not listed & catalog:
        raise PaceError(
            "Region membership shares no element with the run's catalog; map the labels "
            "onto this run's elements"
        )
    mapped, rejected = map_labels(labels, membership)
    keys = [(r["element_id"], r["gene_id"]) for r in mapped]
    if len(set(keys)) != len(keys):
        raise PaceError(
            "Repeated functional assays for one edge need a predeclared label aggregation"
        )
    # Labels whose region maps to no element, or to several, cannot be scored, but a tested
    # positive that is not found is a missed positive: it stays in the evaluation with a
    # missing score, so end-to-end recall keeps every tested positive in its denominator.
    unresolved = []
    for row in rejected:
        if row["reason"] in ("unmapped_region", "ambiguous_region"):
            label, _ = classify_label(row)
            if label is not None:
                unresolved.append(
                    {**row, "element_id": f"unresolved:{row['label_id']}", "label": label}
                )
    evaluated = mapped + unresolved
    if not evaluated:
        raise PaceError(
            "No usable functional labels (every label is excluded for its direction or status); "
            "biological benchmarking unavailable"
        )
    keys = [(r["element_id"], r["gene_id"]) for r in evaluated]
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
    # Every formula baseline uses the run's missing-value policy, so differences between
    # methods come from the model and not from which genes each one is allowed to score.
    policy = run_cfg["scoring"]["partial_policy"]
    methods["ABC_style_single_TSS"] = score(abc_edges, eta=0, partial_policy=policy)[0]
    # Distinguish allocation from simply steepening the contact exponent.
    power_edges = [{**r, "reason": ""} for r in resolved["scores"]]
    methods["contact_power_2"] = score(power_edges, eta=0, partial_policy=policy, contact_power=2)[
        0
    ]
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
    y = [r["label"] for r in evaluated]
    scores_by_method = {
        name: {key: r["pace_score"] for key, r in rows.items()} for name, rows in lookup.items()
    }
    scores_by_method["negative_distance"] = {
        (r["element_id"], r["gene_id"]): -r["distance_bp"] for r in reference_rows
    }
    # Labels every built-in method scores: compared on these, coverage plays no part.
    common_labels = [
        i
        for i, k in enumerate(keys)
        if all(math.isfinite(m.get(k, math.nan)) for m in scores_by_method.values())
    ]

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
    # Labels every method scores, external ones included: the fair set for comparisons
    # with other software. common_labels covers the built-in methods only.
    all_common = [
        i
        for i, k in enumerate(keys)
        if all(math.isfinite(m.get(k, math.nan)) for m in scores_by_method.values())
    ]
    report = []
    strata = cfg.get("stratify", [])
    for field in strata:
        if any(field not in r for r in evaluated):
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
        report.append(
            {
                "method": method,
                "scope": "common_labels",
                **binary_metrics(
                    [y[i] for i in common_labels],
                    [predictions.get(keys[i], math.nan) for i in common_labels],
                    threshold=value,
                ),
            }
        )
        if any(e["status"] == "available" for e in external_status):
            report.append(
                {
                    "method": method,
                    "scope": "common_labels_all_methods",
                    **binary_metrics(
                        [y[i] for i in all_common],
                        [predictions.get(keys[i], math.nan) for i in all_common],
                        threshold=value,
                    ),
                }
            )
        if method in common_scores:
            common = binary_metrics(y, [common_scores[method].get(k, math.nan) for k in keys])
            report.append({"method": method, "scope": "common_normalization", **common})
        for field in strata:
            for level in sorted({r[field] for r in evaluated}, key=str):
                indices = [i for i, r in enumerate(evaluated) if r[field] == level]
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
                "n_labels": len(evaluated),
                "n_labels_mapped": len(mapped),
                "n_positive_unmapped": sum(
                    r["label"] == 1 and r["reason"] == "unmapped_region" for r in unresolved
                ),
                "n_positive_ambiguous": sum(
                    r["label"] == 1 and r["reason"] == "ambiguous_region" for r in unresolved
                ),
                "n_common_labels": len(common_labels),
                "n_common_labels_all_methods": len(all_common),
                "n_membership_elements_outside_catalog": len(outside),
                "partial_policy": run_cfg["scoring"]["partial_policy"],
                "n_excluded": len(rejected),
                "input_profile": run_cfg["execution_profile"],
                "eta_calibration": calibration,
                "biological_validation_claim": "none_for_test_data"
                if run_cfg["execution_profile"] == "demonstration"
                else "task_and_dataset_specific_only",
            },
        )
    return report
