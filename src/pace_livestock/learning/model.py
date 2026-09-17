"""Auditable JSON elastic-net logistic regression with training-only preprocessing."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from ..config import operation_config
from ..errors import PaceError
from ..evaluation.metrics import binary_metrics
from ..io.tables import number, read_table, write_table
from ..provenance import file_hash, output_directory, read_json, write_json

CORE = ["log1p_A", "log1p_C", "log1p_distance", "pace_score"]


def sigmoid(x):
    return np.exp(-np.logaddexp(0, -np.asarray(x)))


def fit_preprocessor(x):
    x = np.asarray(x, dtype=float)
    keep = np.any(np.isfinite(x), axis=0)
    if not keep.any():
        raise PaceError("No finite training features")
    used = x[:, keep]
    medians = np.nanmedian(used, axis=0)
    missing_columns = np.any(~np.isfinite(used), axis=0)
    filled = np.where(np.isfinite(used), used, medians)
    q = np.percentile(filled, [25, 75], axis=0)
    scale = q[1] - q[0]
    scale[scale == 0] = 1
    return {
        "keep": keep.tolist(),
        "medians": medians.tolist(),
        "scales": scale.tolist(),
        "missing_columns": missing_columns.tolist(),
    }


def transform(x, p):
    used = np.asarray(x, dtype=float)[:, np.array(p["keep"], dtype=bool)]
    missing = ~np.isfinite(used)
    filled = np.where(missing, np.array(p["medians"]), used)
    scaled = (filled - np.array(p["medians"])) / np.array(p["scales"])
    return np.column_stack(
        [scaled, missing[:, np.array(p["missing_columns"], dtype=bool)].astype(float)]
    )


def fit_elastic_net(
    x, y, *, lambda1=0.01, lambda2=0.01, weights=None, max_iter=10000, tolerance=1e-8
):
    """Minimize mean weighted logistic loss + lambda1*|beta| + lambda2*||beta||²/2.

    Proximal gradient uses a global Lipschitz step. The intercept is not penalized.
    No library-specific C rescaling is involved.
    """
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if x.ndim != 2 or y.shape != (len(x),) or not np.all(np.isfinite(x)) or set(y) != {0, 1}:
        raise PaceError("Classifier requires finite features and both binary classes")
    if min(lambda1, lambda2) < 0 or not all(math.isfinite(v) for v in (lambda1, lambda2)):
        raise PaceError("Elastic-net penalties must be finite and nonnegative")
    w = np.ones(len(y)) if weights is None else np.asarray(weights, dtype=float)
    if w.shape != y.shape or np.any(~np.isfinite(w)) or np.any(w < 0) or w.sum() <= 0:
        raise PaceError("Classifier sample weights must be nonnegative with positive total")
    w = w / w.sum()
    if any(w[y == c].sum() == 0 for c in (0, 1)):
        raise PaceError("Both classes require positive sample weight")
    design = np.column_stack([np.ones(len(x)), x])
    lipschitz = 0.25 * np.linalg.norm(design * np.sqrt(w[:, None]), 2) ** 2 + lambda2
    step = 1 / lipschitz
    beta = np.zeros(design.shape[1])
    converged = False
    for iteration in range(max_iter):
        gradient = design.T @ (w * (sigmoid(design @ beta) - y))
        gradient[1:] += lambda2 * beta[1:]
        new = beta - step * gradient
        new[1:] = np.sign(new[1:]) * np.maximum(np.abs(new[1:]) - step * lambda1, 0)
        if np.max(np.abs(new - beta)) <= tolerance * (1 + np.max(np.abs(beta))):
            beta, converged = new, True
            break
        beta = new
    z = design @ beta
    objective = float(
        np.sum(w * (np.logaddexp(0, z) - y * z))
        + lambda1 * np.abs(beta[1:]).sum()
        + lambda2 * np.dot(beta[1:], beta[1:]) / 2
    )
    return {
        "intercept": float(beta[0]),
        "coefficients": beta[1:].tolist(),
        "lambda1": lambda1,
        "lambda2": lambda2,
        "converged": converged,
        "iterations": iteration + 1,
        "objective": objective,
    }


def decision(x, model):
    return np.asarray(x) @ np.array(model["coefficients"]) + model["intercept"]


def feature_matrix(rows, feature_names):
    x, valid = [], []
    for r in rows:
        base = [
            number(r.get(k), k, missing=True, minimum=0)
            for k in ("A_used", "Cbar", "distance_bp", "pace_score")
        ]
        valid.append(all(math.isfinite(v) for v in base))
        if math.isfinite(base[3]) and base[3] > 1:
            raise PaceError("pace_score must be in [0,1]")
        features = dict(zip(CORE, [math.log1p(v) for v in base[:3]] + base[3:], strict=True))
        x.append(
            [
                features[name] if name in features else number(r.get(name), name, missing=True)
                for name in feature_names
            ]
        )
    return np.array(x, dtype=float).reshape(len(rows), len(feature_names)), np.array(
        valid, dtype=bool
    )


def training_labels(rows):
    usable, excluded = [], []
    for row in rows:
        status = row["label_status"]
        if str(row.get("mapping_count", "1")) != "1":
            reason = "ambiguous_region_mapping"
        elif row["effect_direction"] == "up":
            reason = "potential_repressive_or_complex"
        elif status == "enhancing_positive" and row["effect_direction"] == "down":
            usable.append({**row, "label": 1})
            continue
        elif status == "powered_negative":
            usable.append({**row, "label": 0})
            continue
        else:
            reason = "unusable_label_status"
        excluded.append({**row, "reason": reason})
    return usable, excluded


def train_classifier(config_path, out):
    cfg = operation_config(
        config_path,
        allowed={
            "data",
            "extra_features",
            "penalties",
            "folds",
            "seed",
            "model_id",
            "calibrate",
            "is_synthetic",
            "context",
        },
        required=["data", "model_id", "is_synthetic", "context"],
        paths=["data"],
    )
    rows = read_table(
        cfg["data"],
        required=[
            "element_id",
            "gene_id",
            "assayed_region_id",
            "mapping_count",
            "group_id",
            "split",
            "label_status",
            "effect_direction",
            "A_used",
            "Cbar",
            "distance_bp",
            "pace_score",
            "regime",
            "activity_sources",
            "contact_sources",
        ],
    )
    rows, excluded = training_labels(rows)
    extras = cfg.get("extra_features", [])
    if len(extras) != len(set(extras)) or set(extras) & set(CORE):
        raise PaceError("extra_features must be unique and exclude core feature names")
    names = CORE + extras
    x, valid = feature_matrix(rows, names)
    excluded += [
        {**r, "reason": "core_missing"} for r, ok in zip(rows, valid, strict=True) if not ok
    ]
    rows, x = [r for r, ok in zip(rows, valid, strict=True) if ok], x[valid]
    if not rows:
        raise PaceError("No unambiguous, core-scorable functional labels")
    groups, entities, regions = {}, {}, {}
    for r in rows:
        split = r["split"]
        if split not in ("train", "calibration", "test"):
            raise PaceError("Learning split must be train, calibration or test")
        for store, key in (
            (groups, r["group_id"]),
            (entities, (r["element_id"], r["gene_id"])),
            (regions, r["assayed_region_id"]),
        ):
            if key in store and store[key] != split:
                raise PaceError("Group, edge or perturbation region leaks across learning splits")
            store[key] = split
    train = np.array([i for i, r in enumerate(rows) if r["split"] == "train"])
    calibration = np.array(
        [i for i, r in enumerate(rows) if r["split"] == "calibration"], dtype=int
    )
    test = np.array([i for i, r in enumerate(rows) if r["split"] == "test"], dtype=int)
    y = np.array([r["label"] for r in rows])
    weights = np.array(
        [number(r.get("sample_weight", 1), "sample_weight", minimum=0) for r in rows]
    )
    if not len(train) or set(y[train]) != {0, 1}:
        raise PaceError("Training split requires both classes")
    penalties = cfg.get("penalties", [[0.01, 0.01]])
    if not penalties or any(len(p) != 2 for p in penalties):
        raise PaceError("penalties must be nonempty [lambda1,lambda2] pairs")
    cv_report = []
    chosen = penalties[0]
    if len(penalties) > 1:
        folds = cfg.get("folds", 3)
        unique_groups = sorted({rows[i]["group_id"] for i in train})
        if folds < 2 or len(unique_groups) < folds:
            raise PaceError(
                "Insufficient independent groups for tuning; predeclare a single penalty pair"
            )
        rng = np.random.default_rng(cfg.get("seed", 17))
        rng.shuffle(unique_groups)
        assignment = {g: i % folds for i, g in enumerate(unique_groups)}
        for l1, l2 in penalties:
            aps = []
            for fold in range(folds):
                tr = np.array([i for i in train if assignment[rows[i]["group_id"]] != fold])
                va = np.array([i for i in train if assignment[rows[i]["group_id"]] == fold])
                if set(y[tr]) != {0, 1} or set(y[va]) != {0, 1}:
                    raise PaceError(
                        "A grouped tuning fold lacks a class; revise grouping or freeze parameters"
                    )
                preprocessing = fit_preprocessor(x[tr])
                model = fit_elastic_net(
                    transform(x[tr], preprocessing),
                    y[tr],
                    lambda1=l1,
                    lambda2=l2,
                    weights=weights[tr],
                )
                aps.append(
                    binary_metrics(
                        y[va], sigmoid(decision(transform(x[va], preprocessing), model))
                    )["average_precision"]
                )
            cv_report.append({"lambda1": l1, "lambda2": l2, "mean_ap": float(np.mean(aps))})
        selected = max(cv_report, key=lambda r: r["mean_ap"])
        chosen = selected["lambda1"], selected["lambda2"]
    preprocessing = fit_preprocessor(x[train])
    model = fit_elastic_net(
        transform(x[train], preprocessing),
        y[train],
        lambda1=chosen[0],
        lambda2=chosen[1],
        weights=weights[train],
    )
    base_preprocessing = fit_preprocessor(x[train, : len(CORE)])
    base_model = fit_elastic_net(
        transform(x[train, : len(CORE)], base_preprocessing),
        y[train],
        lambda1=chosen[0],
        lambda2=chosen[1],
        weights=weights[train],
    )
    calibrator = None
    if cfg.get("calibrate", False):
        if len(calibration) < 4 or set(y[calibration]) != {0, 1}:
            raise PaceError(
                "Sigmoid calibration requires independent calibration data with both classes (>=4 rows)"
            )
        z = decision(transform(x[calibration], preprocessing), model)
        # Fit scalar logistic calibration using the same explicit loss with small, declared L2.
        calibrator = fit_elastic_net(z[:, None], y[calibration], lambda1=0, lambda2=1e-6)
    scope = {
        "regimes": sorted({r["regime"] for r in rows if r["split"] == "train"}),
        "activity_sources": sorted({r["activity_sources"] for r in rows if r["split"] == "train"}),
        "contact_sources": sorted({r["contact_sources"] for r in rows if r["split"] == "train"}),
    }
    asset = {
        "model_id": cfg["model_id"],
        "kind": "classifier",
        "is_synthetic": cfg["is_synthetic"],
        "context": cfg["context"],
        "feature_names": names,
        "preprocessing": preprocessing,
        "model": model,
        "calibrator": calibrator,
        "scope": scope,
        "training_data_sha256": file_hash(cfg["data"]),
        "validation": {},
        "probability_scope": "recorded_perturbation_and_candidate_sampling_distribution"
        if calibrator
        else None,
    }
    report = {
        "cv": cv_report,
        "n_train": len(train),
        "n_calibration": len(calibration),
        "n_test": len(test),
        "n_excluded": len(excluded),
        "optimization_converged": model["converged"],
        "test": binary_metrics(y[test], sigmoid(decision(transform(x[test], preprocessing), model)))
        if len(test)
        else {"reason": "no_independent_test"},
        "base_only_test": binary_metrics(
            y[test],
            sigmoid(decision(transform(x[test, : len(CORE)], base_preprocessing), base_model)),
        )
        if len(test)
        else {"reason": "no_independent_test"},
    }
    with output_directory(out) as dest:
        write_json(dest / "model.json", asset)
        write_json(
            dest / "base_only.json",
            {
                **asset,
                "feature_names": CORE,
                "preprocessing": base_preprocessing,
                "model": base_model,
                "calibrator": None,
            },
        )
        write_json(dest / "training_report.json", report)
        write_table(
            dest / "excluded_labels.tsv",
            excluded,
            fields=None if excluded else ["element_id", "gene_id", "reason"],
        )
    return report


def predict_score_rows(rows, features, model_path, *, execution_profile=None, context=None):
    path = Path(model_path)
    asset = read_json(path / "model.json" if path.is_dir() else path)
    if asset.get("kind") != "classifier":
        raise PaceError("Expected a classifier JSON artifact")
    if execution_profile != "demonstration" and asset.get("is_synthetic"):
        raise PaceError("Synthetic classifier requires demonstration profile")
    if context is not None and context != asset["context"]:
        raise PaceError("Classifier context mismatch")
    if execution_profile == "validated":
        raise PaceError("Classifier has no audited link-validation report; use research profile")
    lookup = {}
    for f in features:
        key = f["entity_type"], f["entity_id"], f["feature_name"]
        if key in lookup:
            raise PaceError(
                "ML feature input needs a declared replicate aggregation before pivoting"
            )
        lookup[key] = f["value"]
    enriched = []
    for r in rows:
        row = dict(r)
        for feature in asset["feature_names"]:
            if feature not in CORE:
                values = [
                    lookup[k]
                    for k in [
                        ("edge", f"{r['element_id']}|{r['gene_id']}", feature),
                        ("element", r["element_id"], feature),
                        ("gene", r["gene_id"], feature),
                    ]
                    if k in lookup
                ]
                if len(values) > 1:
                    raise PaceError(
                        f"Ambiguous entity level for ML feature {feature}; use distinct feature names"
                    )
                row[feature] = values[0] if values else math.nan
        enriched.append(row)
    # TSV adapters use None for missing numbers, while internal arrays use NaN.
    for row in enriched:
        for k, v in list(row.items()):
            if isinstance(v, float) and math.isnan(v):
                row[k] = None
    x, valid = feature_matrix(enriched, asset["feature_names"])
    for i, row in enumerate(rows):
        scope_ok = all(
            row[field] in asset["scope"][scope]
            for field, scope in [
                ("regime", "regimes"),
                ("activity_sources", "activity_sources"),
                ("contact_sources", "contact_sources"),
            ]
        )
        row["pace_ml_score"], row["pace_ml_probability"] = math.nan, math.nan
        row["ml_model_id"] = asset["model_id"]
        row["ml_features"] = ";".join(
            name
            for name, keep in zip(
                asset["feature_names"], asset["preprocessing"]["keep"], strict=True
            )
            if keep
        )
        row["ml_status"] = (
            "core_missing" if not valid[i] else "out_of_scope" if not scope_ok else "resolved"
        )
        if valid[i] and scope_ok:
            z = float(decision(transform(x[i : i + 1], asset["preprocessing"]), asset["model"])[0])
            row["pace_ml_score"] = float(sigmoid(z))
            if asset["calibrator"]:
                row["pace_ml_probability"] = float(sigmoid(decision([[z]], asset["calibrator"])[0]))
    return rows
