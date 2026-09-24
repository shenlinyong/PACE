"""Chromosome-separated eQTL calibration; weak evidence is not perturbation truth."""

from __future__ import annotations

import copy
import math
from bisect import bisect_left
from collections import defaultdict
from itertools import product

import numpy as np

from .boundary_prior import parameter_grid
from .core import score
from .errors import PaceError
from .io.tables import integer, number
from .provenance import clean, digest, file_hash, read_json


def map_weak_labels(variants, units, candidates, *, aggregation="independent_variants"):
    """Map zero-based variant positions to half-open elements without duplicating PIPs.

    independent_variants uses 1-prod(1-PIP), an explicit independence approximation.
    independent_signals sums mutually exclusive variants within each signal first.
    Candidates for genes present in the fine-mapping table form the assessed universe;
    zero mass is unlabelled background, not an experimentally tested negative.
    """
    if aggregation not in ("independent_variants", "independent_signals"):
        raise PaceError("Unknown PIP aggregation")
    by_gene, seen, signals = defaultdict(list), set(), defaultdict(float)
    unit_map = {r["element_id"]: r for r in units}
    candidate_keys = {(r["element_id"], r["gene_id"]) for r in candidates}
    gene_chrom = defaultdict(set)
    for element, gene in candidate_keys:
        gene_chrom[gene].add(unit_map[element]["chrom"])
    for row in variants:
        gene, chrom, variant = row.get("gene_id"), row.get("chrom"), row.get("variant_id")
        if not gene or not variant or gene not in gene_chrom or gene_chrom[gene] != {chrom}:
            raise PaceError("eQTL genes/coordinates must match the cis candidate catalog")
        position = integer(row.get("pos0"), "variant pos0")
        pip = number(row.get("pip"), "PIP", minimum=0, maximum=1)
        signal = row.get("signal_id") if aggregation == "independent_signals" else variant
        if not signal:
            raise PaceError("independent_signals requires signal_id")
        key = gene, variant
        if key in seen:
            raise PaceError("Duplicate gene/variant PIP; supply joint fine-mapping membership once")
        seen.add(key)
        record = dict(gene_id=gene, chrom=chrom, pos0=position, pip=pip, signal_id=signal)
        by_gene[gene].append(record)
        signals[gene, signal] += pip
    if aggregation == "independent_signals" and any(v > 1 + 1e-8 for v in signals.values()):
        raise PaceError("PIPs within a gene/signal must sum to at most one")
    gene_index = {}
    for gene, rows in by_gene.items():
        rows.sort(key=lambda r: r["pos0"])
        gene_index[gene] = [r["pos0"] for r in rows]
    # Binary searches restrict each interval query to variants inside the element.
    output = []
    for element, gene in sorted(candidate_keys):
        if gene not in by_gene:
            continue
        unit = unit_map[element]
        mass = defaultdict(float)
        positions = gene_index[gene]
        low, high = bisect_left(positions, unit["start"]), bisect_left(positions, unit["end"])
        for variant in by_gene[gene][low:high]:
            mass[variant["signal_id"]] += variant["pip"]
        values = [min(1.0, x) for x in mass.values()]
        q = 1 - math.prod(1 - v for v in values)
        output.append(
            dict(
                element_id=element,
                gene_id=gene,
                chrom=unit["chrom"],
                weak_label=q,
                n_signals=len(values),
                label_status="pip_mass" if q else "unlabelled_background",
            )
        )
    if not output or not any(r["weak_label"] > 0 for r in output):
        raise PaceError("No fine-mapping PIP overlaps the candidate elements")
    return output


def soft_average_precision(labels, predictions):
    """Tie-aware fractional AP, a descriptive ranking metric for soft PIP mass."""
    q, p = np.asarray(labels, float), np.asarray(predictions, float)
    if (
        q.shape != p.shape
        or q.ndim != 1
        or not len(q)
        or not np.all(np.isfinite(q))
        or not np.all(np.isfinite(p))
        or np.any(q < 0)
        or np.any(q > 1)
    ):
        raise PaceError("Soft AP needs matching finite predictions and PIPs in [0, 1]")
    if q.sum() <= 0:
        raise PaceError("Every evaluation chromosome needs positive PIP mass")
    order = np.argsort(-p, kind="stable")
    q, p = q[order], p[order]
    ends = np.r_[np.flatnonzero(p[1:] != p[:-1]) + 1, len(p)]
    cumulative = np.cumsum(q)
    positives = np.diff(np.r_[0.0, cumulative[ends - 1]])
    return float(np.sum(positives * cumulative[ends - 1] / ends) / q.sum())


def fit_chromosome_grid(labels, predictions, *, test_chromosomes, baselines=None):
    """Nested leave-one-chromosome-out selection, then untouched final test chromosomes."""
    chroms = sorted({r["chrom"] for r in labels})
    if (
        not isinstance(test_chromosomes, list)
        or not test_chromosomes
        or len(set(test_chromosomes)) != len(test_chromosomes)
    ):
        raise PaceError("Declare distinct final test_chromosomes before fitting")
    tests = set(test_chromosomes)
    if not tests <= set(chroms):
        raise PaceError("Test chromosomes must have assessed eQTL candidates")
    train = sorted(set(chroms) - tests)
    if len(train) < 3:
        raise PaceError(
            "Weak calibration needs at least three training chromosomes plus final test"
        )
    indices = {c: [i for i, r in enumerate(labels) if r["chrom"] == c] for c in chroms}
    metrics = {}
    for params, values in predictions.items():
        if len(values) != len(labels):
            raise PaceError("Grid predictions do not match label universe")
        metrics[params] = {
            c: soft_average_precision(
                [labels[i]["weak_label"] for i in ix], [values[i] for i in ix]
            )
            for c, ix in indices.items()
        }
    if not metrics:
        raise PaceError("No parameter combinations were supplied")

    def select(chromosomes):
        # Exact ties retain the simpler eta=0, beta=0 model.
        return min(
            metrics,
            key=lambda p: (-float(np.mean([metrics[p][c] for c in chromosomes])), p[2], p[1], p[0]),
        )

    folds = []
    for held in train:
        selected = select([c for c in train if c != held])
        folds.append(
            dict(
                test_chromosome=held,
                training_chromosomes=[c for c in train if c != held],
                selected=dict(zip(("gamma", "beta", "eta"), selected, strict=True)),
                soft_ap=metrics[selected][held],
            )
        )
    best = select(train)
    report = dict(
        selected=dict(zip(("gamma", "beta", "eta"), best, strict=True)),
        training_chromosomes=train,
        test_chromosomes=sorted(tests),
        nested_folds=folds,
        nested_macro_soft_ap=float(np.mean([r["soft_ap"] for r in folds])),
        test_soft_ap={c: metrics[best][c] for c in sorted(tests)},
        grid_results=[
            dict(
                gamma=p[0],
                beta=p[1],
                eta=p[2],
                training_macro_soft_ap=float(np.mean([m[c] for c in train])),
            )
            for p, m in sorted(metrics.items())
        ],
        metric="fractional_AP_of_PIP_mass",
        functional_validation=False,
    )
    report["baselines"] = {}
    for name, values in (baselines or {}).items():
        if len(values) != len(labels):
            raise PaceError("Baseline predictions do not match label universe")
        report["baselines"][name] = {
            c: soft_average_precision(
                [labels[i]["weak_label"] for i in ix], [values[i] for i in ix]
            )
            for c, ix in indices.items()
        }
    return report


def model_scope(cfg, tables, assets, resolved_activity):
    from .learning.allocation import calibration_scope

    scales = sorted(
        {
            (r["assay"], r["unit"], r.get("normalization_id"), r["window_id"])
            for r in resolved_activity
            if r["unit"] is not None
        },
        key=str,
    )
    return clean(calibration_scope(cfg, tables, assets, scales))


def apply_weak_model(cfg, tables, assets, resolved_activity):
    path = cfg["allocation"]["weak_model_path"]
    model = read_json(path)
    if cfg["execution_profile"] == "validated" or model.get("functional_validation") is not False:
        raise PaceError("eQTL weak evidence cannot establish functional validation")
    if model.get("schema_version") != "pace-weak-1" or model.get("kind") != "eqtl_weak_calibrator":
        raise PaceError("Not a PACE eQTL weak calibrator")
    if model.get("scope") != model_scope(cfg, tables, assets, resolved_activity):
        raise PaceError("Weak calibrator context, catalog, policy or source prior differs")
    if (
        model.get("status") != "weak_fitted"
        or len(model.get("training_chromosomes", [])) < 3
        or not model.get("test_chromosomes")
    ):
        raise PaceError("Weak calibrator lacks chromosome-separated fitting evidence")
    if set(model["training_chromosomes"]) & set(model["test_chromosomes"]):
        raise PaceError("Weak calibrator training/test chromosomes overlap")
    parameters = model["selected"]
    for key in ("gamma", "beta", "eta"):
        value = number(parameters[key], key, minimum=0, maximum=1 if key == "eta" else None)
        if key == "gamma" and value == 0:
            raise PaceError("gamma must be positive")
        parameters[key] = value
    prior = {
        **assets["contact_prior"],
        "gamma": parameters["gamma"],
        "beta": parameters["beta"],
        "validation": {},
        "weak_calibrator_sha256": file_hash(path),
    }
    prior["calibration_sources"] = [*prior["calibration_sources"], file_hash(path)]
    prior["model_id"] += ":weak:" + file_hash(path)[:12]
    prior["manifest_sha256"] = digest({k: v for k, v in prior.items() if k != "asset_directory"})
    allocation = {
        **model,
        "eta": parameters["eta"],
        "reuse": True,
        "reason": "eQTL PIP ranking; not perturbation validation",
        "artifact_sha256": file_hash(path),
    }
    return prior, allocation


def calibrate(
    cfg,
    variants,
    *,
    gamma_grid,
    beta_grid,
    eta_grid,
    test_chromosomes,
    aggregation,
    abc_gamma=1.0242386,
):
    from .evidence.assets import builtin_contact_prior, load_asset
    from .evidence.resolve import resolve_activity
    from .pipeline import compute
    from .schemas import load_tables

    if cfg["allocation"]["eta"] != "auto" or any(
        cfg["allocation"][k] for k in ("labels_path", "calibrator_path", "weak_model_path")
    ):
        raise PaceError("fit-labels needs eta=auto and no existing calibration inputs")
    if cfg["execution_profile"] == "validated" or cfg["contact"]["mode"] not in (
        "prior_only",
        "shrinkage",
    ):
        raise PaceError("Weak fitting requires research/demonstration prior_only or shrinkage")
    if cfg["inputs"]["resolved_contacts"]:
        raise PaceError("Refit from observed contacts, not imported resolved contacts")
    prior = (
        load_asset(cfg["contact"]["prior_path"], cfg, kind="contact_prior")
        if cfg["contact"]["prior_path"]
        else builtin_contact_prior(cfg)
        if cfg["contact"]["prior_preset"]
        else None
    )
    if prior is None:
        raise PaceError("fit-labels requires an explicit source prior")
    if (
        cfg["contact"]["reliability"] == "per_pair"
        and cfg["contact"]["kappa"] == "auto"
        and prior.get("kappa") is None
    ):
        raise PaceError(
            "Freeze kappa in the source prior or run config before chromosome validation"
        )
    tables = load_tables(cfg)
    scope = model_scope(cfg, tables, {"contact_prior": prior}, resolve_activity(tables, cfg))
    labels = map_weak_labels(
        variants, tables["units"], tables["candidates"], aggregation=aggregation
    )
    gammas = parameter_grid(gamma_grid, "gamma", positive=True)
    betas = parameter_grid(beta_grid, "beta")
    etas = parameter_grid(eta_grid, "eta", maximum=1)
    if 0 not in etas or 0 not in betas:
        raise PaceError("Include eta=0 and beta=0 in the grid for nested ablations")
    if len(gammas) * len(betas) * len(etas) > 2000:
        raise PaceError("Use at most 2000 grid combinations")

    def aligned(rows, field="pace_score"):
        mapping = {(r["element_id"], r["gene_id"]): r[field] for r in rows}
        values = [mapping[r["element_id"], r["gene_id"]] for r in labels]
        if not all(math.isfinite(v) for v in values):
            raise PaceError("Weak fitting requires complete scores for every assessed candidate")
        return values

    predictions = {}
    reference_rows = None
    for gamma, beta in product(gammas, betas):
        source = {**prior, "gamma": gamma, "beta": beta}
        rows = compute(cfg, contact_prior_override=source)["scores"]
        reference_rows = rows
        for eta in etas:
            rescored, _ = score(rows, eta=eta, partial_policy=cfg["scoring"]["partial_policy"])
            predictions[gamma, beta, eta] = aligned(rescored)
    baseline_cfg = copy.deepcopy(cfg)
    baseline_cfg["contact"].update(mode="prior_only", reliability=None, reliability_source=None)
    baseline = compute(
        baseline_cfg,
        contact_prior_override={
            **prior,
            "gamma": number(abc_gamma, "abc_gamma", minimum=0),
            "beta": 0,
        },
    )["scores"]
    min_distance = defaultdict(lambda: math.inf)
    for row in reference_rows:
        min_distance[row["element_id"]] = min(min_distance[row["element_id"]], row["distance_bp"])
    nearest = [
        {**r, "nearest": float(r["distance_bp"] == min_distance[r["element_id"]])}
        for r in reference_rows
    ]
    report = fit_chromosome_grid(
        labels,
        predictions,
        test_chromosomes=test_chromosomes,
        baselines={
            "ABC_powerlaw_matched_catalog": aligned(baseline),
            "nearest_TSS": aligned(nearest, "nearest"),
        },
    )
    model = {
        **report,
        "schema_version": "pace-weak-1",
        "kind": "eqtl_weak_calibrator",
        "scope": scope,
        "status": "weak_fitted",
        "aggregation": aggregation,
        "abc_gamma": abc_gamma,
        "is_synthetic": cfg["execution_profile"] == "demonstration",
        "label_assumption": "zero PIP mass is unlabelled background; variants/signals are treated as independent only under the selected approximation",
        "test_scope": "eQTL labels held out by chromosome; source contact prior is fixed",
        "training_label_hash": digest([r for r in labels if r["chrom"] not in test_chromosomes]),
        "test_label_hash": digest([r for r in labels if r["chrom"] in test_chromosomes]),
    }
    selected = predictions[tuple(model["selected"][k] for k in ("gamma", "beta", "eta"))]
    return model, [
        {**r, "prediction": p, "split": "test" if r["chrom"] in test_chromosomes else "train"}
        for r, p in zip(labels, selected, strict=True)
    ]
