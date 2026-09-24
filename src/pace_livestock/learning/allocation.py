"""Bounded allocation calibration using within-gene functional ranking contrasts.

The convex objective is a ranking surrogate, never a Bernoulli likelihood of the
compositional PACE score. I/O and applicability checks surround a pure optimizer.
"""

from __future__ import annotations

import math
from collections import Counter, defaultdict

import numpy as np

from ..core import score
from ..errors import PaceError
from ..evaluation.metrics import binary_metrics
from ..io.tables import integer, number, read_table
from ..provenance import digest, file_hash, read_json

LABEL_FIELDS = (
    "label_id element_id gene_id species assembly context_id target_level "
    "perturbation_type effect_direction label_status split group_id assay_id source_id"
).split()
PERTURBATIONS = {"CRISPRi", "deletion", "enhancer_inhibition", "synthetic_inhibition"}
METHOD = "gene_balanced_pairwise_logistic_v1"
VALIDATED_METHOD = "grouped_ap_one_se_pairwise_v2"
SHRINKAGE_GRID = tuple(i / 20 for i in range(21))


def fit_eta(groups: dict, *, minimum_genes: int = 3) -> dict:
    """Minimize mean_gene mean_pos,neg softplus(-(delta_logAC + eta*delta_logB)).

    Each group contains two finite N x 2 arrays (positives, negatives). Bisection
    uses the monotone derivative of this convex, one-dimensional objective.
    Pair blocks bound temporary memory independently of the number of pairs.
    """
    minimum_genes = integer(minimum_genes, "minimum_genes", minimum=1)
    informative = {}
    for gene, (positive, negative) in sorted(groups.items()):
        pos, neg = np.asarray(positive, dtype=float), np.asarray(negative, dtype=float)
        if any(a.ndim != 2 or a.shape[1] != 2 or not np.isfinite(a).all() for a in (pos, neg)):
            raise PaceError("Eta fitting requires finite N x 2 log-evidence arrays")
        if not len(pos) or not len(neg):
            continue
        # An allocation contrast below numerical resolution carries no fitting information.
        if (
            max(abs(pos[:, 1].max() - neg[:, 1].min()), abs(pos[:, 1].min() - neg[:, 1].max()))
            > 1e-12
        ):
            informative[gene] = (pos, neg)
    report = {
        "eta": 0.0,
        "status": "fallback",
        "reason": "insufficient_informative_genes",
        "method": METHOD,
        "minimum_genes": minimum_genes,
        "n_informative_genes": len(informative),
        "n_pairs": sum(len(p) * len(n) for p, n in informative.values()),
        "fit_genes": sorted(informative),
    }
    if len(informative) < minimum_genes:
        return report

    def evaluate(eta):
        loss = derivative = 0.0
        for pos, neg in informative.values():
            weight = 1.0 / (len(informative) * len(pos) * len(neg))
            for i in range(0, len(pos), 256):
                for j in range(0, len(neg), 256):
                    delta = pos[i : i + 256, None, :] - neg[None, j : j + 256, :]
                    margin = delta[..., 0] + eta * delta[..., 1]
                    loss += float(np.logaddexp(0, -margin).sum()) * weight
                    derivative -= (
                        float((delta[..., 1] * np.exp(-np.logaddexp(0, margin))).sum()) * weight
                    )
        return loss, derivative

    baseline, d0 = evaluate(0.0)
    _, d1 = evaluate(1.0)
    if d0 >= 0:
        eta = 0.0
    elif d1 <= 0:
        eta = 1.0
    else:
        low, high = 0.0, 1.0
        for _ in range(40):
            midpoint = (low + high) / 2
            if evaluate(midpoint)[1] < 0:
                low = midpoint
            else:
                high = midpoint
        eta = (low + high) / 2
    loss, derivative = evaluate(eta)
    return {
        **report,
        "eta": eta,
        "status": "fitted",
        "reason": "bounded_ranking_optimum",
        "loss_at_zero": baseline,
        "loss_at_eta": loss,
        "derivative_at_eta": derivative,
        "boundary": eta in (0.0, 1.0),
        "validation": "not_independently_assessed",
    }


def calibration_scope(cfg: dict, tables: dict, assets: dict, scales: list) -> dict:
    from ..evidence.resolve import contact_measurement_contract
    from ..schemas import universe_ids

    return {
        "context": cfg["context"],
        "target_level": cfg["target_level"],
        "estimand": cfg["estimand"],
        "regime": cfg["regime"],
        # File locations are not scientific identities. The catalog universe
        # below binds chromosome-size contents and canonical coordinates.
        "catalog": {k: v for k, v in cfg["catalog"].items() if not k.endswith("_path")},
        "panel": sorted(cfg["activity"]["panel"]),
        "activity_policy": cfg["activity"],
        "promoter_selection": cfg["promoters"],
        "scales": scales,
        "contact_policy": {k: v for k, v in cfg["contact"].items() if not k.endswith("_path")},
        "contact_measurement_contract": contact_measurement_contract(
            tables, cfg, prior_asset=assets.get("contact_prior")
        ),
        "promoter_weights": sorted(
            (p["gene_id"], p["promoter_id"], p["pi"]) for p in tables["promoters"]
        ),
        "universe_ids": universe_ids(tables, cfg),
        "asset_hashes": {k: a["manifest_sha256"] for k, a in assets.items()},
        "allocation_validation_policy": _validation_policy(cfg),
    }


def _eligible_labels(rows: list[dict], cfg: dict) -> tuple[list[dict], list[dict]]:
    usable, excluded, seen = [], [], set()
    splits = {k: {} for k in ("gene_id", "element_id", "group_id")}
    for row in rows:
        if any(row.get(k) is None for k in LABEL_FIELDS):
            raise PaceError("Eta label fields must not be missing; use explicit status values")
        if row["label_id"] in seen:
            raise PaceError(f"Duplicate eta label_id: {row['label_id']}")
        seen.add(row["label_id"])
        if row["split"] not in ("train", "calibration", "test"):
            raise PaceError("Eta label split must be train, calibration or test")
        for key, seen_split in splits.items():
            previous = seen_split.setdefault(row[key], row["split"])
            if previous != row["split"]:
                raise PaceError(f"Eta labels leak {key}={row[key]} across splits")
        reason = None
        if any(row[k] != cfg["context"][k] for k in ("species", "assembly", "context_id")):
            reason = "context_mismatch"
        elif row["target_level"] != cfg["target_level"]:
            reason = "target_level_mismatch"
        elif row["split"] == "test":
            reason = "held_out_test"
        elif row["perturbation_type"] not in PERTURBATIONS:
            reason = "unsupported_perturbation"
        elif (
            row["perturbation_type"] == "synthetic_inhibition"
            and cfg["execution_profile"] != "demonstration"
        ):
            raise PaceError("Synthetic eta labels require execution_profile=demonstration")
        elif row["label_status"] not in ("enhancing_positive", "powered_negative"):
            reason = "unusable_label_status"
        elif row["label_status"] == "enhancing_positive" and row["effect_direction"] != "down":
            reason = "positive_requires_downregulation"
        elif row["label_status"] == "powered_negative" and row["effect_direction"] != "none":
            reason = "negative_requires_no_effect"
        if reason:
            excluded.append({"label_id": row["label_id"], "reason": reason})
        else:
            usable.append(row)
    pairs = [(r["element_id"], r["gene_id"]) for r in usable]
    if len(pairs) != len(set(pairs)):
        raise PaceError("Repeated eta labels for one edge need predeclared aggregation")
    return usable, excluded


def _validation_policy(cfg: dict) -> dict:
    allocation = cfg["allocation"]
    groups = integer(allocation.get("minimum_groups", 3), "minimum_groups", minimum=3)
    folds = integer(allocation.get("validation_folds", 5), "validation_folds", minimum=3)
    positive = number(
        allocation.get("minimum_positive_fraction", 0.8),
        "minimum_positive_fraction",
        minimum=0,
        maximum=1,
    )
    if positive <= 0.5:
        raise PaceError("minimum_positive_fraction must be greater than 0.5")
    return {
        "method": VALIDATED_METHOD,
        "minimum_genes": allocation["minimum_genes"],
        "minimum_groups": groups,
        "validation_folds": folds,
        "minimum_positive_fraction": positive,
        "shrinkage_grid": list(SHRINKAGE_GRID),
        "metric": "independent_group_mean_of_within_gene_average_precision",
        "selection": "smallest_shrinkage_within_one_se_of_best",
        "gain_gate": "paired_group_mean_gain_greater_than_one_se",
    }


def _gene_pairs(records: list[dict]) -> dict:
    pairs = defaultdict(lambda: ([], []))
    for row in records:
        index = 0 if row["label_status"] == "enhancing_positive" else 1
        pairs[row["gene_id"]][index].append(row["log_evidence"])
    return {gene: pair for gene, pair in pairs.items() if pair[0] and pair[1]}


def _independent_components(records: list[dict]) -> list[list[dict]]:
    """Keep overlapping group, gene and element identities in the same fold.

    The submitted group_id must represent the broadest known dependency (for
    example locus, chromosome block or experiment). Connected components prevent
    finer inconsistent IDs from separating a shared gene or physical element.
    """
    parent = list(range(len(records)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    seen = {}
    for i, row in enumerate(records):
        for field in ("group_id", "gene_id", "element_id"):
            key = (field, row[field])
            if key in seen:
                parent[find(i)] = find(seen[key])
            else:
                seen[key] = i
    components = defaultdict(list)
    for i, row in enumerate(records):
        components[find(i)].append(row)
    return sorted(components.values(), key=lambda rows: min(r["label_id"] for r in rows))


def _component_ap(records: list[dict], eta: float) -> float:
    # Within one gene the denominator is constant, so finite log supports give
    # exactly the PACE ordering without comparing unrelated compositional scores.
    per_gene = []
    for positive, negative in _gene_pairs(records).values():
        evidence = np.asarray(positive + negative)
        labels = [1] * len(positive) + [0] * len(negative)
        per_gene.append(
            binary_metrics(labels, evidence[:, 0] + eta * evidence[:, 1])["average_precision"]
        )
    return float(np.mean(per_gene)) if per_gene else math.nan


def _se(values) -> float:
    values = np.asarray(values, dtype=float)
    return float(values.std(ddof=1) / math.sqrt(len(values))) if len(values) > 1 else math.inf


def _gain_summary(baseline, selected, policy) -> dict:
    gain = np.asarray(selected) - np.asarray(baseline)
    mean, uncertainty = float(gain.mean()), _se(gain)
    positive = float(np.mean(gain > 1e-12))
    return {
        "mean_ap_at_zero": float(np.mean(baseline)),
        "mean_ap_selected": float(np.mean(selected)),
        "mean_paired_ap_gain": mean,
        "paired_gain_se": uncertainty,
        "positive_group_fraction": positive,
        "passed": bool(
            mean > uncertainty + 1e-12 and positive >= policy["minimum_positive_fraction"]
        ),
    }


def _validated_fit(
    records: list[dict], confirmation: list[dict], policy: dict, *, require_confirmation=False
) -> dict:
    fitted = fit_eta(_gene_pairs(records), minimum_genes=policy["minimum_genes"])
    result = {
        **{
            k: v
            for k, v in fitted.items()
            if k not in ("loss_at_eta", "derivative_at_eta", "boundary")
        },
        "candidate_fit": fitted,
        "method": VALIDATED_METHOD,
        "candidate_method": METHOD,
        "candidate_eta": fitted["eta"],
        "eta": 0.0,
        "status": "fallback",
        "validation": "not_passed",
        "deployment_validated": False,
        "validation_policy": policy,
    }
    if fitted["status"] != "fitted":
        return result
    components = _independent_components(records)
    result["n_independent_groups"] = len(components)
    if len(components) < policy["minimum_groups"]:
        return {**result, "reason": "insufficient_independent_groups"}
    if fitted["eta"] <= 0:
        return {**result, "reason": "candidate_prefers_zero"}
    n_folds = min(policy["validation_folds"], len(components))
    folds = []
    ap = np.empty((len(components), len(SHRINKAGE_GRID)))
    for fold in range(n_folds):
        held_indices = [i for i in range(len(components)) if i % n_folds == fold]
        train = [r for i, c in enumerate(components) if i % n_folds != fold for r in c]
        candidate = fit_eta(_gene_pairs(train), minimum_genes=policy["minimum_genes"])
        if candidate["status"] != "fitted":
            return {**result, "reason": "insufficient_training_genes_in_cv"}
        folds.append(
            {
                "fold": fold,
                "candidate_eta": candidate["eta"],
                "training_group_ids": sorted({r["group_id"] for r in train}),
                "validation_group_ids": sorted(
                    {r["group_id"] for i in held_indices for r in components[i]}
                ),
            }
        )
        for i in held_indices:
            for j, shrinkage in enumerate(SHRINKAGE_GRID):
                ap[i, j] = _component_ap(components[i], candidate["eta"] * shrinkage)
    if not np.isfinite(ap).all():
        return {**result, "reason": "validation_group_without_testable_gene"}
    means = ap.mean(axis=0)
    best = int(np.argmax(means))
    threshold = float(means[best]) - _se(ap[:, best])
    selected = next(i for i, mean in enumerate(means) if mean >= threshold - 1e-12)
    summary = _gain_summary(ap[:, 0], ap[:, selected], policy)
    cv = {
        "n_folds": n_folds,
        "folds": folds,
        "mean_ap_by_shrinkage": [float(v) for v in means],
        "ap_by_independent_group": ap.tolist(),
        "independent_groups": [
            {
                field: sorted({row[field] for row in component})
                for field in ("group_id", "gene_id", "element_id")
            }
            for component in components
        ],
        "one_se_threshold": threshold,
        "selected_shrinkage": SHRINKAGE_GRID[selected],
        **summary,
    }
    result["cross_validation"] = cv
    if selected == 0 or not summary["passed"]:
        return {**result, "reason": "no_stable_independent_ap_gain"}
    eta = fitted["eta"] * SHRINKAGE_GRID[selected]
    if require_confirmation:
        held = _independent_components(confirmation)
        if len(held) < policy["minimum_groups"]:
            return {**result, "reason": "insufficient_calibration_groups"}
        baseline = [_component_ap(c, 0) for c in held]
        predicted = [_component_ap(c, eta) for c in held]
        if not all(math.isfinite(v) for v in baseline + predicted):
            return {**result, "reason": "calibration_group_without_testable_gene"}
        confirmation_report = _gain_summary(baseline, predicted, policy)
        result["confirmation"] = {"n_independent_groups": len(held), **confirmation_report}
        if not confirmation_report["passed"]:
            return {**result, "reason": "calibration_gain_not_confirmed"}
    return {
        **result,
        "eta": eta,
        "status": "fitted",
        "reason": "independent_ap_gain_validated",
        "validation": "independent_grouped_cv",
        "deployment_validated": True,
    }


def fit_from_labels(edges: list[dict], rows: list[dict], cfg: dict) -> dict:
    usable, excluded = _eligible_labels(rows, cfg)
    at_zero = {(r["element_id"], r["gene_id"]): r for r in score(edges, eta=0)[0]}
    at_one = {(r["element_id"], r["gene_id"]): r for r in score(edges, eta=1)[0]}
    used = []
    for row in usable:
        key = row["element_id"], row["gene_id"]
        base, full = at_zero.get(key), at_one.get(key)
        if base is None or full is None:
            reason = "outside_candidate_universe"
        elif not all(math.isfinite(r["log_support"]) for r in (base, full)):
            reason = "zero_or_unresolved_support"
        else:
            reason = None
        if reason:
            excluded.append({"label_id": row["label_id"], "reason": reason})
        else:
            used.append(
                {
                    **row,
                    "log_evidence": [
                        base["log_support"],
                        full["log_support"] - base["log_support"],
                    ],
                }
            )
    policy = _validation_policy(cfg)
    paired_genes = set(_gene_pairs(used))
    consumed = [r for r in used if r["gene_id"] in paired_genes]
    excluded.extend(
        {"label_id": r["label_id"], "reason": "no_within_gene_pair"}
        for r in used
        if r["gene_id"] not in paired_genes
    )
    has_train = any(r["split"] == "train" for r in usable)
    estimation = [r for r in consumed if r["split"] == "train"] if has_train else consumed
    confirmation = [r for r in consumed if r["split"] == "calibration"] if has_train else []
    fitted = _validated_fit(
        estimation,
        confirmation,
        policy,
        require_confirmation=has_train and any(r["split"] == "calibration" for r in usable),
    )
    # The benchmark exclusion contract includes all labels used to estimate,
    # select or approve the parameter, even when selection returns eta=0.
    return {
        **fitted,
        "n_input_labels": len(rows),
        "n_fit_labels": len(consumed),
        "fit_genes": sorted({r["gene_id"] for r in consumed}),
        "fit_label_ids": sorted(r["label_id"] for r in consumed),
        "fit_element_ids": sorted({r["element_id"] for r in consumed}),
        "fit_group_ids": sorted({r["group_id"] for r in consumed}),
        "fit_splits": sorted({r["split"] for r in consumed}),
        "estimation_label_ids": sorted(r["label_id"] for r in estimation),
        "confirmation_label_ids": sorted(r["label_id"] for r in confirmation),
        "estimation_split": "train" if has_train else "calibration_cross_validation",
        "excluded_counts": dict(Counter(r["reason"] for r in excluded)),
        "excluded_labels": excluded,
    }


def resolve_eta(edges: list[dict], cfg: dict, scope: dict) -> dict:
    allocation = cfg["allocation"]
    requested = allocation["eta"]
    base = {
        "schema_version": "pace-eta-2",
        "kind": "allocation_calibrator",
        "scope": scope,
        "scope_sha256": digest(scope),
        "is_synthetic": cfg["execution_profile"] == "demonstration",
        "requested": requested,
    }
    if requested != "auto":
        if allocation["labels_path"] or allocation["calibrator_path"]:
            raise PaceError("Manual eta cannot also request calibration labels or an artifact")
        return {**base, "eta": number(requested, "eta", minimum=0, maximum=1), "status": "fixed"}
    if allocation["labels_path"] and allocation["calibrator_path"]:
        raise PaceError("Provide eta labels OR an eta calibrator, not both")
    if allocation["calibrator_path"]:
        path = allocation["calibrator_path"]
        fitted = read_json(path)
        if (
            not isinstance(fitted, dict)
            or fitted.get("schema_version") != "pace-eta-2"
            or fitted.get("kind") != "allocation_calibrator"
        ):
            raise PaceError("Invalid eta calibrator schema/kind")
        if fitted.get("status") not in ("fixed", "fallback", "fitted"):
            raise PaceError("Eta calibrator must record a fixed, fallback or fitted status")
        if fitted["status"] == "fitted" and fitted.get("method") != VALIDATED_METHOD:
            raise PaceError("Unsupported eta calibration objective")
        if fitted.get("scope_sha256") != digest(fitted.get("scope")) or fitted.get(
            "scope_sha256"
        ) != digest(scope):
            raise PaceError("Eta calibrator scope differs from this run's scientific contract")
        if fitted.get("is_synthetic") is not False and cfg["execution_profile"] != "demonstration":
            raise PaceError("Synthetic or unspecified eta calibration cannot enter a research run")
        eta = number(fitted.get("eta"), "calibrated eta", minimum=0, maximum=1)
        if eta > 0 and (
            fitted["status"] != "fitted"
            or fitted.get("deployment_validated") is not True
            or fitted.get("validation") != "independent_grouped_cv"
            or fitted.get("validation_policy") != _validation_policy(cfg)
        ):
            raise PaceError(
                "Automatic nonzero eta requires a validated grouped calibration artifact"
            )
        return {
            **fitted,
            "eta": eta,
            "requested": requested,
            "reuse": True,
            "calibrator_sha256": file_hash(path),
        }
    if not allocation["labels_path"]:
        return {**base, "eta": 0.0, "status": "fallback", "reason": "no_functional_labels"}
    path = allocation["labels_path"]
    labels = read_table(path, required=LABEL_FIELDS)
    return {
        **base,
        **fit_from_labels(edges, labels, cfg),
        "labels_sha256": file_hash(path),
        "fitting_evidence_sha256": digest(
            sorted((r["element_id"], r["gene_id"], r["A_used"], r["Cbar"]) for r in edges)
        ),
    }
