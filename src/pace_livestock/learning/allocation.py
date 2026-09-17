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
from ..io.tables import integer, number, read_table
from ..provenance import digest, file_hash, read_json

LABEL_FIELDS = (
    "label_id element_id gene_id species assembly context_id target_level "
    "perturbation_type effect_direction label_status split group_id assay_id source_id"
).split()
PERTURBATIONS = {"CRISPRi", "deletion", "enhancer_inhibition", "synthetic_inhibition"}
METHOD = "gene_balanced_pairwise_logistic_v1"


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
    from ..schemas import universe_ids

    return {
        "context": cfg["context"],
        "target_level": cfg["target_level"],
        "estimand": cfg["estimand"],
        "regime": cfg["regime"],
        "catalog": cfg["catalog"],
        "panel": sorted(cfg["activity"]["panel"]),
        "activity_policy": cfg["activity"],
        "scales": scales,
        "contact_policy": {k: v for k, v in cfg["contact"].items() if not k.endswith("_path")},
        "promoter_weights": sorted(
            (p["gene_id"], p["promoter_id"], p["pi"]) for p in tables["promoters"]
        ),
        "universe_ids": universe_ids(tables, cfg),
        "asset_hashes": {k: a["manifest_sha256"] for k, a in assets.items()},
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


def fit_from_labels(edges: list[dict], rows: list[dict], cfg: dict) -> dict:
    usable, excluded = _eligible_labels(rows, cfg)
    at_zero = {(r["element_id"], r["gene_id"]): r for r in score(edges, eta=0)[0]}
    at_one = {(r["element_id"], r["gene_id"]): r for r in score(edges, eta=1)[0]}
    groups, used = defaultdict(lambda: ([], [])), []
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
            continue
        index = 0 if row["label_status"] == "enhancing_positive" else 1
        groups[row["gene_id"]][index].append(
            [base["log_support"], full["log_support"] - base["log_support"]]
        )
        used.append(row)
    groups = {g: pair for g, pair in groups.items() if pair[0] and pair[1]}
    fitted = fit_eta(groups, minimum_genes=cfg["allocation"]["minimum_genes"])
    fit_genes = set(fitted["fit_genes"])
    consumed = [r for r in used if r["gene_id"] in fit_genes]
    excluded.extend(
        {"label_id": r["label_id"], "reason": "no_informative_within_gene_pair"}
        for r in used
        if r["gene_id"] not in fit_genes
    )
    return {
        **fitted,
        "n_input_labels": len(rows),
        "n_fit_labels": len(consumed),
        "fit_label_ids": sorted(r["label_id"] for r in consumed),
        "fit_element_ids": sorted({r["element_id"] for r in consumed}),
        "fit_group_ids": sorted({r["group_id"] for r in consumed}),
        "fit_splits": sorted({r["split"] for r in consumed}),
        "excluded_counts": dict(Counter(r["reason"] for r in excluded)),
        "excluded_labels": excluded,
    }


def resolve_eta(edges: list[dict], cfg: dict, scope: dict) -> dict:
    allocation = cfg["allocation"]
    requested = allocation["eta"]
    base = {
        "schema_version": "pace-eta-1",
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
            or fitted.get("schema_version") != "pace-eta-1"
            or fitted.get("kind") != "allocation_calibrator"
        ):
            raise PaceError("Invalid eta calibrator schema/kind")
        if fitted.get("status") not in ("fixed", "fallback", "fitted"):
            raise PaceError("Eta calibrator must record a fixed, fallback or fitted status")
        if fitted["status"] == "fitted" and fitted.get("method") != METHOD:
            raise PaceError("Unsupported eta calibration objective")
        if fitted.get("scope_sha256") != digest(fitted.get("scope")) or fitted.get(
            "scope_sha256"
        ) != digest(scope):
            raise PaceError("Eta calibrator scope differs from this run's scientific contract")
        if fitted.get("is_synthetic") is not False and cfg["execution_profile"] != "demonstration":
            raise PaceError("Synthetic or unspecified eta calibration cannot enter a research run")
        eta = number(fitted.get("eta"), "calibrated eta", minimum=0, maximum=1)
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
