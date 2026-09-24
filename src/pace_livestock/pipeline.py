"""Measured activity and contact support on a fixed regulatory catalog."""

from __future__ import annotations

import math
from collections import defaultdict

import yaml

from . import SCHEMA_VERSION, __version__
from .config import load_config
from .core import activity, score, tss_contact
from .errors import PaceError
from .evidence.assets import builtin_contact_prior, capabilities, load_asset
from .evidence.promoters import resolve_promoter_weights
from .evidence.resolve import contact_measurement_contract, resolve_activity, resolve_contacts
from .io.tables import write_table
from .provenance import digest, environment, file_hash, output_directory, software_hash, write_json
from .reporting import evidence_catalog, evidence_summary, multiomics_features
from .schemas import load_tables, universe_ids


def compute(cfg: dict, *, contact_prior_override=None):
    if set(cfg) & {"sequence", "fusion", "genome"} or "predictions" in cfg["inputs"]:
        raise PaceError("Unsupported configuration: PACE requires measured activity inputs")
    if cfg["regime"] != "measured":
        raise PaceError("PACE requires measured activity; unsupported evidence regime")
    tables = load_tables(cfg)
    if cfg["regime"] == "measured":
        available = {r["assay"] for r in tables["observed_activity"] + tables["resolved_activity"]}
        if not set(cfg["activity"]["panel"]) <= available:
            raise PaceError(
                "Measured mode requires tables for every declared activity assay; choose an explicit single-layer panel when appropriate"
            )
    if cfg["contact"]["mode"] == "observed" and not (
        tables["observed_contacts"]
        or tables["resolved_contacts"]
        or (cfg["contact"]["allow_prior_fallback"] and cfg["contact"]["prior_path"])
    ):
        raise PaceError(
            "Observed contact mode needs contact observations or an explicitly allowed prior fallback"
        )
    assets = {}
    if cfg["contact"]["prior_path"]:
        assets["contact_prior"] = load_asset(
            cfg["contact"]["prior_path"], cfg, kind="contact_prior"
        )
    elif cfg["contact"].get("prior_preset"):
        if tables["observed_contacts"] or tables["resolved_contacts"]:
            raise PaceError(
                "Human relative prior cannot be mixed with measured contact tables; omit contact inputs for this explicit baseline"
            )
        assets["contact_prior"] = builtin_contact_prior(cfg)
    if contact_prior_override is not None:
        assets["contact_prior"] = contact_prior_override
    resolved_a = resolve_activity(tables, cfg)
    weak_allocation = None
    if cfg["allocation"].get("weak_model_path"):
        from .weak_labels import apply_weak_model

        assets["contact_prior"], weak_allocation = apply_weak_model(cfg, tables, assets, resolved_a)
    resolved_c = resolve_contacts(tables, cfg, prior_asset=assets.get("contact_prior"))
    prior_metadata = assets.get("contact_prior", {})
    for row in resolved_c:
        if row.get("prior_id"):
            row["prior_source_context"] = prior_metadata.get(
                "source_context_id", prior_metadata.get("context_id")
            )
            row["prior_transfer_status"] = prior_metadata.get("transfer_status", "matched_context")
    tables["promoters"], promoter_summary = resolve_promoter_weights(tables, resolved_c, cfg)
    for row in resolved_a:
        row["activity_pseudocount"] = cfg["activity"]["pseudocounts"].get(row["assay"], 0.0)
    a_by_e, c_by_ep, gene_promoters = defaultdict(list), {}, defaultdict(list)
    for r in resolved_a:
        a_by_e[r["element_id"]].append(r)
    for r in resolved_c:
        c_by_ep[r["element_id"], r["promoter_id"]] = r
    for p in tables["promoters"]:
        gene_promoters[p["gene_id"]].append(p)
    units = {r["element_id"]: r for r in tables["units"]}
    bounds = {(r["element_id"], r["gene_id"]): r for r in tables["support_bounds"]}
    edges = []
    for edge in tables["candidates"]:
        element, gene = edge["element_id"], edge["gene_id"]
        ar = a_by_e[element]
        ps = gene_promoters[gene]
        promoter_info = promoter_summary[gene]
        cs = [c_by_ep[element, p["promoter_id"]] for p in ps]
        distance = min(abs(units[element]["anchor0"] - p["tss0"]) for p in ps)
        structural = sorted({r.get("structural_status", "not_assessed") for r in ar + cs})
        reason = next((r["reason"] for r in ar if r["resolution_status"] != "resolved"), "")
        if not reason:
            reason = next(
                (
                    r["reason"]
                    for r, p in zip(cs, ps, strict=True)
                    if p["pi"] > 0 and r["resolution_status"] != "resolved"
                ),
                "",
            )
        edges.append(
            {
                **edge,
                **bounds.get((element, gene), {}),
                "A_used": activity(
                    [r["resolved_value"] for r in ar],
                    [r["activity_pseudocount"] for r in ar],
                ),
                "Cbar": tss_contact([r["resolved_value"] for r in cs], [p["pi"] for p in ps])
                if promoter_info["tss_policy_status"] != "insufficient_retained_weight"
                else math.nan,
                **promoter_info,
                "distance_bp": distance,
                "n_tss": len(ps),
                "n_tss_used": sum(p["pi"] > 0 for p in ps),
                "n_contact_bins": len({(p["chrom"], p["tss0"] // cs[0]["resolution"]) for p in ps})
                if cs[0].get("resolution")
                else len(ps),
                "reason": "insufficient_retained_tss_weight"
                if promoter_info["tss_policy_status"] == "insufficient_retained_weight"
                else reason,
                "activity_sources": ";".join(sorted({r["evidence_type"] for r in ar})),
                "contact_sources": ";".join(sorted({r["evidence_type"] for r in cs})),
                "regime": cfg["regime"],
                "estimand": cfg["estimand"],
                "target_level": cfg["target_level"],
                "structural_status": ";".join(structural),
                "pace_ml_score": math.nan,
                "pace_ml_probability": math.nan,
            }
        )
    from .learning.allocation import calibration_scope, resolve_eta

    scale_contract = sorted(
        {
            (r["assay"], r["unit"], r.get("normalization_id"), r["window_id"])
            for r in resolved_a
            if r["unit"] is not None
        },
        key=lambda row: tuple("" if v is None else str(v) for v in row),
    )
    contact_contract = contact_measurement_contract(
        tables, cfg, prior_asset=assets.get("contact_prior")
    )
    ml_contract = {
        "target_level": cfg["target_level"],
        "estimand": cfg["estimand"],
        "activity_panel": sorted(cfg["activity"]["panel"]),
        "activity_pseudocounts": cfg["activity"]["pseudocounts"],
        "activity_scales": [list(row) for row in scale_contract],
        "contact_definition": {
            **contact_contract,
            "resolution_bp": contact_contract["resolution"],
            "mode": cfg["contact"]["mode"],
            "near_diagonal_bp": cfg["contact"]["near_diagonal_bp"],
            "near_diagonal_policy": cfg["contact"]["near_diagonal_policy"],
            "allow_prior_fallback": cfg["contact"]["allow_prior_fallback"],
            "pseudocount": cfg["contact"]["pseudocount"],
            "pseudocount_distance_bp": cfg["contact"]["pseudocount_distance_bp"],
            "pseudocount_strength": cfg["contact"]["pseudocount_strength"],
        },
        "candidate_construction": {
            "profile": cfg["catalog"]["profile"],
            "window_bp": cfg["catalog"]["width_bp"],
            "offset_bp": cfg["catalog"]["offset_bp"],
            "radius_bp": cfg["catalog"].get("candidate_radius_bp", 5_000_000),
            "include_promoter_units": cfg["catalog"]["include_promoter_units"],
            "promoter_weights": cfg["promoters"]["weights"],
            "promoter_selection": cfg["promoters"],
        },
        "auxiliary_feature_policy": {
            "replicate_aggregation": cfg["activity"]["replicate_aggregation"],
            "assay_definitions": sorted(
                {
                    tuple(
                        str(r.get(k) or "")
                        for k in ("assay", "unit", "normalization_id", "window_id")
                    )
                    for r in tables["observed_activity"]
                    if r["assay"] not in cfg["activity"]["panel"]
                }
            ),
            "feature_definitions": sorted(
                {
                    tuple(
                        str(r.get(k) or "")
                        for k in (
                            "entity_type",
                            "feature_name",
                            "assay",
                            "unit",
                            "normalization_id",
                            "window_id",
                        )
                    )
                    for r in tables["features"]
                }
            ),
            "methylation_assays": sorted({r["assay"] for r in tables["methylation"]}),
            "methylation": {
                k: v for k, v in cfg["methylation"].items() if k != "reference_cpg_path"
            },
            "reference_cpg_sha256": file_hash(cfg["methylation"]["reference_cpg_path"])
            if cfg["methylation"]["reference_cpg_path"]
            else None,
        },
    }
    allocation = weak_allocation or resolve_eta(
        edges, cfg, calibration_scope(cfg, tables, assets, scale_contract)
    )
    ml_contract["allocation_eta"] = allocation["eta"]
    ml_contract["evidence_policy"] = {
        "regime": cfg["regime"],
        "contact_reliability": cfg["contact"]["reliability"],
        "contact_reliability_source": cfg["contact"]["reliability_source"],
        "contact_kappa": cfg["contact"]["kappa"],
        "minimum_callable_fraction": cfg["activity"]["minimum_callable_fraction"],
        "asset_hashes": {k: a["manifest_sha256"] for k, a in assets.items()},
    }
    ml_contract["scoring_policy"] = cfg["scoring"]
    scores, summary = score(
        edges, eta=allocation["eta"], partial_policy=cfg["scoring"]["partial_policy"]
    )
    for row in scores:
        row["allocation_evidence"] = "eqtl_weak" if weak_allocation else "functional_or_fixed"
    for row in summary:
        row.update(promoter_summary[row["gene_id"]])
    features, roles = multiomics_features(tables, scores, resolved_a, cfg)
    if cfg["multiomics"]["mode"] == "ml":
        if not cfg["multiomics"]["model_path"]:
            raise PaceError("multiomics.mode=ml requires a trained classifier")
        from .learning.model import predict_score_rows

        scores = predict_score_rows(
            scores,
            features,
            cfg["multiomics"]["model_path"],
            execution_profile=cfg["execution_profile"],
            context=cfg["context"],
            feature_contract=ml_contract,
        )
        for row in features:
            if row["role"] == "annotation_only":
                row["role"] = (
                    "active_ml"
                    if any(s.get("ml_status") == "resolved" for s in scores)
                    and row["feature_name"] in scores[0].get("ml_features", "").split(";")
                    else "annotation_only"
                )
        for layer in roles:
            if any(
                r["role"] == "active_ml" and r["feature_name"].startswith(layer) for r in features
            ):
                roles[layer] = "active_ml"
    qc = {
        "allocation": {
            k: allocation[k] for k in ("eta", "status", "reason", "reuse") if k in allocation
        },
        "resolved_evidence_summary": evidence_summary(resolved_a, resolved_c),
        "multiomics_roles": roles,
        "n_candidates": len(scores),
        "n_scoreable": sum(r["scoreable"] for r in scores),
        "n_full_scores": sum(
            math.isfinite(r["pace_score"]) and r["score_scope"] == "full_candidate_set"
            for r in scores
        ),
        "n_partial_genes": sum(r["normalization_status"] == "partial" for r in summary),
        "score_bounds": "Sensitivity ranges conditional on support assumptions; not confidence intervals",
        "contact_processing": {
            "prior_source_context": prior_metadata.get(
                "source_context_id", prior_metadata.get("context_id")
            ),
            "prior_transfer_status": prior_metadata.get("transfer_status"),
            "pseudocount_policy": cfg["contact"]["pseudocount"],
            "compatible_prior_available": bool(assets.get("contact_prior")),
            "n_regularized": sum(r["evidence_type"] == "regularized" for r in resolved_c),
            "n_near_diagonal_unresolved": sum(
                r["reason"] == "near_diagonal_unresolved" for r in resolved_c
            ),
            "n_resolved_zero": sum(
                r["resolution_status"] == "resolved" and r["resolved_value"] == 0
                for r in resolved_c
            ),
        },
        "promoter_processing": {
            "settings": cfg["promoters"],
            "n_filtered_genes": sum(
                s["tss_policy_status"] == "filtered" for s in promoter_summary.values()
            ),
            "n_failed_genes": sum(
                s["tss_policy_status"] == "insufficient_retained_weight"
                for s in promoter_summary.values()
            ),
        },
        "activity_pseudocounts": cfg["activity"]["pseudocounts"],
        "biological_validation": "not_assessed",
        "capabilities": capabilities(cfg),
        "interpretation": "PACE is a composition of regulatory support, not a causal probability or expression effect.",
    }
    ids = universe_ids(tables, cfg)
    contract = {
        **ids,
        "estimand": cfg["estimand"],
        "target_level": cfg["target_level"],
        "context": cfg["context"],
        "panel": sorted(cfg["activity"]["panel"]),
        "activity_pseudocounts": cfg["activity"]["pseudocounts"],
        "promoter_selection": cfg["promoters"],
        "scales": scale_contract,
        "contact_scale": cfg["contact"]["scale"],
        "contact_measurement_contract": contact_contract,
        "eta": allocation["eta"],
        "catalog_profile": cfg["catalog"]["profile"],
        "scoring_policy": cfg["scoring"],
        "contact_near_diagonal_bp": cfg["contact"]["near_diagonal_bp"],
        "contact_near_diagonal_policy": cfg["contact"]["near_diagonal_policy"],
        "evidence_policy": {
            "regime": cfg["regime"],
            "contact_mode": cfg["contact"]["mode"],
            "contact_reliability": cfg["contact"]["reliability"],
            "contact_reliability_source": cfg["contact"]["reliability_source"],
            "contact_kappa": cfg["contact"]["kappa"],
            "allow_prior_fallback": cfg["contact"]["allow_prior_fallback"],
            "pseudocount": cfg["contact"]["pseudocount"],
            "pseudocount_distance_bp": cfg["contact"]["pseudocount_distance_bp"],
            "pseudocount_strength": cfg["contact"]["pseudocount_strength"],
            "minimum_callable_fraction": cfg["activity"]["minimum_callable_fraction"],
            "asset_hashes": {k: a["manifest_sha256"] for k, a in assets.items()},
        },
    }
    manifest = {
        "software_version": __version__,
        "software_sha256": software_hash(),
        "schema_version": SCHEMA_VERSION,
        "seed": cfg["seed"],
        "run_id": cfg["run_id"],
        "execution_profile": cfg["execution_profile"],
        "regime": cfg["regime"],
        "environment": environment(),
        "comparison_contract": contract,
        "input_hashes": {name: file_hash(path) for name, path in cfg["inputs"].items() if path},
        "asset_hashes": {k: a["manifest_sha256"] for k, a in assets.items()},
        "universe_ids": ids,
        "asset_manifests": {
            kind: {key: value for key, value in asset.items() if key != "asset_directory"}
            for kind, asset in assets.items()
        },
        "config_hash": digest(cfg),
        "allocation": allocation,
        "ml_feature_contract": ml_contract,
    }
    evidence, sources = evidence_catalog(
        tables, resolved_a, resolved_c, assets, cfg, features=features
    )
    from .reporting import regional_scores

    return {
        "scores": scores,
        "gene_summary": summary,
        "region_scores": regional_scores(scores, tables["region_membership"]),
        "resolved_activity": resolved_a,
        "resolved_contacts": resolved_c,
        "promoter_weights": tables["promoters"],
        "features": features,
        "evidence": evidence,
        "sources": sources,
        "qc": qc,
        "manifest": manifest,
        "eta_calibration": allocation,
    }


def run(config_path, out):
    cfg = config_path if isinstance(config_path, dict) else load_config(config_path)
    # Validation/computation occurs before publication; failed runs leave no success directory.
    result = compute(cfg)
    with output_directory(out) as dest:
        write_table(dest / "scores.tsv.gz", result["scores"])
        for name in (
            "gene_summary",
            "region_scores",
            "resolved_activity",
            "resolved_contacts",
            "promoter_weights",
            "evidence",
            "sources",
        ):
            write_table(dest / f"{name}.tsv", result[name])
        write_table(dest / "multiomics_features.tsv.gz", result["features"], fields=None)
        write_json(dest / "qc_report.json", result["qc"])
        write_json(dest / "run_manifest.json", result["manifest"])
        write_json(dest / "eta_calibration.json", result["eta_calibration"])
        write_json(dest / "ml_feature_contract.json", result["manifest"]["ml_feature_contract"])
        (dest / "resolved_config.yaml").write_text(
            yaml.safe_dump(cfg, sort_keys=True), encoding="utf-8"
        )
        (dest / "report.md").write_text(
            f"# PACE run: {cfg['run_id']}\n\n"
            f"Measured activity; profile: {cfg['execution_profile']}; target: {cfg['target_level']}.\n\n"
            f"Scorable candidates: {result['qc']['n_scoreable']} / {len(result['scores'])}. "
            "Partial normalization is conditional on the measurable subset.\n\n"
            "Scores are relative support shares. Changes in shares alone do not establish changes "
            "in enhancer activity or gene expression. Biological validation is not supplied by these software checks.\n",
            encoding="utf-8",
        )
    return result
