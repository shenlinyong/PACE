"""Canonical three-regime pipeline sharing one bulk-proxy mathematical kernel."""

from __future__ import annotations

import math
from collections import defaultdict

import yaml

from . import SCHEMA_VERSION, __version__
from .config import load_config
from .core import activity, score, tss_contact
from .errors import PaceError
from .evidence.assets import capabilities, load_asset
from .evidence.resolve import resolve_activity, resolve_contacts
from .io.tables import write_table
from .provenance import digest, environment, file_hash, output_directory, software_hash, write_json
from .reporting import evidence_catalog, evidence_summary, multiomics_features
from .schemas import load_tables, universe_ids


def compute(cfg: dict):
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
    for kind, path in (
        ("sequence", cfg["sequence"]["model_path"]),
        ("contact_prior", cfg["contact"]["prior_path"]),
        ("fusion", cfg["fusion"]["calibrator_path"]),
    ):
        if path:
            assets[kind] = load_asset(path, cfg, kind=kind)
    if cfg["regime"] == "genome_only" and "sequence" not in assets:
        raise PaceError("genome_only requires a quantitative sequence asset")
    predictions, windows, variants = tables["predictions"], [], None
    if assets.get("sequence") and cfg["genome"]["reference_path"]:
        from .sequence.genome import prepare_windows
        from .sequence.model import predict_windows

        if predictions:
            raise PaceError("Provide either sequence inputs or precomputed predictions, not both")
        windows, variants = prepare_windows(
            cfg, tables["units"], input_length=assets["sequence"]["input_length"]
        )
        predictions = predict_windows(
            windows, assets["sequence"], max_n_fraction=cfg["sequence"]["max_n_fraction"]
        )
    elif assets.get("sequence") and not predictions and not tables["resolved_activity"]:
        raise PaceError("Sequence asset needs reference input or manifest-matched predictions")
    if cfg["regime"] == "genome_only" and not (predictions or tables["resolved_activity"]):
        raise PaceError("No genome-only quantitative evidence available")
    resolved_a = resolve_activity(
        tables,
        cfg,
        predictions=predictions,
        sequence_asset=assets.get("sequence"),
        fusion_asset=assets.get("fusion"),
    )
    resolved_c = resolve_contacts(
        tables, cfg, prior_asset=assets.get("contact_prior"), variants=variants
    )
    a_by_e, c_by_ep, gene_promoters = defaultdict(list), {}, defaultdict(list)
    for r in resolved_a:
        a_by_e[r["element_id"]].append(r)
    for r in resolved_c:
        c_by_ep[r["element_id"], r["promoter_id"]] = r
    for p in tables["promoters"]:
        gene_promoters[p["gene_id"]].append(p)
    units = {r["element_id"]: r for r in tables["units"]}
    edges = []
    for edge in tables["candidates"]:
        element, gene = edge["element_id"], edge["gene_id"]
        ar = a_by_e[element]
        ps = gene_promoters[gene]
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
                "A_used": activity([r["resolved_value"] for r in ar]),
                "Cbar": tss_contact([r["resolved_value"] for r in cs], [p["pi"] for p in ps]),
                "distance_bp": distance,
                "n_tss": len(ps),
                "reason": reason,
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
    scores, summary = score(edges, eta=cfg["allocation"]["eta"])
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
        "resolved_evidence_summary": evidence_summary(resolved_a, resolved_c),
        "multiomics_roles": roles,
        "n_candidates": len(scores),
        "n_scoreable": sum(r["scoreable"] for r in scores),
        "genome_windows": [{k: v for k, v in w.items() if k != "sequences"} for w in windows],
        "biological_validation": "not_assessed",
        "sequence_context": "individual_sequence_prediction"
        if cfg["genome"]["variant_path"]
        else "reference_context_prediction"
        if windows
        else "external_quantitative_evidence"
        if predictions
        else "not_used",
        "capabilities": capabilities(cfg),
        "interpretation": "PACE is a composition of regulatory support, not a causal probability or expression effect.",
    }
    ids = universe_ids(tables, cfg)
    scale_contract = sorted(
        {
            (r["assay"], r["unit"], r.get("normalization_id"), r["window_id"])
            for r in resolved_a
            if r["unit"] is not None
        }
    )
    contract = {
        **ids,
        "estimand": cfg["estimand"],
        "target_level": cfg["target_level"],
        "context": cfg["context"],
        "panel": sorted(cfg["activity"]["panel"]),
        "scales": scale_contract,
        "contact_scale": cfg["contact"]["scale"],
        "eta": cfg["allocation"]["eta"],
        "catalog_profile": cfg["catalog"]["profile"],
        "contact_near_diagonal_bp": cfg["contact"]["near_diagonal_bp"],
        "contact_near_diagonal_policy": cfg["contact"]["near_diagonal_policy"],
        "evidence_policy": {
            "regime": cfg["regime"],
            "contact_mode": cfg["contact"]["mode"],
            "contact_reliability": cfg["contact"]["reliability"],
            "contact_reliability_source": cfg["contact"]["reliability_source"],
            "allow_prior_fallback": cfg["contact"]["allow_prior_fallback"],
            "minimum_callable_fraction": cfg["activity"]["minimum_callable_fraction"],
            "quality_stratum": cfg["fusion"]["quality_stratum"],
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
        "genome_hashes": {
            name: file_hash(path)
            for name, path in cfg["genome"].items()
            if name.endswith("_path") and path
        },
        "asset_hashes": {k: a["manifest_sha256"] for k, a in assets.items()},
        "universe_ids": ids,
        "config_hash": digest(cfg),
    }
    evidence, sources = evidence_catalog(tables, resolved_a, resolved_c, assets, cfg)
    return {
        "scores": scores,
        "gene_summary": summary,
        "resolved_activity": resolved_a,
        "resolved_contacts": resolved_c,
        "features": features,
        "evidence": evidence,
        "sources": sources,
        "qc": qc,
        "manifest": manifest,
    }


def run(config_path, out):
    cfg = load_config(config_path)
    # Validation/computation occurs before publication; failed runs leave no success directory.
    result = compute(cfg)
    with output_directory(out) as dest:
        write_table(dest / "scores.tsv.gz", result["scores"])
        for name in (
            "gene_summary",
            "resolved_activity",
            "resolved_contacts",
            "evidence",
            "sources",
        ):
            write_table(dest / f"{name}.tsv", result[name])
        write_table(dest / "multiomics_features.tsv.gz", result["features"], fields=None)
        write_json(dest / "qc_report.json", result["qc"])
        write_json(dest / "run_manifest.json", result["manifest"])
        (dest / "resolved_config.yaml").write_text(
            yaml.safe_dump(cfg, sort_keys=True), encoding="utf-8"
        )
        (dest / "report.md").write_text(
            f"# PACE run: {cfg['run_id']}\n\n"
            f"Regime: {cfg['regime']}; profile: {cfg['execution_profile']}; estimand: bulk_proxy.\n\n"
            f"Scorable candidates: {result['qc']['n_scoreable']} / {len(result['scores'])}. "
            "Partial normalization is conditional on the measurable subset.\n\n"
            "Scores are relative support shares. Changes in shares alone do not establish changes "
            "in enhancer activity or gene expression. Biological validation is not supplied by these software checks.\n",
            encoding="utf-8",
        )
    return result
