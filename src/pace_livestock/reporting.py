"""Auditable nine-omics annotations, independent of the core score."""

import math
from collections import Counter, defaultdict

from .io.methylation import summarize_methylation
from .provenance import digest

LAYERS = (
    "RNA",
    "accessibility",
    "H3K27ac",
    "H3K4me1",
    "H3K4me3",
    "H3K27me3",
    "CTCF",
    "DNA_methylation",
    "Hi-C",
)


def multiomics_features(tables, scores, resolved_activity, cfg):
    rows = []
    for row in tables["features"]:
        rows.append({**row, "role": "annotation_only"})
    panel = set(cfg["activity"]["panel"])
    for row in tables["observed_activity"]:
        if row["assay"] not in panel:
            rows.append(
                {
                    "entity_type": "element",
                    "entity_id": row["element_id"],
                    "feature_name": row["assay"],
                    "value": row["signal"] if row["measurement_status"] == "observed" else math.nan,
                    "evidence_id": f"observation:{row['sample_id']}:{row['element_id']}:{row['assay']}",
                    "status": row["measurement_status"],
                    "role": "annotation_only",
                }
            )
    if tables["methylation"]:
        for row in summarize_methylation(tables["methylation"], tables["units"]):
            for field in (
                "M_site",
                "M_pooled",
                "covered_cpg",
                "total_reads",
                "cpg_coverage_fraction",
            ):
                rows.append(
                    {
                        "entity_type": "element",
                        "entity_id": row["element_id"],
                        "feature_name": f"DNA_methylation:{field}",
                        "value": row[field],
                        "evidence_id": f"methylation:{row['sample_id']}:{row['element_id']}",
                        "status": row["status"],
                        "role": "annotation_only",
                    }
                )
    for row in tables["expression"]:
        rows.append(
            {
                "entity_type": "gene",
                "entity_id": row["gene_id"],
                "feature_name": "RNA:TPM",
                "value": row["tpm"],
                "evidence_id": f"expression:{row['sample_id']}:{row['gene_id']}",
                "status": row["status"],
                "role": "annotation_only",
            }
        )
    for row in scores:
        for feature, value in (
            ("log1p_A", math.log1p(row["A_used"])),
            ("log1p_C", math.log1p(row["Cbar"])),
            ("log1p_distance", math.log1p(row["distance_bp"])),
            ("pace_score", row["pace_score"]),
        ):
            rows.append(
                {
                    "entity_type": "edge",
                    "entity_id": f"{row['element_id']}|{row['gene_id']}",
                    "element_id": row["element_id"],
                    "gene_id": row["gene_id"],
                    "feature_name": feature,
                    "value": value,
                    "evidence_id": "core_formula",
                    "status": "resolved" if math.isfinite(value) else "unresolved",
                    "role": "active_core",
                }
            )
    # Promoter annotations are aggregated with original pi; missing pi mass is explicit.
    pvalues = {
        (r["entity_id"], r["feature_name"]): r for r in rows if r["entity_type"] == "promoter"
    }
    features = {feature for _, feature in pvalues}
    genes = defaultdict(list)
    for p in tables["promoters"]:
        genes[p["gene_id"]].append(p)
    for gene, promoters in genes.items():
        for feature in sorted(features):
            missing, value = 0.0, 0.0
            for p in promoters:
                f = pvalues.get((p["promoter_id"], feature))
                if p["pi"] > 0 and (not f or not math.isfinite(f["value"])):
                    missing += p["pi"]
                elif p["pi"] > 0:
                    value += p["pi"] * f["value"]
            rows.append(
                {
                    "entity_type": "gene",
                    "entity_id": gene,
                    "feature_name": f"promoter:{feature}",
                    "value": value if missing == 0 else math.nan,
                    "evidence_id": "pi_weighted_annotation",
                    "status": "resolved" if missing == 0 else "unresolved",
                    "missing_pi_mass": missing,
                    "role": "annotation_only",
                }
            )
    roles = {layer: "disabled" for layer in LAYERS}
    roles["Hi-C"] = "active_core"
    if panel & {"ATAC", "DNase"}:
        roles["accessibility"] = "active_core"
    if "H3K27ac" in panel:
        roles["H3K27ac"] = "active_core"
    for row in rows:
        for layer in LAYERS:
            if row["feature_name"].startswith(layer) and roles[layer] == "disabled":
                roles[layer] = "annotation_only"
    return rows, roles


def evidence_summary(activity, contacts):
    return {
        "activity": dict(
            Counter(
                (r["evidence_type"] if r["resolution_status"] == "resolved" else "unresolved")
                for r in activity
            )
        ),
        "contact": dict(
            Counter(
                (r["evidence_type"] if r["resolution_status"] == "resolved" else "unresolved")
                for r in contacts
            )
        ),
    }


def evidence_catalog(tables, resolved_activity, resolved_contacts, assets, cfg):
    """Publish resolvable parent IDs for original observations and derived estimates."""
    sources = {r["source_id"]: dict(r) for r in tables["sources"]}
    evidence = {r["evidence_id"]: dict(r) for r in tables["evidence"]}
    samples = {r["sample_id"]: r for r in tables["samples"]}
    for r in tables["observed_activity"]:
        key = f"observation:{r['element_id']}:{r['sample_id']}:{r['assay']}"
        evidence[key] = {
            "evidence_id": key,
            "evidence_type": "observed",
            "source_id": samples[r["sample_id"]]["source_id"],
            "parent_evidence_ids": None,
            "model_id": None,
            "unit": r["unit"],
            "processing_method": "input_normalized_observation",
            "checksum": digest(r),
        }
    for r in tables["observed_contacts"]:
        key = f"contact:{r['element_id']}:{r['promoter_id']}:{r['sample_id']}"
        evidence[key] = {
            "evidence_id": key,
            "evidence_type": "observed",
            "source_id": r["source_id"],
            "parent_evidence_ids": None,
            "model_id": None,
            "unit": r["scale"],
            "processing_method": "input_contact_observation",
            "checksum": digest(r),
        }
    for kind, asset in assets.items():
        key = f"model:{asset['model_id']}"
        sources[key] = {
            "source_id": key,
            "path_or_accession": asset["model_id"],
            "source_type": kind,
            "assembly": cfg["context"]["assembly"],
            "processing_method": "manifest_checked_model",
            "normalization_id": asset.get("normalization_id"),
            "checksum": asset["manifest_sha256"],
        }
    generated_source = "pace:derived_evidence"
    sources[generated_source] = {
        "source_id": generated_source,
        "path_or_accession": cfg["run_id"],
        "source_type": "derived",
        "assembly": cfg["context"]["assembly"],
        "processing_method": "PACE_bulk_proxy_resolution",
        "normalization_id": None,
        "checksum": digest(cfg),
    }
    for r in resolved_activity + resolved_contacts:
        if r["evidence_id"] in evidence:
            continue
        model = r.get("model_id") or r.get("prior_id")
        parents = r.get("parent_evidence_ids")
        # Derived rows may combine observations and a model; model identity remains a separate field.
        evidence[r["evidence_id"]] = {
            "evidence_id": r["evidence_id"],
            "evidence_type": r["evidence_type"],
            "source_id": f"model:{model}" if model else generated_source,
            "parent_evidence_ids": parents,
            "model_id": model,
            "unit": r.get("unit", r.get("scale")),
            "processing_method": r.get("aggregation_method", r.get("resolved_mode", "imported")),
            "checksum": digest(r),
        }
    return sorted(evidence.values(), key=lambda r: r["evidence_id"]), sorted(
        sources.values(), key=lambda r: r["source_id"]
    )
