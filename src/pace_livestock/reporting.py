"""Auditable multi-omics annotations, independent of the core score."""

import math
from collections import Counter, defaultdict

from .errors import PaceError
from .evidence.resolve import aggregate_observations
from .io.methylation import load_reference_cpg, promoter_methylation_regions, summarize_methylation
from .provenance import digest

LAYERS = (
    "RNA",
    "accessibility",
    "H3K27ac",
    "H3K4me1",
    "H3K4me3",
    "H3K27me3",
    "H3K9me3",
    "CTCF",
    "DNA_methylation",
    "Hi-C",
)


def observation_evidence_id(row):
    return f"observation:{row['element_id']}:{row['sample_id']}:{row['assay']}"


def expression_evidence_id(row):
    return f"expression:{row['sample_id']}:{row['gene_id']}"


def methylation_evidence_id(sample, assay, entity_type, entity_id):
    return f"methylation:{sample}:{assay}:{entity_type}:{entity_id}"


def methylation_regions(tables, cfg):
    settings = cfg.get("methylation", {})
    return {
        "element": tables["units"],
        "promoter": promoter_methylation_regions(
            tables["promoters"],
            upstream_bp=settings.get("promoter_upstream_bp", 2000),
            downstream_bp=settings.get("promoter_downstream_bp", 500),
        ),
    }


def _aggregate_annotations(rows, samples):
    """Use the same technical/biological-replicate/donor hierarchy as core activity.

    Retain failed observations as provenance without letting their stored values
    enter the mean. Imported features have no implicit sample pooling contract.
    """
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["entity_type"], row["entity_id"], row["feature_name"]].append(row)
    output = []
    for key, group in grouped.items():
        if len(group) == 1:
            output.append(group[0])
            continue
        if any(not row.get("sample_id") for row in group):
            raise PaceError(
                f"Duplicate annotation without an explicit sample aggregation contract: {key}"
            )
        contracts = {
            (r.get("assay"), r.get("unit"), r.get("normalization_id"), r.get("window_id"))
            for r in group
        }
        if len(contracts) != 1:
            raise PaceError(f"Incompatible measurement contracts for annotation {key}")
        value, used_samples = aggregate_observations(
            [
                {
                    **row,
                    "measurement_status": "observed"
                    if row["status"] in ("observed", "resolved")
                    else "invalid",
                }
                for row in group
            ],
            samples,
            "value",
        )
        donors = {samples[sample]["donor_id"] for sample in used_samples}
        biological_replicates = {
            (samples[sample]["donor_id"], samples[sample]["biological_replicate"])
            for sample in used_samples
        }
        parents = sorted({r["evidence_id"] for r in group})
        output.append(
            {
                "entity_type": key[0],
                "entity_id": key[1],
                "feature_name": key[2],
                "value": value,
                "status": "resolved" if math.isfinite(value) else "unresolved",
                "evidence_id": "annotation_aggregate:" + digest([key, parents]),
                "parent_evidence_ids": ";".join(parents),
                "aggregation_method": "equal_donor_mean",
                "n_observed_samples": len(used_samples),
                "n_observed_biological_replicates": len(biological_replicates),
                "n_observed_donors": len(donors),
                "role": "annotation_only",
            }
        )
    return output


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
                    "evidence_id": observation_evidence_id(row),
                    "status": row["measurement_status"],
                    "role": "annotation_only",
                    **{
                        k: row[k]
                        for k in ("sample_id", "assay", "unit", "normalization_id", "window_id")
                    },
                }
            )
    if tables["methylation"]:
        settings = cfg.get("methylation", {})
        denominators = load_reference_cpg(settings.get("reference_cpg_path"))
        for entity_type, regions in methylation_regions(tables, cfg).items():
            reference_cpg = {
                identifier: count
                for (kind, identifier), count in denominators.items()
                if kind == entity_type
            }
            for row in summarize_methylation(
                tables["methylation"],
                regions,
                minimum_coverage=settings.get("minimum_coverage", 1),
                reference_cpg=reference_cpg,
            ):
                for field in (
                    "M_site",
                    "M_pooled",
                    "covered_cpg",
                    "total_reads",
                    "cpg_coverage_fraction",
                ):
                    prefix = f"DNA_methylation:{row['assay']}"
                    rows.append(
                        {
                            "entity_type": entity_type,
                            "entity_id": row["element_id"],
                            "feature_name": f"{prefix}:{field}",
                            "value": row[field],
                            "evidence_id": methylation_evidence_id(
                                row["sample_id"], row["assay"], entity_type, row["element_id"]
                            ),
                            "status": row["status"] if math.isfinite(row[field]) else "unresolved",
                            "role": "annotation_only",
                            "sample_id": row["sample_id"],
                            "assay": row["assay"],
                        }
                    )
    for row in tables["expression"]:
        rows.append(
            {
                "entity_type": "gene",
                "entity_id": row["gene_id"],
                "feature_name": "RNA:TPM",
                "value": row["tpm"] if row["status"] == "observed" else math.nan,
                "evidence_id": expression_evidence_id(row),
                "status": row["status"],
                "role": "annotation_only",
                "sample_id": row["sample_id"],
                "assay": "RNA",
            }
        )
    rows = _aggregate_annotations(rows, {r["sample_id"]: r for r in tables["samples"]})
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
            parents = []
            for p in promoters:
                f = pvalues.get((p["promoter_id"], feature))
                if f:
                    parents.append(f["evidence_id"])
                if p["pi"] > 0 and (
                    not f
                    or f["status"] not in ("observed", "resolved")
                    or not math.isfinite(f["value"])
                ):
                    missing += p["pi"]
                elif p["pi"] > 0:
                    value += p["pi"] * f["value"]
            rows.append(
                {
                    "entity_type": "gene",
                    "entity_id": gene,
                    "feature_name": f"promoter:{feature}",
                    "value": value if missing == 0 else math.nan,
                    "evidence_id": "pi_weighted_annotation:" + digest([gene, feature]),
                    "parent_evidence_ids": ";".join(sorted(set(parents))) or None,
                    "aggregation_method": "fixed_pi_weighted_annotation",
                    "status": "resolved" if missing == 0 else "unresolved",
                    "missing_pi_mass": missing,
                    "role": "annotation_only",
                }
            )
    roles = {layer: "disabled" for layer in LAYERS}
    if cfg["contact"]["mode"] != "prior_only" and (
        tables["observed_contacts"]
        or any(
            r["evidence_type"] in ("observed", "aggregate", "fused")
            for r in tables["resolved_contacts"]
        )
    ):
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


def evidence_catalog(tables, resolved_activity, resolved_contacts, assets, cfg, *, features=None):
    """Publish resolvable parent IDs for original observations and derived estimates."""
    sources = {r["source_id"]: dict(r) for r in tables["sources"]}
    evidence = {r["evidence_id"]: dict(r) for r in tables["evidence"]}
    samples = {r["sample_id"]: r for r in tables["samples"]}
    for r in tables["observed_activity"]:
        key = observation_evidence_id(r)
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
    for r in tables["expression"]:
        key = expression_evidence_id(r)
        evidence[key] = {
            "evidence_id": key,
            "evidence_type": "observed",
            "source_id": samples[r["sample_id"]]["source_id"],
            "parent_evidence_ids": None,
            "model_id": None,
            "unit": "TPM",
            "processing_method": "input_expression_with_explicit_status",
            "checksum": digest(r),
        }
    if tables["methylation"]:
        count_groups = defaultdict(list)
        for row in tables["methylation"]:
            count_groups[row["sample_id"], row["assay"]].append(row)
        regions = methylation_regions(tables, cfg)
        reference_cpg = load_reference_cpg(cfg.get("methylation", {}).get("reference_cpg_path"))
        for (sample, assay), counts in count_groups.items():
            parent = f"methylation_counts:{sample}:{assay}"
            evidence[parent] = {
                "evidence_id": parent,
                "evidence_type": "observed",
                "source_id": samples[sample]["source_id"],
                "parent_evidence_ids": None,
                "model_id": None,
                "unit": "CpG_dyad_read_counts",
                "processing_method": "input_merged_dyad_counts",
                "checksum": digest(
                    sorted(counts, key=lambda r: (r["chrom"], int(r["dyad_start0"])))
                ),
            }
            for entity_type, windows in regions.items():
                for region in windows:
                    key = methylation_evidence_id(sample, assay, entity_type, region["element_id"])
                    evidence[key] = {
                        "evidence_id": key,
                        "evidence_type": "aggregate",
                        "source_id": generated_source,
                        "parent_evidence_ids": parent,
                        "model_id": None,
                        "unit": "methylation_summary",
                        "processing_method": "CpG_site_and_pooled_proportion",
                        "checksum": digest(
                            {
                                "counts": evidence[parent]["checksum"],
                                "region": region,
                                "minimum_coverage": cfg.get("methylation", {}).get(
                                    "minimum_coverage", 1
                                ),
                                "reference_cpg": reference_cpg.get(
                                    (entity_type, region["element_id"])
                                ),
                            }
                        ),
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
    feature_groups = defaultdict(list)
    for row in features or []:
        feature_groups[row["evidence_id"]].append(row)
    for key, rows in feature_groups.items():
        if key in evidence:
            continue
        parents = sorted(
            {
                parent
                for row in rows
                for parent in (row.get("parent_evidence_ids") or "").split(";")
                if parent
            }
        )
        if key == "core_formula":
            parents = sorted({r["evidence_id"] for r in resolved_activity + resolved_contacts})
        evidence[key] = {
            "evidence_id": key,
            "evidence_type": "aggregate",
            "source_id": generated_source,
            "parent_evidence_ids": ";".join(parents) or None,
            "model_id": None,
            "unit": "derived_feature",
            "processing_method": rows[0].get("aggregation_method", "core_formula"),
            "checksum": digest(rows),
        }
    for row in evidence.values():
        if row["source_id"] not in sources:
            raise PaceError(f"Generated evidence references an unknown source: {row['source_id']}")
        for parent in (row.get("parent_evidence_ids") or "").split(";"):
            if parent and parent not in evidence:
                raise PaceError(f"Generated evidence references an unknown parent: {parent}")
    return sorted(evidence.values(), key=lambda r: r["evidence_id"]), sorted(
        sources.values(), key=lambda r: r["source_id"]
    )
