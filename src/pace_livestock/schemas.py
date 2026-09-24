"""Standard tables and cross-table scientific integrity constraints."""

from __future__ import annotations

import math
from collections import defaultdict

from .errors import PaceError
from .io.tables import integer, number, read_table, unique
from .provenance import digest, file_hash

SCHEMAS = {
    "units": "element_id chrom start end anchor0 element_roles canonical_catalog_id",
    "region_membership": "region_id element_id source_id membership_rule",
    "promoters": "gene_id promoter_id chrom tss0 strand pi pi_source",
    "candidates": "element_id gene_id candidate_universe_id",
    "samples": "sample_id donor_id assay biological_replicate technical_replicate species assembly context_id source_id",
    "observed_activity": "element_id sample_id assay signal measurement_status callable_fraction unit normalization_id window_id",
    "observed_contacts": "element_id promoter_id sample_id contact_value measurement_status bin_pair_id scale resolution source_id",
    "resolved_activity": "element_id assay observed_value resolved_value evidence_id evidence_type observation_sample_id parent_evidence_ids resolution_status reason unit normalization_id window_id",
    "resolved_contacts": "element_id promoter_id resolved_value evidence_id evidence_type observation_sample_id prior_id reliability resolved_mode bin_pair_id resolution_status reason scale",
    "features": "entity_type entity_id feature_name value evidence_id status",
    "methylation": "chrom dyad_start0 methylated_count total_count sample_id assay",
    "expression": "gene_id sample_id tpm status",
    "labels": "label_id assayed_region_id gene_id context_id perturbation_type effect_direction effect_size label_status assay_id group_id source_id",
    "evidence": "evidence_id evidence_type source_id parent_evidence_ids model_id unit processing_method checksum",
    "sources": "source_id path_or_accession source_type assembly processing_method normalization_id checksum",
    "support_bounds": "element_id gene_id support_lower support_upper bound_source",
}
STATUSES = {"observed", "unmeasured", "low_coverage", "unmappable", "invalid", "not_applicable"}
EVIDENCE_TYPES = {"observed", "contact_prior", "fused", "aggregate", "regularized"}
CONTACT_ASSAYS = {
    "Hi-C",
    "HiC",
    "Micro-C",
    "Capture-C",
    "Capture-Hi-C",
    "PCHi-C",
    "HiChIP",
    "PLAC-seq",
}


def load_tables(cfg: dict) -> dict[str, list[dict]]:
    tables = {
        name: read_table(path, required=SCHEMAS[name].split()) if path else []
        for name, path in cfg["inputs"].items()
    }
    for name in ("units", "promoters", "candidates"):
        if not tables[name]:
            raise PaceError(f"inputs.{name}: a nonempty table is required")
    tables["inferred_metadata"] = infer_metadata(tables, cfg)
    validate_tables(tables, cfg)
    return tables


def infer_metadata(tables: dict, cfg: dict) -> dict[str, list[dict]]:
    """Fill omitted samples/sources tables for the common single-animal case.

    Every sample becomes a separate biological replicate of one animal. Nothing is
    inferred when the table is supplied, and population-level runs must declare
    which animal each library came from.
    """
    inferred = {}
    measured = [
        (name, row)
        for name in ("observed_activity", "observed_contacts", "expression", "methylation")
        for row in tables[name]
    ]
    if not cfg["inputs"]["samples"] and measured:
        if cfg["target_level"] != "individual":
            raise PaceError(
                "target_level=population_mean needs --samples: a table giving the donor_id "
                "(animal) of every sample_id"
            )
        assays, sources = {}, {}
        for name, row in measured:
            assay = {
                "observed_contacts": "HiC",
                "expression": "RNA-seq",
            }.get(name, row.get("assay"))
            sample = row["sample_id"]
            if sample in assays and assays[sample] != assay:
                raise PaceError(f"sample_id {sample} is used for both {assays[sample]} and {assay}")
            assays[sample] = assay
            sources.setdefault(sample, row.get("source_id") or sample)
        replicate = defaultdict(int)
        rows = []
        for sample in sorted(assays):
            replicate[assays[sample]] += 1
            rows.append(
                {
                    "sample_id": sample,
                    "donor_id": "individual_1",
                    "assay": assays[sample],
                    "biological_replicate": str(replicate[assays[sample]]),
                    "technical_replicate": "1",
                    "species": cfg["context"]["species"],
                    "assembly": cfg["context"]["assembly"],
                    "context_id": cfg["context"]["context_id"],
                    "source_id": sources[sample],
                    "metadata_origin": "inferred_single_individual",
                }
            )
        tables["samples"] = rows
        inferred["samples"] = rows
    if not cfg["inputs"]["sources"]:
        referenced = {}
        origin = {
            "samples": None,
            "observed_contacts": "observed_contacts",
            "region_membership": "region_membership",
            "evidence": "evidence",
        }
        for name, input_name in origin.items():
            for row in tables[name]:
                if row.get("source_id"):
                    referenced.setdefault(row["source_id"], input_name)
        normalizations = defaultdict(set)
        by_sample = {r["sample_id"]: r["source_id"] for r in tables["samples"]}
        for name in ("observed_activity", "observed_contacts"):
            for row in tables[name]:
                source = row.get("source_id") or by_sample.get(row["sample_id"])
                if source:
                    normalizations[source].add(row.get("normalization_id"))
        rows = []
        for source_id, input_name in sorted(referenced.items()):
            norms = normalizations[source_id] - {None}
            path = cfg["inputs"].get(input_name) if input_name else None
            if input_name is None:
                sample_tables = [
                    n
                    for n in ("observed_activity", "expression", "methylation")
                    if any(r["sample_id"] == source_id for r in tables[n])
                ]
                path = cfg["inputs"].get(sample_tables[0]) if sample_tables else None
            rows.append(
                {
                    "source_id": source_id,
                    "path_or_accession": path or "not_declared",
                    "source_type": "not_declared",
                    "assembly": cfg["context"]["assembly"],
                    "processing_method": "not_declared",
                    "normalization_id": next(iter(norms)) if len(norms) == 1 else None,
                    "checksum": None,
                }
            )
        tables["sources"] = rows
        if rows:
            inferred["sources"] = rows
    return inferred


def validate_tables(t: dict, cfg: dict) -> None:
    keys = {
        "units": ("element_id",),
        "promoters": ("gene_id", "promoter_id"),
        "candidates": ("element_id", "gene_id"),
        "samples": ("sample_id",),
        "sources": ("source_id",),
        "evidence": ("evidence_id",),
        "observed_activity": ("element_id", "sample_id", "assay"),
        "observed_contacts": ("element_id", "promoter_id", "sample_id"),
        "resolved_activity": ("element_id", "assay"),
        "resolved_contacts": ("element_id", "promoter_id"),
        "expression": ("gene_id", "sample_id"),
        "methylation": ("chrom", "dyad_start0", "sample_id", "assay"),
        "features": ("entity_type", "entity_id", "feature_name"),
        "labels": ("label_id",),
        "support_bounds": ("element_id", "gene_id"),
    }
    for name, key in keys.items():
        unique(t[name], key, name)
    units = {r["element_id"]: r for r in t["units"]}
    genes = {r["gene_id"] for r in t["promoters"]}
    samples = {r["sample_id"]: r for r in t["samples"]}
    sources = {r["source_id"] for r in t["sources"]}
    evidence_rows = {r["evidence_id"]: r for r in t["evidence"]}
    evidence = set(evidence_rows)
    chrom_sizes = {}
    if cfg["catalog"].get("chrom_sizes_path"):
        size_rows = read_table(cfg["catalog"]["chrom_sizes_path"], required=["chrom", "length"])
        unique(size_rows, ("chrom",), "chromosome sizes")
        chrom_sizes = {
            r["chrom"]: integer(r["length"], "chromosome length", minimum=1) for r in size_rows
        }
    promoter_coords = {}
    gene_coords, pi_sums = defaultdict(set), defaultdict(float)
    for row in t["units"]:
        for field in ("start", "end", "anchor0"):
            row[field] = integer(row[field], f"units.{field}")
        start, end = row["start"], row["end"]
        if start >= end or row["anchor0"] != (start + end - 1) // 2:
            raise PaceError(f"units {row['element_id']}: invalid half-open interval or anchor")
        if chrom_sizes and (row["chrom"] not in chrom_sizes or end > chrom_sizes[row["chrom"]]):
            raise PaceError("Unit lies outside the declared chromosome bounds")
        if cfg["catalog"]["profile"] == "canonical_grid":
            width, offset = cfg["catalog"]["width_bp"], cfg["catalog"]["offset_bp"]
            if end - start != width or (start - offset) % width:
                raise PaceError(
                    f"units {row['element_id']}: interval does not match canonical grid"
                )
    unique(t["units"], ("chrom", "start", "end"), "unit coordinates")
    if len({r["canonical_catalog_id"] for r in t["units"]}) != 1:
        raise PaceError("units: mixed canonical_catalog_id values")
    for row in t["promoters"]:
        row["tss0"] = integer(row["tss0"], "promoters.tss0")
        if chrom_sizes and (
            row["chrom"] not in chrom_sizes or row["tss0"] >= chrom_sizes[row["chrom"]]
        ):
            raise PaceError("Promoter TSS lies outside the declared chromosome bounds")
        row["pi"] = number(row["pi"], "promoters.pi", minimum=0, maximum=1)
        if row["strand"] not in ("+", "-"):
            raise PaceError("promoters.strand must be + or -")
        coord = row["chrom"], row["tss0"], row["strand"]
        if row["promoter_id"] in promoter_coords and promoter_coords[row["promoter_id"]] != coord:
            raise PaceError("A physical promoter_id has inconsistent coordinates")
        promoter_coords[row["promoter_id"]] = coord
        if coord in gene_coords[row["gene_id"]]:
            raise PaceError("Duplicate physical TSS within a gene; deduplicate transcripts first")
        gene_coords[row["gene_id"]].add(coord)
        pi_sums[row["gene_id"]] += row["pi"]
    if cfg["promoters"]["weights"] == "equal":
        for row in t["promoters"]:
            row["pi"], row["pi_source"] = 1 / len(gene_coords[row["gene_id"]]), "equal_physical_tss"
    elif any(not math.isclose(x, 1, abs_tol=1e-10) for x in pi_sums.values()):
        raise PaceError(
            "Promoter pi values must sum to one per gene before any missing-data filtering"
        )
    if any(len({c[0] for c in coords}) != 1 for coords in gene_coords.values()):
        raise PaceError("Gene promoters on multiple chromosomes are unsupported")
    for row in t["candidates"]:
        if row["element_id"] not in units or row["gene_id"] not in genes:
            raise PaceError(f"Unknown candidate unit or gene: {row}")
        if units[row["element_id"]]["chrom"] != next(iter(gene_coords[row["gene_id"]]))[0]:
            raise PaceError("Only cis candidate links are supported")
        distance = min(
            abs(units[row["element_id"]]["anchor0"] - coord[1])
            for coord in gene_coords[row["gene_id"]]
        )
        if distance > cfg["catalog"].get("candidate_radius_bp", 5_000_000):
            raise PaceError("Candidate lies outside catalog.candidate_radius_bp")
    if len({r["candidate_universe_id"] for r in t["candidates"]}) != 1:
        raise PaceError("candidates: mixed candidate_universe_id values")
    candidate_keys = {(r["element_id"], r["gene_id"]) for r in t["candidates"]}
    for row in t["support_bounds"]:
        if (row["element_id"], row["gene_id"]) not in candidate_keys or not row["bound_source"]:
            raise PaceError("Support bounds require a known candidate and bound_source")
        row["support_lower"] = number(row["support_lower"], "support_lower", minimum=0)
        row["support_upper"] = number(
            row["support_upper"], "support_upper", minimum=0, missing=True
        )
        if row["support_upper"] < row["support_lower"]:
            raise PaceError("support_upper must be at least support_lower")
    for row in t["samples"] + t["sources"]:
        if row["source_id"] not in sources:
            raise PaceError(f"Unknown source_id: {row['source_id']}")
        for field in ("species", "assembly", "context_id"):
            if field in row and row[field] != cfg["context"][field]:
                raise PaceError(f"{field} mismatch for sample/source {row}")
    for name in ("observed_activity", "observed_contacts", "methylation", "expression"):
        for row in t[name]:
            if row["sample_id"] not in samples:
                raise PaceError(f"{name}: unknown real sample_id {row['sample_id']}")
            if "element_id" in row and row["element_id"] not in units:
                raise PaceError(f"{name}: unknown element_id {row['element_id']}")
            if "assay" in row and row["assay"] != samples[row["sample_id"]]["assay"]:
                raise PaceError(f"{name}: assay differs from samples table")
    for name, value_col in (
        ("observed_activity", "signal"),
        ("observed_contacts", "contact_value"),
    ):
        for row in t[name]:
            if row["measurement_status"] not in STATUSES:
                raise PaceError(f"{name}: invalid measurement_status")
            row[value_col] = number(row[value_col], f"{name}.{value_col}", missing=True, minimum=0)
            if row["measurement_status"] == "observed" and math.isnan(row[value_col]):
                raise PaceError(f"{name}: observed row requires an explicit finite value")
            if name == "observed_activity":
                row["callable_fraction"] = number(
                    row["callable_fraction"], "callable_fraction", minimum=0, maximum=1
                )
                if not all(row[k] for k in ("unit", "normalization_id", "window_id")):
                    raise PaceError("Activity requires unit, normalization_id and window_id")
            else:
                for field in (
                    "raw_count",
                    "count_to_contact",
                    "balance_weight_1",
                    "balance_weight_2",
                ):
                    if row.get(field) is not None:
                        row[field] = number(row[field], field, minimum=0)
                if row.get("raw_count") is not None and row["raw_count"] != math.floor(
                    row["raw_count"]
                ):
                    raise PaceError("raw_count must contain integer counts")
                if (
                    cfg["contact"]["reliability"] == "per_pair"
                    and row["measurement_status"] == "observed"
                ):
                    if row.get("raw_count") is None or not row.get("count_to_contact"):
                        raise PaceError("per_pair requires raw_count and positive count_to_contact")
                    if not math.isclose(
                        row[value_col],
                        row["raw_count"] * row["count_to_contact"],
                        rel_tol=1e-8,
                        abs_tol=1e-12,
                    ):
                        raise PaceError("Raw counts/factor do not reproduce contact_value")
                if row.get("near_diagonal_value") is not None:
                    row["near_diagonal_value"] = number(
                        row["near_diagonal_value"], "near_diagonal_value", minimum=0, missing=True
                    )
                    if (
                        math.isfinite(row["near_diagonal_value"])
                        and row.get("near_diagonal_method") != "neighbor_max"
                    ):
                        raise PaceError(
                            "Near-diagonal replacement requires a declared neighbor_max method"
                        )
                if samples[row["sample_id"]]["assay"] not in CONTACT_ASSAYS:
                    raise PaceError(
                        "Contact observation requires a chromosome-contact assay sample"
                    )
                if row["promoter_id"] not in promoter_coords or row["source_id"] not in sources:
                    raise PaceError("Contact has unknown promoter/source")
                row["resolution"] = integer(row["resolution"], "contact.resolution", minimum=1)
                if row["scale"] != cfg["contact"]["scale"]:
                    raise PaceError("Contact scale differs from run configuration")
    for name in ("resolved_activity", "resolved_contacts"):
        for row in t[name]:
            if row["element_id"] not in units or row["evidence_id"] not in evidence:
                raise PaceError(f"{name}: unknown unit or evidence_id")
            if row["evidence_type"] not in EVIDENCE_TYPES:
                raise PaceError(f"{name}: invalid evidence_type")
            declared = evidence_rows[row["evidence_id"]]
            if declared["evidence_type"] != row["evidence_type"]:
                raise PaceError("Resolved evidence type disagrees with evidence catalog")
            for parent in (row.get("parent_evidence_ids") or "").split(";"):
                if parent and parent not in evidence:
                    raise PaceError(f"Unknown resolved parent evidence: {parent}")
            if row["evidence_type"] == "aggregate" and not row.get("parent_evidence_ids"):
                raise PaceError("Imported aggregate evidence requires its parent evidence IDs")
            if name == "resolved_contacts" and row["promoter_id"] not in promoter_coords:
                raise PaceError("Imported contact has unknown promoter_id")
            sample = row["observation_sample_id"]
            if sample and sample not in samples:
                raise PaceError(f"{name}: unknown observation_sample_id")
            if row["evidence_type"] == "observed" and not sample:
                raise PaceError("Observed evidence requires a real sample")
            if row["evidence_type"] == "regularized" and not (
                sample or row.get("parent_evidence_ids")
            ):
                raise PaceError("Regularized contacts require measured parent evidence")
            if row["evidence_type"] == "contact_prior" and sample:
                raise PaceError("Contact prior evidence cannot impersonate a sample")
            if name == "resolved_activity":
                if row["evidence_type"] not in ("observed", "aggregate"):
                    raise PaceError("Activity imports must contain measured evidence")
                if sample and samples[sample]["assay"] != row["assay"]:
                    raise PaceError("Imported activity assay differs from samples table")
            elif sample and samples[sample]["assay"] not in CONTACT_ASSAYS:
                raise PaceError("Imported contact requires a chromosome-contact assay sample")
            if row["resolution_status"] not in ("resolved", "unresolved", "invalid"):
                raise PaceError("Invalid resolution_status")
            row["resolved_value"] = number(
                row["resolved_value"], "resolved_value", missing=True, minimum=0
            )
            if row["resolution_status"] == "resolved" and math.isnan(row["resolved_value"]):
                raise PaceError("Resolved evidence requires a finite value")
    for row in t["evidence"]:
        if row["source_id"] not in sources or row["evidence_type"] not in EVIDENCE_TYPES:
            raise PaceError("Evidence requires a known source and a legal type")
        for parent in (row["parent_evidence_ids"] or "").split(";"):
            if parent and parent not in evidence:
                raise PaceError(f"Unknown parent_evidence_id {parent}")
    for row in t["expression"]:
        if row["gene_id"] not in genes:
            raise PaceError("Expression has unknown gene_id")
        if samples[row["sample_id"]]["assay"] not in ("RNA", "RNA-seq"):
            raise PaceError("Expression sample requires assay RNA or RNA-seq")
        if row["status"] not in STATUSES:
            raise PaceError("Expression has invalid status")
        row["tpm"] = number(row["tpm"], "expression.tpm", minimum=0, missing=True)
        if row["status"] == "observed" and not math.isfinite(row["tpm"]):
            raise PaceError("Observed expression requires a finite TPM")
        if row["status"] != "observed":
            row["tpm"] = math.nan
    known_features = {
        "element": set(units),
        "promoter": set(promoter_coords),
        "gene": genes,
        "edge": {f"{r['element_id']}|{r['gene_id']}" for r in t["candidates"]},
    }
    for row in t["features"]:
        row["value"] = number(row["value"], "features.value", missing=True)
        if row["evidence_id"] not in evidence:
            raise PaceError("Feature has unknown evidence_id")
        if row["entity_type"] not in ("element", "promoter", "gene", "edge"):
            raise PaceError("Feature entity_type must be element, promoter, gene or edge")
        if row["status"] not in STATUSES | {"resolved", "unresolved"}:
            raise PaceError("Invalid feature status")
        if row["status"] not in ("observed", "resolved"):
            row["value"] = math.nan
        elif math.isnan(row["value"]):
            raise PaceError("Observed feature requires a finite value")
        if row["entity_id"] not in known_features[row["entity_type"]]:
            raise PaceError("Feature entity_id is outside the declared catalog")
    for row in t["region_membership"]:
        if row["element_id"] not in units or row["source_id"] not in sources:
            raise PaceError("Region membership has an unknown element or source")
    if cfg["catalog"]["include_promoter_units"]:
        cells = {(r["chrom"], r["start"], r["end"]) for r in t["units"]}
        width, offset = cfg["catalog"]["width_bp"], cfg["catalog"]["offset_bp"]
        if cfg["catalog"]["profile"] == "canonical_grid":
            from .catalog import canonical_cell_bounds

            for row in t["promoters"]:
                start, end = canonical_cell_bounds(row["tss0"], width=width, offset=offset)
                if row["chrom"] in chrom_sizes and (start < 0 or end > chrom_sizes[row["chrom"]]):
                    continue
                if (row["chrom"], start, end) not in cells:
                    raise PaceError(
                        "Promoter unit absent: prepare the catalog and provide catalog.chrom_sizes_path for valid boundary exclusions, or explicitly set include_promoter_units=false"
                    )


def universe_ids(t: dict, cfg: dict) -> dict:
    return {
        "canonical_catalog_id": digest(
            {
                "definition": {k: v for k, v in cfg["catalog"].items() if k != "chrom_sizes_path"},
                "reference_sizes_sha256": file_hash(cfg["catalog"]["chrom_sizes_path"])
                if cfg["catalog"].get("chrom_sizes_path")
                else None,
                "units": sorted(
                    (r["element_id"], r["chrom"], r["start"], r["end"]) for r in t["units"]
                ),
            }
        ),
        "candidate_universe_id": digest(
            sorted((r["element_id"], r["gene_id"]) for r in t["candidates"])
        ),
        "promoter_universe_id": digest(
            sorted(
                (r["gene_id"], r["promoter_id"], r["chrom"], r["tss0"], r["strand"], r["pi"])
                for r in t["promoters"]
            )
        ),
    }
