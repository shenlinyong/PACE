"""Configuration-driven fitting and genomic preparation commands."""

from __future__ import annotations

from pathlib import Path

import yaml

from .catalog import candidate_edges, canonical_units
from .config import operation_config
from .errors import PaceError
from .evidence.contact import distance_prior, fit_distance_prior
from .io.tables import integer, number, read_table, write_table
from .provenance import file_hash, output_directory, write_json

ASSET_KEYS = {"model_id", "species", "assembly", "context_id", "target_level", "is_synthetic"}


def asset_metadata(cfg, kind):
    if type(cfg["is_synthetic"]) is not bool:
        raise PaceError("is_synthetic must be an explicit YAML boolean")
    return {
        **{k: cfg[k] for k in ASSET_KEYS},
        "kind": kind,
        "training_sources": [file_hash(cfg["data"])],
        "calibration_sources": [],
        "test_sources": [],
        "validation": {},
    }


def fit_contact_command(path, out):
    required = ASSET_KEYS | {"data", "scale", "resolution", "bin_edges", "d_ref", "d_min"}
    cfg = operation_config(
        path,
        allowed=required | {"normalization_id", "balancing", "window_id"},
        required=required,
        paths=["data"],
    )
    rows = read_table(
        cfg["data"], required=["bin_pair_id", "distance_bp", "contact_value", "split", "region_id"]
    )
    seen, region_splits = set(), {}
    for row in rows:
        if row["bin_pair_id"] in seen:
            raise PaceError("Prior fitting must use each measured bin pair once")
        seen.add(row["bin_pair_id"])
        if row["split"] not in ("train", "test"):
            raise PaceError("Contact prior split must be train or test")
        if row["region_id"] in region_splits and region_splits[row["region_id"]] != row["split"]:
            raise PaceError("Contact prior training region leaks into held-out evaluation")
        region_splits[row["region_id"]] = row["split"]
    train = [r for r in rows if r["split"] == "train"]
    d = [number(r["distance_bp"], "distance_bp", minimum=0) for r in train]
    c = [number(r["contact_value"], "contact_value", minimum=0) for r in train]
    fitted = fit_distance_prior(
        d, c, bin_edges=cfg["bin_edges"], d_ref=cfg["d_ref"], d_min=cfg["d_min"]
    )
    manifest = {
        **asset_metadata(cfg, "contact_prior"),
        **fitted,
        "scale": cfg["scale"],
        "resolution": integer(cfg["resolution"], "resolution", minimum=1),
        **{key: cfg.get(key) for key in ("normalization_id", "balancing", "window_id")},
        "training_regions": sorted(r for r, split in region_splits.items() if split == "train"),
    }
    residuals = [
        {
            "bin_pair_id": r["bin_pair_id"],
            "observed": number(r["contact_value"], "contact_value", minimum=0),
            "predicted": distance_prior(
                number(r["distance_bp"], "distance_bp", minimum=0), manifest
            ),
        }
        for r in rows
        if r["split"] == "test"
    ]
    with output_directory(out) as dest:
        write_json(dest / "manifest.json", manifest)
        write_table(
            dest / "held_out_residuals.tsv",
            [{**r, "residual": r["observed"] - r["predicted"]} for r in residuals],
            fields=["bin_pair_id", "observed", "predicted", "residual"],
        )
    return manifest


def prepare_command(path, out):
    # A separate operation schema prevents silently ignoring options for another adapter.
    from .config import load_yaml

    kind = load_yaml(path).get("kind")
    if kind == "catalog":
        cfg = operation_config(
            path,
            allowed={
                "kind",
                "bed",
                "gtf",
                "chrom_sizes",
                "source_id",
                "width",
                "offset",
                "radius",
                "include_promoters",
                "aliases",
            },
            required=["bed", "gtf", "chrom_sizes", "source_id"],
            paths=["bed", "gtf", "chrom_sizes", "aliases"],
        )
        from .io.bed_gtf import read_bed, read_gtf

        aliases = (
            {
                r["alias"]: r["canonical"]
                for r in read_table(cfg["aliases"], required=["alias", "canonical"])
            }
            if cfg.get("aliases")
            else {}
        )
        sizes = {
            aliases.get(r["chrom"], r["chrom"]): integer(r["length"], "chrom length", minimum=1)
            for r in read_table(cfg["chrom_sizes"], required=["chrom", "length"])
        }
        promoters, transcripts = read_gtf(cfg["gtf"], aliases=aliases)
        regions = read_bed(cfg["bed"], source_id=cfg["source_id"], aliases=aliases)
        units, memberships, excluded = canonical_units(
            regions,
            sizes,
            promoters,
            width=cfg.get("width", 500),
            offset=cfg.get("offset", 0),
            include_promoters=cfg.get("include_promoters", True),
        )
        candidates = candidate_edges(units, promoters, radius=cfg.get("radius", 5_000_000))
        with output_directory(out) as dest:
            for name, rows in (
                ("units", units),
                ("promoters", promoters),
                ("region_membership", memberships),
                ("candidates", candidates),
                ("transcript_mapping", transcripts),
            ):
                write_table(dest / f"{name}.tsv", rows)
            write_table(
                dest / "chrom_sizes.tsv",
                [{"chrom": chrom, "length": length} for chrom, length in sorted(sizes.items())],
            )
            (dest / "run_catalog_config.yaml").write_text(
                yaml.safe_dump(
                    {
                        "catalog": {
                            "profile": "canonical_grid",
                            "width_bp": cfg.get("width", 500),
                            "offset_bp": cfg.get("offset", 0),
                            "include_promoter_units": cfg.get("include_promoters", True),
                            "chrom_sizes_path": "chrom_sizes.tsv",
                            "candidate_radius_bp": cfg.get("radius", 5_000_000),
                        },
                        "inputs": {
                            name: f"{name}.tsv"
                            for name in ("units", "promoters", "region_membership", "candidates")
                        },
                    },
                    sort_keys=False,
                ),
                encoding="utf-8",
            )
            write_json(
                dest / "preparation_report.json",
                {
                    "excluded_units": excluded,
                    "assembly_inferred": False,
                    "reference_sizes_sha256": file_hash(cfg["chrom_sizes"]),
                    "gtf_sha256": file_hash(cfg["gtf"]),
                    "bed_sha256": file_hash(cfg["bed"]),
                },
            )
        return {"units": len(units), "candidates": len(candidates)}
    if kind == "bigwig":
        cfg = operation_config(
            path,
            allowed={
                "kind",
                "track",
                "units",
                "sample_id",
                "assay",
                "unit",
                "normalization_id",
                "window_id",
                "missing_is_measured_zero",
                "minimum_callable_fraction",
            },
            required=[
                "track",
                "units",
                "sample_id",
                "assay",
                "unit",
                "normalization_id",
                "window_id",
            ],
            paths=["track", "units"],
        )
        from .io.bigwig import quantify_bigwig

        units = read_table(cfg["units"], required=["element_id", "chrom", "start", "end"])
        rows = quantify_bigwig(
            cfg["track"],
            units,
            missing_is_measured_zero=cfg.get("missing_is_measured_zero", False),
            minimum_callable_fraction=cfg.get("minimum_callable_fraction", 0),
        )
        rows = [
            {
                **r,
                **{
                    k: cfg[k]
                    for k in ("sample_id", "assay", "unit", "normalization_id", "window_id")
                },
            }
            for r in rows
        ]
        name = "observed_activity.tsv"
    elif kind == "cooler":
        cfg = operation_config(
            path,
            allowed={
                "kind",
                "contact",
                "pairs",
                "resolution",
                "balanced",
                "missing_pixels_are_zero",
                "sample_id",
                "source_id",
                "scale",
                "normalization_id",
            },
            required=[
                "contact",
                "pairs",
                "resolution",
                "balanced",
                "missing_pixels_are_zero",
                "sample_id",
                "source_id",
                "scale",
            ],
            paths=["pairs"],
        )
        from .io.cooler import query_contacts

        file, *group = cfg["contact"].split("::", 1)
        uri = str((Path(path).resolve().parent / file).resolve()) + (
            "::" + group[0] if group else ""
        )
        pairs = read_table(
            cfg["pairs"], required=["element_id", "promoter_id", "chrom", "anchor0", "tss0"]
        )
        for r in pairs:
            r["anchor0"], r["tss0"] = integer(r["anchor0"], "anchor0"), integer(r["tss0"], "tss0")
        rows = query_contacts(
            uri,
            pairs,
            resolution=cfg["resolution"],
            balanced=cfg["balanced"],
            missing_pixels_are_zero=cfg["missing_pixels_are_zero"],
        )
        rows = [
            {
                **r,
                **{k: cfg[k] for k in ("sample_id", "source_id", "scale")},
                "normalization_id": cfg.get("normalization_id"),
                "balancing": "balanced" if cfg["balanced"] else "unbalanced",
                "window_id": "bin_pair",
            }
            for r in rows
        ]
        name = "observed_contacts.tsv"
    elif kind == "bed_features":
        cfg = operation_config(
            path,
            allowed={
                "kind",
                "bed",
                "units",
                "source_id",
                "evidence_id",
                "feature_prefix",
                "entity_type",
                "motif_strands",
            },
            required=["bed", "units", "source_id", "evidence_id", "feature_prefix"],
            paths=["bed", "units"],
        )
        from .io.bed_gtf import read_bed
        from .io.intervals import interval_features

        units = read_table(cfg["units"], required=["element_id", "chrom", "start", "end"])
        peaks = read_bed(cfg["bed"], source_id=cfg["source_id"])
        if type(cfg.get("motif_strands", False)) is not bool:
            raise PaceError("motif_strands must be a YAML boolean")
        rows = interval_features(
            units,
            peaks,
            evidence_id=cfg["evidence_id"],
            feature_prefix=cfg["feature_prefix"],
            entity_type=cfg.get("entity_type", "element"),
            motif_strands=cfg.get("motif_strands", False),
        )
        name = "features.tsv"
    elif kind == "methylation":
        cfg = operation_config(
            path,
            allowed={"kind", "counts", "units", "minimum_coverage", "reference_cpg"},
            required=["counts", "units"],
            paths=["counts", "units", "reference_cpg"],
        )
        from .io.methylation import summarize_methylation

        units = read_table(cfg["units"], required=["element_id", "chrom", "start", "end"])
        for r in units:
            r["start"], r["end"] = integer(r["start"], "start"), integer(r["end"], "end")
        counts = read_table(
            cfg["counts"],
            required=[
                "chrom",
                "dyad_start0",
                "methylated_count",
                "total_count",
                "sample_id",
                "assay",
            ],
        )
        reference_cpg = (
            {
                r["element_id"]: integer(r["n_cpg"], "n_cpg")
                for r in read_table(cfg["reference_cpg"], required=["element_id", "n_cpg"])
            }
            if cfg.get("reference_cpg")
            else None
        )
        rows = summarize_methylation(
            counts,
            units,
            minimum_coverage=cfg.get("minimum_coverage", 1),
            reference_cpg=reference_cpg,
        )
        name = "methylation_summary.tsv"
    elif kind == "rna":
        cfg = operation_config(
            path,
            allowed={"kind", "expression", "transcript_mapping"},
            required=["expression", "transcript_mapping"],
            paths=["expression", "transcript_mapping"],
        )
        from .io.rna import transcript_tpm_to_gene

        rows = transcript_tpm_to_gene(
            read_table(cfg["expression"], required=["transcript_id", "sample_id", "tpm"]),
            read_table(
                cfg["transcript_mapping"], required=["transcript_id", "gene_id", "promoter_id"]
            ),
        )
        name = "expression.tsv"
    else:
        raise PaceError(
            "prepare.kind must be catalog, bigwig, cooler, bed_features, methylation or rna"
        )
    with output_directory(out) as dest:
        write_table(dest / name, rows)
        write_json(
            dest / "preparation_report.json",
            {"kind": kind, "rows": len(rows), "configuration": cfg},
        )
    return {"rows": len(rows)}
