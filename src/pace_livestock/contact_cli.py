"""Boundary preparation, prior fitting, contact fusion and eQTL calibration commands."""

from __future__ import annotations

import math

from .boundary_prior import (
    BoundaryIndex,
    contact_prior,
    fit_contact_grid,
    motif_boundaries,
    scan_motifs,
    unique_count_pairs,
)
from .config import load_config, operation_config
from .errors import PaceError
from .io.tables import integer, number, read_table, write_table
from .operations import ASSET_KEYS
from .provenance import file_hash, output_directory, write_json

BOUNDARY_FIELDS = ["chrom", "position0", "strength"]
MEASUREMENT_FIELDS = {"scale", "resolution", "normalization_id", "balancing", "window_id"}


def metadata(cfg):
    if type(cfg["is_synthetic"]) is not bool:
        raise PaceError("is_synthetic must be an explicit YAML boolean")
    for key in ASSET_KEYS - {"is_synthetic"}:
        if not isinstance(cfg[key], str) or not cfg[key].strip():
            raise PaceError(f"{key} must be a nonempty string")
    if cfg["target_level"] not in ("individual", "population_mean"):
        raise PaceError("Unknown target_level")
    if cfg["scale"] in ("oe", "O/E", "log_oe", "pvalue", "correlation"):
        raise PaceError("Contact prior scale must retain the distance background")
    return {
        **{k: cfg[k] for k in ASSET_KEYS | MEASUREMENT_FIELDS},
        "resolution": integer(cfg["resolution"], "resolution", minimum=1),
        "kind": "contact_prior",
        "training_sources": [],
        "calibration_sources": [],
        "test_sources": [],
        "validation": {},
    }


def boundaries_command(path, out):
    cfg = operation_config(
        path,
        allowed={"motifs", "fasta", "pwm", "threshold", "chip", "max_gap_bp"},
        paths=["motifs", "fasta", "pwm", "chip"],
    )
    if bool(cfg.get("motifs")) == bool(cfg.get("fasta")):
        raise PaceError("Supply motifs OR fasta plus pwm and threshold")
    if cfg.get("motifs"):
        if cfg.get("pwm") is not None or cfg.get("threshold") is not None:
            raise PaceError("PWM/threshold are only used with fasta")
        motifs = read_table(cfg["motifs"], required=["chrom", "start", "end", "strand", "strength"])
    else:
        if not cfg.get("pwm") or cfg.get("threshold") is None:
            raise PaceError("FASTA scanning requires a CTCF PWM and log2-odds threshold")
        motifs = scan_motifs(
            cfg["fasta"], read_table(cfg["pwm"], required=list("ACGT")), threshold=cfg["threshold"]
        )
    chip = read_table(cfg["chip"], required=["chrom", "start", "end"]) if cfg.get("chip") else None
    boundaries = motif_boundaries(motifs, chip=chip, max_gap_bp=cfg.get("max_gap_bp", 1_000_000))
    with output_directory(out) as dest:
        write_table(
            dest / "motifs.tsv",
            motifs,
            fields=None if motifs else ["chrom", "start", "end", "strand", "strength"],
        )
        write_table(
            dest / "boundaries.tsv", boundaries, fields=None if boundaries else BOUNDARY_FIELDS
        )
        write_json(
            dest / "boundary_report.json",
            dict(
                n_motifs=len(motifs),
                n_boundaries=len(boundaries),
                boundary_status="candidate_unvalidated",
                parameters={k: cfg.get(k) for k in ("threshold", "max_gap_bp")},
                source_hashes={
                    k: file_hash(cfg[k]) for k in ("motifs", "fasta", "pwm", "chip") if cfg.get(k)
                },
            ),
        )
    return {"output": str(out), "n_boundaries": len(boundaries), "status": "candidate_unvalidated"}


def prior_command(path, out, *, fit=False):
    required = ASSET_KEYS | MEASUREMENT_FIELDS | {"d_ref", "d_min"}
    required |= {"data", "gamma_grid", "beta_grid"} if fit else {"a", "gamma", "beta"}
    optional = {"boundaries", "pairs", "test_chromosomes", "kappa"}
    cfg = operation_config(
        path, allowed=required | optional, required=required, paths=["data", "boundaries", "pairs"]
    )
    manifest = metadata(cfg)
    boundaries = (
        read_table(cfg["boundaries"], required=BOUNDARY_FIELDS) if cfg.get("boundaries") else []
    )
    index = BoundaryIndex(boundaries)
    d_ref = number(cfg["d_ref"], "d_ref", minimum=0)
    d_min = number(cfg["d_min"], "d_min", minimum=0)
    if min(d_ref, d_min) <= 0:
        raise PaceError("d_ref and d_min must be positive")
    residuals, fitting_report = [], None
    if fit:
        rows = read_table(
            cfg["data"],
            required=[
                "chrom",
                "anchor0",
                "tss0",
                "raw_count",
                "count_to_contact",
                "contact_value",
                "sample_id",
                "bin_pair_id",
                "resolution",
                "measurement_status",
            ],
        )
        for row in rows:
            for key in MEASUREMENT_FIELDS - {"resolution"}:
                if row.get(key) != cfg[key]:
                    raise PaceError(f"Contact fitting {key} differs from configured measurement")
        rows = unique_count_pairs(rows, resolution=manifest["resolution"], boundaries=index)
        tests = cfg.get("test_chromosomes", [])
        if (
            not isinstance(tests, list)
            or len(set(tests)) != len(tests)
            or not set(tests) <= {r["chrom"] for r in rows}
        ):
            raise PaceError("test_chromosomes must be distinct chromosomes present in contact data")
        train = [r for r in rows if r["chrom"] not in tests]
        fitting_report = fit_contact_grid(
            train,
            gamma_grid=cfg["gamma_grid"],
            beta_grid=cfg["beta_grid"],
            d_ref=d_ref,
            d_min=d_min,
        )
        manifest.update(
            {k: fitting_report[k] for k in ("a", "gamma", "beta", "d_ref", "d_min", "kappa")}
        )
        manifest.update(
            training_sources=[file_hash(cfg["data"])],
            training_chromosomes=sorted({r["chrom"] for r in train}),
            test_chromosomes=sorted(tests),
            fitting_coordinate_policy="bin_centers",
        )
        for row in rows:
            if row["chrom"] in tests:
                expected = (
                    manifest["a"]
                    * (max(row["distance_bp"], d_min) / d_ref) ** -manifest["gamma"]
                    * math.exp(-manifest["beta"] * row["boundary_strength"])
                    / row["count_to_contact"]
                )
                residuals.append(
                    dict(
                        sample_id=row["sample_id"],
                        bin_pair_id=row["bin_pair_id"],
                        chrom=row["chrom"],
                        raw_count=row["raw_count"],
                        expected_raw_count=expected,
                        residual=row["raw_count"] - expected,
                    )
                )
    else:
        if cfg.get("test_chromosomes"):
            raise PaceError("test_chromosomes are only used by fit-hic")
        manifest.update({k: number(cfg[k], k, minimum=0) for k in ("a", "gamma", "beta")})
        if min(manifest["a"], manifest["gamma"]) <= 0:
            raise PaceError("Prior amplitude and gamma must be positive")
        manifest.update(
            d_ref=d_ref,
            d_min=d_min,
            parameter_source="user_supplied",
            prior_coordinate_policy="genomic_anchors",
        )
        if cfg.get("kappa") is not None:
            manifest["kappa"] = number(cfg["kappa"], "kappa", minimum=0)
            if manifest["kappa"] == 0:
                raise PaceError("kappa must be positive")
    if fit and cfg.get("kappa") is not None:
        raise PaceError("fit-hic estimates kappa; do not supply a manual value")
    if (
        manifest["beta"] > 0 or fit and any(float(b) > 0 for b in cfg["beta_grid"])
    ) and not cfg.get("boundaries"):
        raise PaceError("Positive beta requires an explicit boundary table")
    manifest["boundary_status"] = "candidate_unvalidated" if cfg.get("boundaries") else "not_used"
    predictions = []
    if cfg.get("pairs"):
        for row in read_table(cfg["pairs"], required=["chrom", "anchor0", "tss0"]):
            predictions.append(
                {
                    **row,
                    "prior_value": contact_prior(
                        row["chrom"],
                        integer(row["anchor0"], "anchor0"),
                        integer(row["tss0"], "tss0"),
                        manifest,
                        index,
                    ),
                }
            )
    with output_directory(out) as dest:
        if cfg.get("boundaries"):
            write_table(
                dest / "boundaries.tsv", boundaries, fields=None if boundaries else BOUNDARY_FIELDS
            )
            manifest.update(
                boundary_table="boundaries.tsv",
                boundary_sha256=file_hash(dest / "boundaries.tsv"),
                boundary_source_sha256=file_hash(cfg["boundaries"]),
            )
        write_json(dest / "manifest.json", manifest)
        if fitting_report:
            write_json(dest / "fit_report.json", fitting_report)
            write_table(
                dest / "held_out_residuals.tsv",
                residuals,
                fields=None
                if residuals
                else [
                    "sample_id",
                    "bin_pair_id",
                    "chrom",
                    "raw_count",
                    "expected_raw_count",
                    "residual",
                ],
            )
        if cfg.get("pairs"):
            write_table(dest / "prior_contacts.tsv", predictions)
    return {
        "output": str(out),
        "model_id": manifest["model_id"],
        "gamma": manifest["gamma"],
        "beta": manifest["beta"],
    }


def labels_command(path, out):
    from .weak_labels import calibrate

    required = {
        "run_config",
        "labels",
        "gamma_grid",
        "beta_grid",
        "eta_grid",
        "test_chromosomes",
        "aggregation",
    }
    cfg = operation_config(
        path,
        allowed=required | {"abc_gamma", "species", "assembly", "context_id"},
        required=required,
        paths=["run_config", "labels"],
    )
    run_cfg = load_config(cfg["run_config"])
    # The run fixes the context; an explicitly stated context must agree with it.
    if any(
        cfg.get(k) is not None and cfg[k] != run_cfg["context"][k]
        for k in ("species", "assembly", "context_id")
    ):
        raise PaceError("eQTL species/assembly/context differs from the run")
    variants = read_table(cfg["labels"], required=["variant_id", "chrom", "pos0", "gene_id", "pip"])
    model, labels = calibrate(
        run_cfg,
        variants,
        **{
            k: cfg[k]
            for k in ("gamma_grid", "beta_grid", "eta_grid", "test_chromosomes", "aggregation")
        },
        abc_gamma=cfg.get("abc_gamma", 1.0242386),
    )
    model["labels_sha256"] = file_hash(cfg["labels"])
    with output_directory(out) as dest:
        write_json(dest / "weak_model.json", model)
        write_table(dest / "weak_labels.tsv", labels)
    return {
        "output": str(out),
        "selected": model["selected"],
        "status": "weak_fitted",
        "functional_validation": False,
    }


def contact_command(command, path, out):
    if command == "boundaries":
        return boundaries_command(path, out)
    if command in ("prior", "fit-hic"):
        return prior_command(path, out, fit=command == "fit-hic")
    if command == "fit-labels":
        return labels_command(path, out)
    from .pipeline import run

    cfg = load_config(path)
    if cfg["contact"]["mode"] != "shrinkage" or cfg["contact"]["reliability"] != "per_pair":
        raise PaceError("fuse needs contact.mode=shrinkage and reliability=per_pair")
    result = run(cfg, out)
    return {
        "output": str(out),
        "n_candidates": len(result["scores"]),
        "eta": result["eta_calibration"]["eta"],
    }
