"""Direct file options for the public PACE executable; YAML remains optional."""

from __future__ import annotations

import argparse
from pathlib import Path

from .config import DEFAULTS, load_config, load_yaml
from .errors import PaceError

MODES = {
    "measured": "measured",
}
PATH_OPTIONS = {
    "units": ("inputs", "units"),
    "promoters": ("inputs", "promoters"),
    "candidates": ("inputs", "candidates"),
    "samples": ("inputs", "samples"),
    "sources": ("inputs", "sources"),
    "evidence": ("inputs", "evidence"),
    "activity": ("inputs", "observed_activity"),
    "contacts": ("inputs", "observed_contacts"),
    "resolved_activity": ("inputs", "resolved_activity"),
    "resolved_contacts": ("inputs", "resolved_contacts"),
    "features": ("inputs", "features"),
    "expression": ("inputs", "expression"),
    "methylation": ("inputs", "methylation"),
    "contact_prior": ("contact", "prior_path"),
    "ml_model": ("multiomics", "model_path"),
    "eta_labels": ("allocation", "labels_path"),
    "eta_model": ("allocation", "calibrator_path"),
    "chrom_sizes": ("catalog", "chrom_sizes_path"),
    "reference_cpg": ("methylation", "reference_cpg_path"),
    "support_bounds": ("inputs", "support_bounds"),
}
VALUE_OPTIONS = {
    "species": ("context", "species"),
    "assembly": ("context", "assembly"),
    "tissue": ("context", "context_id"),
    "eta": ("allocation", "eta"),
    "eta_min_genes": ("allocation", "minimum_genes"),
    "eta_min_groups": ("allocation", "minimum_groups"),
    "eta_validation_folds": ("allocation", "validation_folds"),
    "panel": ("activity", "panel"),
    "contact_mode": ("contact", "mode"),
    "contact_scale": ("contact", "scale"),
    "contact_resolution": ("contact", "resolution"),
    "contact_normalization": ("contact", "normalization_id"),
    "contact_balancing": ("contact", "balancing"),
    "contact_window": ("contact", "window_id"),
    "contact_reliability": ("contact", "reliability"),
    "reliability_source": ("contact", "reliability_source"),
    "allow_prior_fallback": ("contact", "allow_prior_fallback"),
    "allow_cross_context_prior": ("contact", "allow_cross_context_prior"),
    "minimum_tss_weight": ("promoters", "minimum_weight"),
    "missing_tss_policy": ("promoters", "missing_policy"),
    "minimum_retained_tss_weight": ("promoters", "minimum_retained_weight"),
    "prior_preset": ("contact", "prior_preset"),
    "pseudocount": ("contact", "pseudocount"),
    "partial_policy": ("scoring", "partial_policy"),
    "catalog_profile": ("catalog", "profile"),
    "unit_width": ("catalog", "width_bp"),
    "grid_offset": ("catalog", "offset_bp"),
    "include_promoters": ("catalog", "include_promoter_units"),
    "candidate_radius": ("catalog", "candidate_radius_bp"),
    "methylation_min_coverage": ("methylation", "minimum_coverage"),
    "promoter_upstream": ("methylation", "promoter_upstream_bp"),
    "promoter_downstream": ("methylation", "promoter_downstream_bp"),
}


def add_run_options(
    parser: argparse.ArgumentParser, *, output: bool = True, help_mode: str | None = None
) -> None:
    parser.add_argument(
        "--config", help="Optional YAML; relative paths are relative to the YAML file"
    )
    parser.add_argument(
        "--mode",
        "--regime",
        choices=list(MODES),
        help=(
            argparse.SUPPRESS
            if help_mode
            else "Measured activity (the default and only supported mode)"
        ),
    )
    if output:
        parser.add_argument(
            "-o",
            "--out",
            required=True,
            help="New output directory; --force retains a backup of recognized PACE outputs",
        )
    context = parser.add_argument_group("sample and biological context")
    for flag in ("species", "assembly", "tissue", "run-id"):
        context.add_argument("--" + flag)
    context.add_argument("--target-level", choices=["individual", "population_mean"])
    context.add_argument(
        "--profile",
        choices=["research", "validated", "demonstration"],
        help="Default: research; synthetic assets require demonstration",
    )
    context.add_argument("--seed", type=int)
    files = parser.add_argument_group("canonical input tables and genomic assets")
    files.add_argument(
        "--catalog-dir",
        help="Read available <table>.tsv[.gz] files from this directory; explicit file options take precedence",
    )
    for flag, (section, key) in PATH_OPTIONS.items():
        aliases = ["--" + flag.replace("_", "-")]
        if flag == "activity":
            aliases.append("--observed-activity")
        if flag == "contacts":
            aliases.append("--observed-contacts")
        descriptions = {
            "activity": "Measured activity table (ATAC-seq, DNase-seq and/or H3K27ac; one row per element/sample/assay)",
            "contacts": "Measured promoter contact table (Hi-C/Prom-Hi-C; one row per element/promoter/sample)",
            "samples": "Sample metadata table; list every biological/technical replicate with its sample_id",
            "expression": "Optional RNA-seq gene-expression table used as an annotation",
            "methylation": "Optional WGBS/RRBS methylation table used as an annotation",
        }
        files.add_argument(
            *aliases,
            help=(
                descriptions.get(
                    flag, f"{section}.{key}; paths are relative to the working directory"
                )
            ),
        )
    model = parser.add_argument_group("scoring and calibration")
    model.add_argument("--eta", help="auto (default; falls back to 0) or a fixed number in [0,1]")
    model.add_argument(
        "--eta-min-genes",
        type=int,
        help="Minimum informative genes for automatic fitting (default: 3)",
    )
    model.add_argument("--panel", nargs="+", choices=["ATAC", "DNase", "H3K27ac"])
    model.add_argument("--eta-min-groups", type=int)
    model.add_argument("--eta-validation-folds", type=int)
    model.add_argument("--contact-mode", choices=["observed", "prior_only", "shrinkage"])
    model.add_argument(
        "--prior-preset",
        choices=["abc_human"],
        help="Explicit unvalidated human contact-shape baseline; requires prior_only",
    )
    model.add_argument("--pseudocount", choices=["auto", "none", "powerlaw"])
    model.add_argument("--partial-policy", choices=["withhold", "conditional"])
    model.add_argument("--minimum-tss-weight", type=float)
    model.add_argument("--missing-tss-policy", choices=["strict", "drop_missing"])
    model.add_argument("--minimum-retained-tss-weight", type=float)
    model.add_argument(
        "--allow-cross-context-prior",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Use a same-species, same-assembly prior from another tissue; records the transfer",
    )
    model.add_argument("--contact-scale")
    model.add_argument(
        "--contact-resolution",
        type=int,
        help="Common contact resolution in bp; inputs and prior must match",
    )
    model.add_argument("--contact-normalization")
    model.add_argument("--contact-balancing")
    model.add_argument("--contact-window")
    model.add_argument("--contact-reliability", type=float)
    model.add_argument("--reliability-source")
    model.add_argument(
        "--allow-prior-fallback", action=argparse.BooleanOptionalAction, default=None
    )
    catalog = parser.add_argument_group("candidate universe")
    catalog.add_argument("--catalog-profile", choices=["canonical_grid", "provided_regions"])
    catalog.add_argument("--unit-width", type=int)
    catalog.add_argument("--grid-offset", type=int)
    catalog.add_argument(
        "--candidate-radius",
        type=int,
        help="Declared cis candidate radius in bp (default: 5000000)",
    )
    catalog.add_argument("--include-promoters", action=argparse.BooleanOptionalAction, default=None)
    omics = parser.add_argument_group("methylation annotation")
    omics.add_argument("--methylation-min-coverage", type=int)
    omics.add_argument(
        "--promoter-upstream",
        type=int,
        help="Upstream methylation window in transcription direction",
    )
    omics.add_argument(
        "--promoter-downstream",
        type=int,
        help="Downstream methylation window in transcription direction",
    )


def config_from_args(args: argparse.Namespace) -> dict:
    overrides = {}

    def put(section, key, value):
        overrides.setdefault(section, {})[key] = value

    if args.catalog_dir:
        root = Path(args.catalog_dir).resolve()
        if not root.is_dir():
            raise PaceError(f"--catalog-dir is not a directory: {root}")
        metadata = root / "run_catalog_config.yaml"
        if metadata.is_file():
            for key, value in load_yaml(metadata).get("catalog", {}).items():
                if key.endswith("_path") and value:
                    value = str((root / value).resolve())
                put("catalog", key, value)
        for name in DEFAULTS["inputs"]:
            plain, zipped = root / f"{name}.tsv", root / f"{name}.tsv.gz"
            if plain.exists() and zipped.exists():
                raise PaceError(f"Ambiguous catalog table: both {plain.name} and {zipped.name}")
            if plain.exists() or zipped.exists():
                put("inputs", name, str(plain if plain.exists() else zipped))
    for option, (section, key) in PATH_OPTIONS.items():
        value = getattr(args, option)
        if value is not None:
            put(section, key, str(Path(value).resolve()))
    for option, (section, key) in VALUE_OPTIONS.items():
        value = getattr(args, option)
        if value is not None:
            put(section, key, value)
    for option, key in [
        ("run_id", "run_id"),
        ("target_level", "target_level"),
        ("profile", "execution_profile"),
        ("seed", "seed"),
    ]:
        if getattr(args, option) is not None:
            overrides[key] = getattr(args, option)
    if args.mode:
        overrides["regime"] = MODES[args.mode]
    if args.ml_model:
        put("multiomics", "mode", "ml")
    if args.eta_labels or args.eta_model:
        if args.eta is None:
            put("allocation", "eta", "auto")
        # A CLI calibration source explicitly replaces the other source from YAML.
        if args.eta_labels and not args.eta_model:
            put("allocation", "calibrator_path", None)
        if args.eta_model and not args.eta_labels:
            put("allocation", "labels_path", None)
    if args.eta is not None and args.eta != "auto" and not (args.eta_labels or args.eta_model):
        put("allocation", "labels_path", None)
        put("allocation", "calibrator_path", None)
    return load_config(args.config, overrides=overrides)
