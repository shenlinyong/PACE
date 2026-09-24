"""Validated YAML settings and configuration-relative paths."""

from __future__ import annotations

import copy
from pathlib import Path

import yaml

from .errors import PaceError
from .io.tables import integer, number

DEFAULTS = {
    "schema_version": "pace-1",
    "run_id": "pace",
    "regime": "measured",
    "execution_profile": "research",
    "estimand": "bulk_proxy",
    "target_level": "individual",
    "context": {"species": None, "assembly": None, "context_id": None},
    "inputs": {
        k: None
        for k in (
            "units",
            "region_membership",
            "promoters",
            "candidates",
            "samples",
            "observed_activity",
            "observed_contacts",
            "resolved_activity",
            "resolved_contacts",
            "features",
            "methylation",
            "expression",
            "labels",
            "evidence",
            "sources",
            "support_bounds",
        )
    },
    "catalog": {
        "profile": "canonical_grid",
        "width_bp": 500,
        "offset_bp": 0,
        "include_promoter_units": True,
        "chrom_sizes_path": None,
        "candidate_radius_bp": 5_000_000,
    },
    "activity": {
        "panel": ["ATAC", "H3K27ac"],
        "combine": "geometric_equal",
        "missing_policy": "unresolved",
        "minimum_callable_fraction": 0.0,
        "replicate_aggregation": "equal_donor_mean",
        "pseudocounts": {},
    },
    "contact": {
        "mode": "observed",
        "scale": "depth_normalized_contact",
        "prior_path": None,
        "prior_preset": None,
        "near_diagonal_policy": "prior_or_neighbor",
        "near_diagonal_bp": 0,
        "allow_prior_fallback": False,
        "allow_cross_context_prior": False,
        "reliability": None,
        "kappa": "auto",
        "reliability_source": None,
        "resolution": None,
        "normalization_id": None,
        "balancing": None,
        "window_id": None,
        "pseudocount": "auto",
        "pseudocount_distance_bp": 5000,
        "pseudocount_strength": 1.0,
    },
    "scoring": {"partial_policy": "withhold"},
    "promoters": {
        "weights": "provided",
        "minimum_weight": 0.0,
        "missing_policy": "strict",
        "minimum_retained_weight": 0.9,
    },
    "allocation": {
        "eta": "auto",
        "missing_policy": "fixed_gene_set",
        "labels_path": None,
        "calibrator_path": None,
        "weak_model_path": None,
        "minimum_genes": 3,
        "minimum_groups": 3,
        "validation_folds": 5,
        "minimum_positive_fraction": 0.8,
    },
    "multiomics": {"mode": "annotate", "model_path": None},
    "methylation": {
        "minimum_coverage": 1,
        "promoter_upstream_bp": 2000,
        "promoter_downstream_bp": 500,
        "reference_cpg_path": None,
    },
    "comparison": {
        "full_delta_requires_complete": True,
        "allow_conditional_intersection": True,
        "minimum_common_units": 2,
    },
    "output": {"format": "tsv_gz", "retain_all_candidates": True},
    "seed": 17,
}


class StrictLoader(yaml.SafeLoader):
    """Reject duplicate keys rather than silently using the last spelling."""


def _mapping(loader, node, deep=False):
    mapping = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=deep)
        if key in mapping:
            raise PaceError(f"Duplicate YAML key: {key}")
        mapping[key] = loader.construct_object(value_node, deep=deep)
    return mapping


StrictLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _mapping)


def load_yaml(path: str | Path) -> dict:
    result = yaml.load(Path(path).read_text(encoding="utf-8"), Loader=StrictLoader)
    if not isinstance(result, dict):
        raise PaceError(f"{path}: expected a YAML mapping")
    return result


def strict_keys(data: dict, allowed, name: str) -> None:
    if not isinstance(data, dict):
        raise PaceError(f"{name}: expected a mapping")
    extra = set(data) - set(allowed)
    if extra:
        raise PaceError(f"{name}: unknown keys {sorted(extra)}")


def load_config(path: str | Path | None = None, *, overrides: dict | None = None) -> dict:
    if path is not None and Path(path).is_dir():
        # A previous result folder stores its complete settings.
        path = Path(path) / "resolved_config.yaml"
        if not path.is_file():
            raise PaceError(f"{path.parent} is not a PACE result folder (no resolved_config.yaml)")
    base = Path(path).resolve().parent if path is not None else Path.cwd()
    user = load_yaml(path) if path is not None else {}
    for key, value in (overrides or {}).items():
        if isinstance(value, dict) and isinstance(user.get(key, {}), dict):
            user.setdefault(key, {}).update(value)
        else:
            user[key] = value
    strict_keys(user, DEFAULTS, "config")
    cfg = copy.deepcopy(DEFAULTS)
    for key, value in user.items():
        if isinstance(DEFAULTS[key], dict):
            strict_keys(value, DEFAULTS[key], key)
            cfg[key].update(value)
        else:
            cfg[key] = value
    for key in ("species", "assembly", "context_id"):
        if not isinstance(cfg["context"][key], str) or not cfg["context"][key]:
            flag = {"context_id": "tissue"}.get(key, key)
            raise PaceError(f"--{flag} is required (context.{key})")
    enums = {
        "schema_version": {"pace-1"},
        "regime": {"measured"},
        "execution_profile": {"demonstration", "research", "validated"},
        "estimand": {"bulk_proxy"},
        "target_level": {"individual", "population_mean"},
    }
    for key, values in enums.items():
        if cfg[key] not in values:
            raise PaceError(f"{key} must be one of {sorted(values)}")
    panels = [{"ATAC"}, {"DNase"}, {"H3K27ac"}, {"ATAC", "H3K27ac"}, {"DNase", "H3K27ac"}]
    panel = cfg["activity"]["panel"]
    if not isinstance(panel, list) or len(panel) != len(set(panel)) or set(panel) not in panels:
        raise PaceError("activity.panel must be ATAC, DNase, H3K27ac, or accessibility + H3K27ac")
    choices = [
        ("catalog", "profile", {"canonical_grid", "provided_regions"}),
        ("activity", "combine", {"geometric_equal"}),
        ("activity", "missing_policy", {"unresolved"}),
        ("activity", "replicate_aggregation", {"equal_donor_mean"}),
        ("allocation", "missing_policy", {"fixed_gene_set"}),
        ("contact", "mode", {"observed", "prior_only", "shrinkage"}),
        (
            "contact",
            "near_diagonal_policy",
            {"prior_or_neighbor", "prior_or_unresolved", "unresolved"},
        ),
        ("contact", "pseudocount", {"auto", "none", "powerlaw"}),
        ("scoring", "partial_policy", {"withhold", "conditional"}),
        ("promoters", "weights", {"provided", "equal"}),
        ("promoters", "missing_policy", {"strict", "drop_missing"}),
        ("multiomics", "mode", {"annotate", "ml"}),
        ("output", "format", {"tsv_gz"}),
    ]
    for section, key, values in choices:
        if cfg[section][key] not in values:
            raise PaceError(f"{section}.{key}: expected one of {sorted(values)}")
    pseudocounts = cfg["activity"]["pseudocounts"]
    strict_keys(pseudocounts, panel, "activity.pseudocounts")
    cfg["activity"]["pseudocounts"] = {
        assay: number(value, f"activity.pseudocounts.{assay}", minimum=0)
        for assay, value in pseudocounts.items()
    }
    for key in ("minimum_weight", "minimum_retained_weight"):
        cfg["promoters"][key] = number(
            cfg["promoters"][key], f"promoters.{key}", minimum=0, maximum=1
        )
    if cfg["promoters"]["minimum_retained_weight"] == 0:
        raise PaceError("promoters.minimum_retained_weight must be greater than zero")
    allocation = cfg["allocation"]
    if allocation["eta"] != "auto":
        allocation["eta"] = number(allocation["eta"], "allocation.eta", minimum=0, maximum=1)
        if allocation["labels_path"] or allocation["calibrator_path"]:
            raise PaceError("Use allocation.eta=auto with functional labels or a frozen calibrator")
    if allocation["weak_model_path"]:
        if (
            allocation["eta"] != "auto"
            or allocation["labels_path"]
            or allocation["calibrator_path"]
        ):
            raise PaceError("weak_model_path requires eta=auto and no functional calibrator/labels")
        if cfg["execution_profile"] == "validated":
            raise PaceError("eQTL weak calibration is not functional validation")
    if allocation["labels_path"] and allocation["calibrator_path"]:
        raise PaceError("Provide eta labels OR a frozen eta calibrator, not both")
    allocation["minimum_genes"] = integer(
        allocation["minimum_genes"], "allocation.minimum_genes", minimum=2
    )
    for key in ("minimum_groups", "validation_folds"):
        allocation[key] = integer(allocation[key], f"allocation.{key}", minimum=3)
    allocation["minimum_positive_fraction"] = number(
        allocation["minimum_positive_fraction"],
        "allocation.minimum_positive_fraction",
        minimum=0.5,
        maximum=1,
    )
    if allocation["minimum_positive_fraction"] <= 0.5:
        raise PaceError("allocation.minimum_positive_fraction must be greater than 0.5")
    cfg["catalog"]["width_bp"] = integer(cfg["catalog"]["width_bp"], "catalog.width_bp", minimum=1)
    cfg["catalog"]["offset_bp"] = integer(cfg["catalog"]["offset_bp"], "catalog.offset_bp")
    if cfg["catalog"]["offset_bp"] >= cfg["catalog"]["width_bp"]:
        raise PaceError("catalog.offset_bp must be smaller than catalog.width_bp")
    cfg["catalog"]["candidate_radius_bp"] = integer(
        cfg["catalog"]["candidate_radius_bp"], "catalog.candidate_radius_bp", minimum=1
    )
    cfg["seed"] = integer(cfg["seed"], "seed")
    cfg["contact"]["near_diagonal_bp"] = integer(
        cfg["contact"]["near_diagonal_bp"], "contact.near_diagonal_bp"
    )
    cfg["contact"]["pseudocount_distance_bp"] = integer(
        cfg["contact"]["pseudocount_distance_bp"], "contact.pseudocount_distance_bp", minimum=1
    )
    cfg["contact"]["pseudocount_strength"] = number(
        cfg["contact"]["pseudocount_strength"], "contact.pseudocount_strength", minimum=0
    )
    preset = cfg["contact"]["prior_preset"]
    if preset not in (None, "abc_human"):
        raise PaceError("contact.prior_preset must be null or abc_human")
    if preset and (cfg["contact"]["prior_path"] or cfg["contact"]["mode"] != "prior_only"):
        raise PaceError(
            "abc_human is an explicit prior_only baseline; do not mix its arbitrary scale with measured contacts"
        )
    if preset and cfg["execution_profile"] != "research":
        raise PaceError("abc_human is an unvalidated transferred prior for research only")
    if preset:
        cfg["contact"].update(
            scale="relative_distance_contact",
            resolution=5000,
            normalization_id="abc_human_relative_shape",
            balancing="not_applicable",
            window_id="bin_pair",
        )
    if cfg["contact"]["resolution"] is not None:
        cfg["contact"]["resolution"] = integer(
            cfg["contact"]["resolution"], "contact.resolution", minimum=1
        )
    for key in ("scale", "normalization_id", "balancing", "window_id", "reliability_source"):
        value = cfg["contact"][key]
        if value is not None and (not isinstance(value, str) or not value.strip()):
            raise PaceError(f"contact.{key} must be a nonempty string or null")
    if cfg["contact"]["scale"] is None:
        raise PaceError("contact.scale must name the contact measurement scale")
    if cfg["contact"]["kappa"] != "auto":
        cfg["contact"]["kappa"] = number(cfg["contact"]["kappa"], "contact.kappa", minimum=0)
        if cfg["contact"]["kappa"] == 0:
            raise PaceError("contact.kappa must be positive")
    if cfg["contact"]["reliability"] == "per_pair":
        if cfg["contact"]["mode"] != "shrinkage":
            raise PaceError("per_pair reliability requires shrinkage mode")
        cfg["contact"]["reliability_source"] = "gamma_poisson"
    elif cfg["contact"]["reliability"] is not None:
        cfg["contact"]["reliability"] = number(
            cfg["contact"]["reliability"], "contact.reliability", minimum=0, maximum=1
        )
    cfg["activity"]["minimum_callable_fraction"] = number(
        cfg["activity"]["minimum_callable_fraction"],
        "minimum_callable_fraction",
        minimum=0,
        maximum=1,
    )
    for key in ("minimum_coverage", "promoter_upstream_bp", "promoter_downstream_bp"):
        cfg["methylation"][key] = integer(
            cfg["methylation"][key],
            f"methylation.{key}",
            minimum=1 if key == "minimum_coverage" else 0,
        )
    if (
        not cfg["methylation"]["promoter_upstream_bp"]
        + cfg["methylation"]["promoter_downstream_bp"]
    ):
        raise PaceError("Methylation promoter window must have a positive width")
    cfg["comparison"]["minimum_common_units"] = integer(
        cfg["comparison"]["minimum_common_units"], "comparison.minimum_common_units", minimum=2
    )
    for section, key in [
        ("catalog", "include_promoter_units"),
        ("contact", "allow_prior_fallback"),
        ("contact", "allow_cross_context_prior"),
        ("output", "retain_all_candidates"),
        ("comparison", "full_delta_requires_complete"),
        ("comparison", "allow_conditional_intersection"),
    ]:
        if type(cfg[section][key]) is not bool:
            raise PaceError(f"{section}.{key} must be a YAML boolean")
    if (
        not cfg["output"]["retain_all_candidates"]
        or not cfg["comparison"]["full_delta_requires_complete"]
    ):
        raise PaceError("Dropping candidates or reporting incomplete full deltas is unsupported")
    if cfg["contact"]["scale"] in ("oe", "O/E", "log_oe", "pvalue", "correlation"):
        raise PaceError(
            "Contact scale must retain distance background; convert O/E using matching expected contacts"
        )
    for section in (
        "inputs",
        "contact",
        "multiomics",
        "allocation",
        "catalog",
        "methylation",
    ):
        for key, value in cfg[section].items():
            if value is not None and (section == "inputs" or key.endswith("_path")):
                if not isinstance(value, str):
                    raise PaceError(f"{section}.{key}: expected a path string")
                cfg[section][key] = str((base / value).resolve())
    return cfg


# Command-line spelling of operation keys whose option name differs from the key.
OPTION_NAMES = {
    "context_id": "tissue",
    "run_config": "run",
    "data": "contacts",
    "bed": "peaks",
    "track": "bigwig",
    "contact": "hic",
    "region_membership": "membership",
    "is_synthetic": "synthetic",
}


def operation_base(source) -> Path:
    """Directory that relative paths of an operation resolve against."""
    return Path.cwd() if isinstance(source, dict) else Path(source).resolve().parent


def operation_config(source, *, allowed, required=(), paths=()) -> dict:
    """Validate operation settings given as command-line options or a legacy YAML file.

    A mapping (built from command-line options) resolves relative paths against the
    working directory; a YAML path resolves them against the YAML file's folder.
    Unset (None) command-line options are dropped so defaults apply.
    """
    if isinstance(source, dict):
        cfg = {k: v for k, v in source.items() if v is not None}
    else:
        cfg = load_yaml(source)
    strict_keys(cfg, allowed, "options")
    for key in sorted(required):
        if cfg.get(key) is None:
            flag = OPTION_NAMES.get(key, key.replace("_", "-"))
            raise PaceError(f"Required option missing: --{flag}")
    base = operation_base(source)
    for key in paths:
        if cfg.get(key):
            cfg[key] = str((base / cfg[key]).resolve())
    return cfg
