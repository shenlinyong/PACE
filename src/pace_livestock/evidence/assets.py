"""Model asset contracts, scope checking and capability reports."""

from pathlib import Path

from ..errors import PaceError
from ..provenance import file_hash, read_json

TASKS = ("locus_prediction", "individual_delta", "link_prediction")


def load_asset(path, cfg: dict, *, kind: str) -> dict:
    root = Path(path)
    manifest_path = root / "manifest.json" if root.is_dir() else root
    if not manifest_path.exists():
        raise PaceError(f"Missing {kind} asset: {manifest_path}")
    m = read_json(manifest_path)
    required = {
        "model_id",
        "kind",
        "is_synthetic",
        "species",
        "assembly",
        "context_id",
        "target_level",
        "training_sources",
        "calibration_sources",
        "test_sources",
        "validation",
    }
    missing = required - set(m)
    if missing:
        raise PaceError(f"{kind} manifest lacks {sorted(missing)}")
    if m["kind"] != kind or type(m["is_synthetic"]) is not bool:
        raise PaceError(f"Invalid kind or is_synthetic in {manifest_path}")
    if cfg["execution_profile"] != "demonstration" and m["is_synthetic"]:
        raise PaceError(
            f"Synthetic {kind} weights cannot be used in {cfg['execution_profile']} profile"
        )
    if cfg["execution_profile"] == "demonstration" and not m["is_synthetic"]:
        raise PaceError("demonstration requires explicitly synthetic assets")
    for field, expected in {**cfg["context"], "target_level": cfg["target_level"]}.items():
        if m[field] != expected:
            raise PaceError(f"{kind} {field} mismatch: {m[field]!r} versus {expected!r}")
    if not isinstance(m["validation"], dict):
        raise PaceError("Asset validation must be a task/report mapping")
    for task, report in m["validation"].items():
        if task not in TASKS or not isinstance(report, dict):
            raise PaceError(f"Invalid validation task {task}")
        if report.get("validated"):
            p = manifest_path.parent / report.get("report_path", "")
            if not p.is_file() or file_hash(p) != report.get("sha256"):
                raise PaceError(
                    f"{task}: validated claim requires an existing report and matching sha256"
                )
            if report.get("context_id") != cfg["context"]["context_id"]:
                raise PaceError(f"{task}: validation report context mismatch")
    if cfg["execution_profile"] == "validated":
        required_task = "locus_prediction" if kind == "sequence" else "link_prediction"
        if not m["validation"].get(required_task, {}).get("validated", False):
            raise PaceError(f"validated profile requires {required_task} evidence for {kind}")
    if kind == "sequence":
        for key in (
            "assays",
            "input_length",
            "output_window",
            "output_type",
            "signal_unit",
            "normalization_id",
            "adapter",
            "weights_sha256",
        ):
            if key not in m:
                raise PaceError(f"Sequence manifest requires {key}")
        if m["output_type"] != "normalized_signal":
            raise PaceError("Sequence peak probabilities/logits are not quantitative activity")
        if m["output_window"] != cfg["catalog"]["width_bp"]:
            raise PaceError("Sequence output window does not match scoring units")
        if not set(cfg["activity"]["panel"]) <= set(m["assays"]):
            raise PaceError("Sequence model does not contain all requested assay heads")
        if m["adapter"] in ("cnn", "fixed_linear"):
            weights = manifest_path.parent / m.get("weights_file", "weights.safetensors")
            if not weights.is_file() or file_hash(weights) != m["weights_sha256"]:
                raise PaceError("Missing sequence weights or SHA256 mismatch")
        elif m["adapter"] != "external_table":
            raise PaceError(f"Unknown sequence adapter {m['adapter']}")
    m["asset_directory"] = str(manifest_path.parent.resolve())
    m["manifest_sha256"] = file_hash(manifest_path)
    return m


def capabilities(cfg: dict) -> dict:
    assets, blocks = {}, []
    for kind, path in (
        ("sequence", cfg["sequence"]["model_path"]),
        ("contact_prior", cfg["contact"]["prior_path"]),
        ("fusion", cfg["fusion"]["calibrator_path"]),
    ):
        if not path:
            assets[kind] = {"weights_status": "absent"}
            if (kind == "sequence" and cfg["regime"] == "genome_only") or (
                kind == "contact_prior" and cfg["contact"]["mode"] in ("prior_only", "shrinkage")
            ):
                blocks.append(f"{kind}: required asset absent")
            continue
        try:
            m = load_asset(path, cfg, kind=kind)
            assets[kind] = {
                "weights_status": "synthetic" if m["is_synthetic"] else "real",
                "model_id": m["model_id"],
                "task_validation_status": m["validation"],
            }
        except (PaceError, OSError, ValueError) as exc:
            assets[kind] = {"weights_status": "invalid", "reason": str(exc)}
            blocks.append(str(exc))
    return {
        "implementation_status": "implemented",
        "execution_profile": cfg["execution_profile"],
        "scope": cfg["context"],
        "regime": cfg["regime"],
        "assets": assets,
        "task_validation_status": {
            task: "asset_report_available"
            if any(
                a.get("task_validation_status", {}).get(task, {}).get("validated", False)
                for a in assets.values()
            )
            else "not_assessed"
            for task in TASKS
        },
        "blocking_reasons": blocks,
        "ready": not blocks,
        "research_extensions": {
            k: "not_implemented"
            for k in [
                "genome_wide_discovery",
                "copy_resolved_support",
                "complex_cnv_reconstruction",
                "sequence_contact_training",
                "poisson_gamma_shrinkage",
                "posterior_sampling",
            ]
        },
    }
