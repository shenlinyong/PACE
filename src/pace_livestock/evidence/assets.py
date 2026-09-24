"""Model asset contracts, scope checking and capability reports."""

from pathlib import Path

from ..errors import PaceError
from ..provenance import digest, file_hash, read_json

TASKS = ("link_prediction",)


def builtin_contact_prior(cfg):
    """Explicit human-derived shape on a relative scale, never a calibrated animal asset."""
    if cfg["contact"].get("prior_preset") != "abc_human":
        return None
    if cfg["contact"]["mode"] != "prior_only" or cfg["execution_profile"] != "research":
        raise PaceError("abc_human requires research profile and prior_only contact mode")
    m = {
        "model_id": "abc_human_shape_1.024238616787792",
        "kind": "contact_prior",
        "is_synthetic": False,
        **cfg["context"],
        "target_level": cfg["target_level"],
        "prior_source": "abc_human_default",
        "source_species": "Homo sapiens",
        "transfer_status": "unvalidated_for_target_context",
        "reference": "https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction/blob/main/config/config.yaml",
        "training_sources": ["ABC published human reference gamma"],
        "calibration_sources": [],
        "test_sources": [],
        "validation": {},
        "a": 1.0,
        "gamma": 1.024238616787792,
        "d_ref": 5000.0,
        "d_min": 5000.0,
        "scale": "relative_distance_contact",
        "resolution": 5000,
        "normalization_id": "abc_human_relative_shape",
        "balancing": "not_applicable",
        "window_id": "bin_pair",
    }
    m["manifest_sha256"] = digest(m)
    return m


def load_asset(path, cfg: dict, *, kind: str) -> dict:
    if kind != "contact_prior":
        raise PaceError("Only contact-prior assets are supported by this activity model")
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
    transferred = m["context_id"] != cfg["context"]["context_id"]
    for field, expected in {**cfg["context"], "target_level": cfg["target_level"]}.items():
        if m[field] != expected:
            if (
                field == "context_id"
                and cfg["contact"].get("allow_cross_context_prior", False)
                and cfg["execution_profile"] != "validated"
            ):
                continue
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
            if report.get("context_id") != m["context_id"]:
                raise PaceError(f"{task}: validation report context mismatch")
    if cfg["execution_profile"] == "validated":
        required_task = "link_prediction"
        if not m["validation"].get(required_task, {}).get("validated", False):
            raise PaceError(f"validated profile requires {required_task} evidence for {kind}")
    m["asset_directory"] = str(manifest_path.parent.resolve())
    if "beta" in m or "boundary_table" in m:
        from ..boundary_prior import load_boundaries
        from ..io.tables import number

        number(m.get("beta", 0), "contact prior beta", minimum=0)
        load_boundaries(m)
        if "kappa" in m and number(m["kappa"], "prior kappa", minimum=0) == 0:
            raise PaceError("Prior kappa must be positive")
    m["manifest_sha256"] = file_hash(manifest_path)
    if transferred:
        m["source_context_id"] = m["context_id"]
        m["target_context_id"] = cfg["context"]["context_id"]
        m["transfer_status"] = "cross_context_unvalidated"
        m["source_validation"] = m["validation"]
        m["validation"] = {}
    return m


def capabilities(cfg: dict) -> dict:
    assets, blocks = {}, []
    for kind, path in (("contact_prior", cfg["contact"]["prior_path"]),):
        if cfg["contact"].get("prior_preset"):
            m = builtin_contact_prior(cfg)
            assets[kind] = {
                "weights_status": "transferred_unvalidated",
                "model_id": m["model_id"],
                "task_validation_status": {},
            }
            continue
        if not path:
            assets[kind] = {"weights_status": "absent"}
            if cfg["contact"]["mode"] in ("prior_only", "shrinkage"):
                blocks.append(f"{kind}: required asset absent")
            continue
        try:
            m = load_asset(path, cfg, kind=kind)
            assets[kind] = {
                "weights_status": "transferred_unvalidated"
                if m.get("transfer_status")
                else "synthetic"
                if m["is_synthetic"]
                else "real",
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
    }
