"""Measured-only scope, experimental imports and honest evaluation reports."""

import copy
import json
import math
from pathlib import Path

import pytest
from test_canonical_learning import learning_fixture

from pace_livestock.cli import main
from pace_livestock.config import load_config
from pace_livestock.demo import create_example
from pace_livestock.errors import PaceError
from pace_livestock.evidence.resolve import resolve_contacts
from pace_livestock.io.tables import read_table, write_table
from pace_livestock.learning.model import train_classifier
from pace_livestock.pipeline import compute, run
from pace_livestock.schemas import load_tables, validate_tables


@pytest.mark.parametrize("mode", ["hybrid", "genome", "genome_only"])
def test_retired_modes_fail_before_writing(tmp_path, mode):
    with pytest.raises(SystemExit) as error:
        main([mode, "--out", str(tmp_path / "out")])
    assert error.value.code == 2
    assert not (tmp_path / "out").exists()
    path = create_example(tmp_path / "input")
    with pytest.raises(PaceError, match="regime"):
        load_config(path, overrides={"regime": mode})


@pytest.mark.parametrize(
    "command",
    [
        "train-sequence",
        "predict-sequence",
        "prepare-genome",
        "fit-fusion",
        "variant-effects",
    ],
)
def test_retired_commands_are_unavailable(command):
    with pytest.raises(SystemExit) as error:
        main([command, "--help"])
    assert error.value.code == 2


@pytest.mark.parametrize("section", ["sequence", "fusion", "genome"])
def test_retired_yaml_and_dictionary_sections_fail(tmp_path, section):
    path = create_example(tmp_path / "input")
    with pytest.raises(PaceError, match="unknown keys"):
        load_config(path, overrides={section: {}})
    cfg = load_config(path)
    cfg[section] = {}
    with pytest.raises(PaceError, match="Unsupported configuration"):
        compute(cfg)


@pytest.mark.parametrize(
    "flag",
    [
        "--reference",
        "--vcf",
        "--sequence-model",
        "--fusion-model",
        "--predictions",
        "--callable",
        "--ploidy",
        "--unrecorded-site-policy",
    ],
)
def test_retired_file_flags_fail(tmp_path, flag):
    path = create_example(tmp_path / "input")
    with pytest.raises(SystemExit) as error:
        main(["run", "--config", str(path), flag, "retired", "--out", str(tmp_path / "out")])
    assert error.value.code == 2
    assert not (tmp_path / "out").exists()


def test_measured_output_round_trip_and_required_activity(tmp_path):
    path = create_example(tmp_path / "input")
    original = run(path, tmp_path / "first")
    cfg = load_config(path)
    cfg["inputs"]["observed_activity"] = None
    with pytest.raises(PaceError, match="every declared activity assay"):
        compute(cfg)
    for name in ("resolved_activity", "resolved_contacts", "sources", "evidence"):
        cfg["inputs"][name] = str(tmp_path / "first" / f"{name}.tsv")
    cfg["inputs"]["observed_contacts"] = None
    restored = run(cfg, tmp_path / "restored")
    assert [r["pace_score"] for r in restored["scores"]] == pytest.approx(
        [r["pace_score"] for r in original["scores"]]
    )
    assert "predictions" not in restored
    assert "genome_hashes" not in restored["manifest"]
    assert not (tmp_path / "first/predictions.tsv").exists()


def test_same_bin_import_uses_declared_resolution(tmp_path):
    cfg = load_config(create_example(tmp_path / "input"))
    tables = load_tables(cfg)
    exported = compute(cfg)
    tables["resolved_contacts"] = exported["resolved_contacts"]
    tables["observed_contacts"] = []
    tables["evidence"], tables["sources"] = exported["evidence"], exported["sources"]
    for row in tables["resolved_contacts"]:
        row.update(resolution=50000, bin_pair_id="0:0", resolved_value=3, observed_value=3)
    validate_tables(tables, cfg)
    actual = resolve_contacts(tables, cfg)
    assert all(math.isnan(r["resolved_value"]) for r in actual)
    assert all(r["reason"] == "near_diagonal_unresolved" for r in actual)


@pytest.mark.parametrize("kind", ["activity", "contacts"])
def test_imported_observations_require_the_correct_assay(tmp_path, kind):
    cfg = load_config(create_example(tmp_path / "input"))
    tables, exported = load_tables(cfg), compute(cfg)
    name = f"resolved_{kind}"
    tables[name] = copy.deepcopy(exported[name])
    tables["evidence"], tables["sources"] = exported["evidence"], exported["sources"]
    tables[name][0]["observation_sample_id"] = "S_Hi-C" if kind == "activity" else "S_ATAC"
    with pytest.raises(PaceError, match="assay"):
        validate_tables(tables, cfg)


@pytest.mark.parametrize("source", ["sequence_prediction", "fused", "contact_prior"])
def test_imported_activity_rejects_nonexperimental_sources(tmp_path, source):
    cfg = load_config(create_example(tmp_path / "input"))
    tables, exported = load_tables(cfg), compute(cfg)
    tables["resolved_activity"] = exported["resolved_activity"]
    tables["evidence"], tables["sources"] = exported["evidence"], exported["sources"]
    row = tables["resolved_activity"][0]
    row["evidence_type"] = source
    row["observation_sample_id"] = None
    for e in tables["evidence"]:
        if e["evidence_id"] == row["evidence_id"]:
            e["evidence_type"] = source
    with pytest.raises(PaceError):
        validate_tables(tables, cfg)


def test_auxiliary_measurement_definition_changes_contract(tmp_path):
    cfg = load_config(create_example(tmp_path / "input"))
    samples = read_table(cfg["inputs"]["samples"])
    samples.append({**samples[0], "sample_id": "AUX", "assay": "H3K4me1"})
    write_table(cfg["inputs"]["samples"], samples)
    observations = read_table(cfg["inputs"]["observed_activity"])
    observations.append({**observations[0], "sample_id": "AUX", "assay": "H3K4me1"})
    write_table(cfg["inputs"]["observed_activity"], observations)
    first = compute(cfg)
    observations[-1].update(unit="other_unit", normalization_id="other_protocol")
    write_table(cfg["inputs"]["observed_activity"], observations)
    second = compute(cfg)
    assert first["manifest"]["ml_feature_contract"] != second["manifest"]["ml_feature_contract"]
    assert [r["pace_score"] for r in first["scores"]] == [r["pace_score"] for r in second["scores"]]


def test_ml_report_evaluates_scope_and_final_calibration(tmp_path):
    path, rows, _ = learning_fixture(tmp_path)
    for row in rows:
        if row["split"] == "calibration":
            row["H3K4me1"] = 10 - row["H3K4me1"]
    write_table(tmp_path / "training.tsv", rows)
    report = train_classifier(path, tmp_path / "calibrated")
    assert report["test"]["auroc"] == 1
    assert report["test_probability"]["auroc"] == 0
    for row in rows:
        if row["split"] == "test":
            row["contact_sources"] = "contact_prior"
    write_table(tmp_path / "training.tsv", rows)
    report = train_classifier(path, tmp_path / "out_of_scope")
    assert report["test"]["coverage"] == 0
    assert report["test_status_counts"] == {"out_of_scope": 4}


def test_retired_modules_and_dependencies_are_absent():
    root = Path(__file__).resolve().parents[1]
    for relative in ("src/pace_livestock/io/variants.py", "src/pace_livestock/evidence/fusion.py"):
        assert not (root / relative).exists()
    assert not list((root / "src/pace_livestock/sequence").glob("*.py"))
    metadata = (root / "pyproject.toml").read_text()
    assert all(package not in metadata for package in ("torch", "safetensors", "pysam"))
    assert (
        json.loads(json.dumps(load_config(root / "examples/measured/config.yaml")))["regime"]
        == "measured"
    )
