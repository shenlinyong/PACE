"""Exercise the public command with real files, direct options and transactional failures."""

import json
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

from pace_livestock.config import load_config
from pace_livestock.demo import create_example
from pace_livestock.io.tables import read_table
from pace_livestock.pipeline import compute


def invoke(*args, cwd=None):
    return subprocess.run(
        [sys.executable, "-m", "pace_livestock", *map(str, args)],
        text=True,
        capture_output=True,
        cwd=cwd,
    )


@pytest.mark.parametrize(
    "mode,regime", [("measured", "measured"), ("hybrid", "hybrid"), ("genome", "genome_only")]
)
def test_direct_modes_without_yaml(tmp_path, mode, regime):
    root = tmp_path / "files with spaces"
    config = create_example(root, regime)
    reference = compute(load_config(config))
    argv = [
        "--mode",
        mode,
        "--catalog-dir",
        root,
        "--species",
        "synthetic",
        "--assembly",
        "toy_assembly",
        "--tissue",
        "toy_tissue",
        "--profile",
        "demonstration",
        "--contact-scale",
        "toy_contact",
        "--no-include-promoters",
        "--out",
        tmp_path / "result",
    ]
    if regime != "measured":
        argv += [
            "--sequence-model",
            root / "models/sequence",
            "--reference",
            root / "genome.fa",
            "--vcf",
            root / "sample.vcf",
            "--callable",
            root / "callable.bed",
            "--ploidy",
            root / "ploidy.tsv",
            "--sample-id",
            "toy_animal",
            "--individual-id",
            "toy_animal",
        ]
        if regime == "hybrid":
            argv += ["--fusion-model", root / "models/fusion"]
        else:
            argv += ["--contact-prior", root / "models/contact"]
    # Remove the YAML: the command must use only the flags and supplied data files.
    config.unlink()
    result = invoke(*argv)
    assert result.returncode == 0, result.stderr
    info = json.loads(result.stdout)
    assert info["eta"] == 0 and info["eta_status"] == "fallback"
    rows = read_table(tmp_path / "result/scores.tsv.gz")
    # Direct flags and the independently validated YAML pipeline describe the same inputs.
    assert [float(r["pace_score"]) for r in rows] == pytest.approx(
        [r["pace_score"] for r in reference["scores"]]
    )
    manifest = json.loads((tmp_path / "result/run_manifest.json").read_text())
    assert manifest["regime"] == regime


def test_alias_override_paths_and_version(tmp_path):
    root = tmp_path / "input"
    config = create_example(root, "measured")
    new_activity = tmp_path / "override.tsv"
    new_activity.write_bytes((root / "observed_activity.tsv").read_bytes())
    result = invoke(
        "measured",
        "--config",
        config,
        "--activity",
        "override.tsv",
        "--eta",
        "0.35",
        "--out",
        "out",
        cwd=tmp_path,
    )
    assert result.returncode == 0, result.stderr
    resolved = yaml.safe_load((tmp_path / "out/resolved_config.yaml").read_text())
    assert resolved["inputs"]["observed_activity"] == str(new_activity)
    assert resolved["inputs"]["units"] == str(root / "units.tsv")
    assert json.loads(result.stdout)["eta"] == 0.35
    from pace_livestock import __version__

    assert invoke("--version").stdout.strip() == f"PACE {__version__}"


def test_cli_errors_and_output_protection(tmp_path):
    config = create_example(tmp_path / "input", "measured")
    for extra in [["--eta", "1.5"], ["--eta", "nan"], ["--species", "wrong"]]:
        result = invoke("--config", config, "--out", tmp_path / "failed", *extra)
        assert result.returncode == 2, result.stdout
        assert "PACE error" in result.stderr
        assert not (tmp_path / "failed").exists()
    assert invoke("run", "--config", config, "--out", tmp_path / "valid").returncode == 0
    protected = invoke("run", "--config", config, "--out", tmp_path / "valid")
    assert protected.returncode == 2 and "already exists" in protected.stderr
    assert invoke("fit-eta", "--config", config, "--out", tmp_path / "no_labels").returncode == 2


def test_installed_executable():
    # This test runs against a pip-installed package in CI, as well as the local editable install.
    executable = Path(sys.executable).parent / "PACE"
    assert executable.is_file(), "Install this revision with pip install -e . to create PACE"
    result = subprocess.run([str(executable), "--version"], capture_output=True, text=True)
    assert result.returncode == 0
    from pace_livestock import __version__

    assert result.stdout.strip() == f"PACE {__version__}"
