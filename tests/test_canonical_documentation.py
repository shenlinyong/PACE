"""Prevent recurrence of conflicting public formulas and executable entry points."""

import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("public_docs", ROOT / "scripts/check_public_docs.py")
CHECKER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CHECKER)


def test_public_documentation_contract():
    assert CHECKER.audit(ROOT) == []


@pytest.mark.parametrize("wrong", ["S/(sum(S)+U(G))", "The default eta=1.", "[gone](absent.md)"])
def test_audit_rejects_known_documentation_failures(tmp_path, wrong):
    page = tmp_path / "README.md"
    page.write_text(wrong)
    assert CHECKER.check_markdown(page, tmp_path)


def test_compatibility_launcher_uses_current_model(tmp_path):
    out = tmp_path / "run"
    command = [
        sys.executable,
        str(ROOT / "scripts/pace.py"),
        "measured",
        "--config",
        str(ROOT / "examples/measured/config.yaml"),
        "--eta",
        "0.5",
        "--out",
        str(out),
    ]
    result = subprocess.run(command, capture_output=True, text=True, cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    summary = json.loads(result.stdout)
    assert summary["eta"] == 0.5
    manifest = json.loads((out / "run_manifest.json").read_text())
    assert manifest["comparison_contract"]["eta"] == 0.5
    rejected = subprocess.run(
        [sys.executable, str(ROOT / "scripts/pace.py"), "--competition-power", "1"],
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    assert rejected.returncode != 0


def test_setup_uses_current_installer():
    result = subprocess.run(
        ["bash", str(ROOT / "setup.sh"), "--help"], capture_output=True, text=True
    )
    assert result.returncode == 0, result.stderr
    assert "--python" in result.stdout and "--extras" in result.stdout
