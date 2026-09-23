"""Behavioral contracts spanning strict ingestion, evidence resolution and final outputs."""

import copy
import math

import numpy as np
import pytest
import yaml

from pace_livestock.config import load_config, load_yaml
from pace_livestock.demo import create_example, demo
from pace_livestock.errors import PaceError
from pace_livestock.evaluation.compare import compare_rows, compare_runs
from pace_livestock.evaluation.metrics import binary_metrics
from pace_livestock.io.tables import read_table, write_table
from pace_livestock.pipeline import compute, run
from pace_livestock.provenance import output_directory
from pace_livestock.schemas import load_tables


@pytest.fixture
def measured(tmp_path):
    return create_example(tmp_path / "inputs", "measured")


@pytest.mark.parametrize("regime", ["measured"])
def test_complete_offline_modes(tmp_path, regime):
    demo(regime, tmp_path / regime)
    root = tmp_path / regime
    rows = read_table(root / "results/scores.tsv.gz")
    # Each gene has exactly the three declared candidates, and support shares sum to one.
    assert len(rows) == 6
    for gene in ("G1", "G2"):
        assert sum(float(r["pace_score"]) for r in rows if r["gene_id"] == gene) == pytest.approx(
            1, abs=1e-10
        )
    # The persisted inputs remain rerunnable after the demo transaction moves into place.
    cfg = load_config(root / "results/resolved_config.yaml")
    assert len(compute(cfg)["scores"]) == 6


def test_strict_config(measured):
    cfg = load_yaml(measured)
    cfg["activity"]["panle"] = ["ATAC"]
    measured.write_text(yaml.safe_dump(cfg))
    with pytest.raises(PaceError, match="unknown keys"):
        load_config(measured)
    measured.write_text("regime: measured\nregime: hybrid\n")
    with pytest.raises(PaceError, match="Duplicate YAML"):
        load_config(measured)


def test_missing_panel_not_replaced(measured):
    path = measured.parent / "observed_activity.tsv"
    rows = read_table(path)
    rows[0]["signal"], rows[0]["measurement_status"] = None, "unmeasured"
    write_table(path, rows)
    result = compute(load_config(measured))
    bad = [r for r in result["scores"] if r["element_id"] == "E1"]
    assert all(math.isnan(r["pace_score"]) for r in bad)
    assert all(r["normalization_status"] == "partial" for r in result["scores"])


def test_true_zero_retained(measured):
    path = measured.parent / "observed_activity.tsv"
    rows = read_table(path)
    rows[0]["signal"] = 0
    write_table(path, rows)
    result = compute(load_config(measured))
    assert all(
        r["support"] == 0 and r["pace_score"] == 0
        for r in result["scores"]
        if r["element_id"] == "E1"
    )


def test_order_invariance(measured):
    expected = compute(load_config(measured))
    for path in measured.parent.glob("*.tsv"):
        rows = read_table(path)
        if rows:
            write_table(path, rows[::-1])
    result = compute(load_config(measured))
    # Permutation of tables must not change numerical results or universe identities.
    np.testing.assert_allclose(
        [r["pace_score"] for r in result["scores"]],
        [r["pace_score"] for r in expected["scores"]],
        atol=1e-12,
        rtol=0,
    )
    assert result["manifest"]["universe_ids"] == expected["manifest"]["universe_ids"]


@pytest.mark.parametrize(
    "table,column,value",
    [
        ("samples", "assembly", "wrong"),
        ("observed_activity", "sample_id", "fake_sample"),
        ("observed_activity", "signal", "-1"),
        ("observed_activity", "signal", "nan"),
        ("observed_activity", "window_id", "grid:1000:mean"),
        ("units", "end", "5501"),
        ("promoters", "pi", ".5"),
    ],
)
def test_invalid_scientific_input(measured, table, column, value):
    path = measured.parent / f"{table}.tsv"
    rows = read_table(path)
    rows[0][column] = value
    write_table(path, rows)
    with pytest.raises(PaceError):
        compute(load_config(measured))


def test_shared_bin_integrity(measured):
    path = measured.parent / "observed_contacts.tsv"
    rows = read_table(path)
    rows[1]["bin_pair_id"] = rows[0]["bin_pair_id"]
    write_table(path, rows)
    with pytest.raises(PaceError, match="bin pair"):
        compute(load_config(measured))


def test_same_bin_contact_needs_explicit_correction(measured):
    path = measured.parent / "observed_contacts.tsv"
    rows = read_table(path)
    # All coordinates now fall in one huge observed contact bin; no prior is supplied.
    for row in rows:
        row["resolution"] = 100000
    write_table(path, rows)
    result = compute(load_config(measured))
    assert all(math.isnan(r["pace_score"]) for r in result["scores"])
    assert all(r["reason"] == "near_diagonal_unresolved" for r in result["resolved_contacts"])


def test_published_evidence_parents_are_resolvable(measured):
    result = compute(load_config(measured))
    ids = {r["evidence_id"] for r in result["evidence"]}
    sources = {r["source_id"] for r in result["sources"]}
    for row in result["evidence"]:
        assert row["source_id"] in sources
        assert all(p in ids for p in (row["parent_evidence_ids"] or "").split(";") if p)


def test_annotation_does_not_change_score(measured):
    cfg = load_config(measured)
    original = compute(cfg)
    samples = read_table(measured.parent / "samples.tsv")
    samples.append({**samples[0], "sample_id": "RNA", "assay": "RNA"})
    write_table(measured.parent / "samples.tsv", samples)
    expr = measured.parent / "expression.tsv"
    write_table(expr, [{"gene_id": "G1", "sample_id": "RNA", "tpm": 1e9, "status": "observed"}])
    cfg["inputs"]["expression"] = str(expr)
    result = compute(cfg)
    # The model contract leaves RNA outside the main multiplicative score.
    np.testing.assert_allclose(
        [r["pace_score"] for r in original["scores"]],
        [r["pace_score"] for r in result["scores"]],
        atol=1e-12,
        rtol=0,
    )
    assert result["qc"]["multiomics_roles"]["RNA"] == "annotation_only"


def test_transaction_and_compare(measured, tmp_path):
    original = run(measured, tmp_path / "a")
    with pytest.raises(PaceError, match="already exists"):
        run(measured, tmp_path / "a")
    run(measured, tmp_path / "b")
    rows, _, _ = compare_runs(tmp_path / "a", tmp_path / "b")
    assert all(r["full_delta_pace"] == pytest.approx(0, abs=1e-12) for r in rows)
    assert len(original["scores"]) == len(rows)
    with pytest.raises(RuntimeError):
        with output_directory(tmp_path / "failed") as directory:
            (directory / "partial").write_text("bad")
            raise RuntimeError("interrupt")
    assert not (tmp_path / "failed").exists()


def test_common_denominator_counterexample():
    from pace_livestock.core import score

    a, _ = score(
        [dict(element_id=f"E{i}", gene_id="G", A_used=s, Cbar=1) for i, s in enumerate([1, 1, 2])]
    )
    b, _ = score(
        [
            dict(element_id=f"E{i}", gene_id="G", A_used=s, Cbar=1)
            for i, s in enumerate([1, 1, math.nan])
        ]
    )
    rows = compare_rows(a, b)
    # Review §3.2: .25 vs .5 are conditional on different denominators; common E0/E1 both .5.
    assert rows[0]["original_score_a"] == 0.25
    assert math.isnan(rows[0]["original_score_b"])
    assert rows[0]["common_score_a"] == rows[0]["common_score_b"] == 0.5
    assert rows[0]["conditional_delta_pace"] == 0
    assert math.isnan(rows[0]["full_delta_pace"])


def test_composition_change_is_not_activity_change():
    from pace_livestock.core import score

    a, _ = score(
        [dict(element_id=e, gene_id="G", A_used=s, Cbar=1) for e, s in [("E1", 1), ("E2", 1)]]
    )
    b, _ = score(
        [dict(element_id=e, gene_id="G", A_used=s, Cbar=1) for e, s in [("E1", 1), ("E2", 2)]]
    )
    row = compare_rows(a, b)[0]
    # Review §3.3: unchanged support in E1, but share changes from 1/2 to 1/3.
    assert row["delta_support"] == 0
    assert row["delta_A"] == 0
    assert row["full_delta_pace"] == pytest.approx(-1 / 6, abs=1e-12)


def test_ap_ties_and_missing_predictions():
    # Tie-aware AP: first tied block contains 1 TP / 2 calls, final recall adds .5 at precision 2/3.
    result = binary_metrics([1, 0, 1, 1], [0.9, 0.9, 0.1, math.nan], threshold=0.5)
    assert result["average_precision"] == pytest.approx(7 / 12)
    assert result["auroc"] == 0.25
    assert result["missed_positive"] == 1
    assert result["recall_end_to_end"] == pytest.approx(1 / 3)


def test_schema_duplicate_key(measured):
    path = measured.parent / "units.tsv"
    rows = read_table(path)
    write_table(path, rows + [copy.deepcopy(rows[0])])
    with pytest.raises(PaceError, match="duplicate"):
        load_tables(load_config(measured))
