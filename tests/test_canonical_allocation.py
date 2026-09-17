"""Independent fractional-score and convex calibration checks, including leakage."""

import copy
import math

import numpy as np
import pytest
import yaml

from pace_livestock.config import load_config
from pace_livestock.core import score
from pace_livestock.demo import create_example
from pace_livestock.errors import PaceError
from pace_livestock.evaluation.benchmark import benchmark_command
from pace_livestock.io.tables import write_table
from pace_livestock.learning.allocation import fit_eta, fit_from_labels
from pace_livestock.pipeline import compute, run
from pace_livestock.provenance import read_json


def edges():
    return [
        dict(element_id=f"E{i + 1}", gene_id=g, A_used=a, Cbar=c)
        for g, cs in [("G1", [3, 2, 2]), ("G2", [1, 2, 6])]
        for i, (a, c) in enumerate(zip([4, 2, 1], cs, strict=True))
    ]


def labels():
    return [
        dict(
            label_id=f"{g}-{e}",
            element_id=e,
            gene_id=g,
            species="synthetic",
            assembly="toy_assembly",
            context_id="toy_tissue",
            target_level="individual",
            perturbation_type="synthetic_inhibition",
            effect_direction="down" if positive else "none",
            label_status="enhancing_positive" if positive else "powered_negative",
            split="calibration",
            group_id=g,
            assay_id="toy_assay",
            source_id="toy_source",
        )
        for g, positive_e in [("G1", "E1"), ("G2", "E3")]
        for e, positive in [(positive_e, True), ("E2", False)]
    ]


@pytest.fixture
def cfg(tmp_path):
    cfg = load_config(create_example(tmp_path / "inputs", "measured"))
    cfg["allocation"].update(eta="auto", minimum_genes=2)
    return cfg


def test_fractional_formula():
    rows, _ = score(edges(), eta=0.5)
    # Hand calculation: B_G1=(3/4,1/2,1/4), B_G2=(1/4,1/2,3/4).
    # AC*sqrt(B) gives (6sqrt(3),2sqrt(2),1) and (2,2sqrt(2),3sqrt(3)).
    a = np.array([6 * math.sqrt(3), 2 * math.sqrt(2), 1])
    b = np.array([2, 2 * math.sqrt(2), 3 * math.sqrt(3)])
    np.testing.assert_allclose(
        [r["pace_score"] for r in rows], np.r_[a / a.sum(), b / b.sum()], rtol=1e-13
    )


@pytest.mark.parametrize("eta", [-0.1, 1.1, math.nan, math.inf, True, "auto"])
def test_invalid_eta(eta):
    with pytest.raises(PaceError, match="eta"):
        score(edges(), eta=eta)


def test_fractional_zero_missing_and_extremes():
    data = edges()
    data[0]["A_used"] = 0
    assert score(data, eta=0.5)[0][0]["pace_score"] == 0
    data[3]["Cbar"] = math.nan
    assert math.isnan(score(data, eta=0.5)[0][0]["pace_score"])
    assert score(data, eta=0)[0][0]["pace_score"] == 0
    for row in data:
        row["Cbar"] = 0
    assert all(math.isnan(r["pace_score"]) for r in score(data, eta=0.5)[0])
    for row in data:
        row["Cbar"], row["A_used"] = 1e300, 1e300
    # Six identical supports across two genes: each gene has three equal shares.
    assert all(r["pace_score"] == pytest.approx(1 / 3) for r in score(data, eta=0.37)[0])


def test_known_interior_optimum_and_permutation():
    groups = {"G1": ([[-0.5, 0]], [[0, -1]]), "G2": ([[0.5, -1]], [[0, 0]])}
    fitted = fit_eta(groups, minimum_genes=2)
    # L = [softplus(0.5-eta)+softplus(eta-0.5)]/2 has its unique minimum at 0.5.
    assert fitted["eta"] == pytest.approx(0.5, abs=1e-10)
    assert fitted["loss_at_eta"] == pytest.approx(math.log(2), abs=1e-12)
    assert fitted == fit_eta(dict(reversed(list(groups.items()))), minimum_genes=2)


def test_endpoints_and_unidentifiable():
    # A single positive contrast makes loss strictly decreasing; reversing it increases loss.
    assert fit_eta({"G": ([[0, 0]], [[0, -1]])}, minimum_genes=1)["eta"] == 1
    assert fit_eta({"G": ([[0, -1]], [[0, 0]])}, minimum_genes=1)["eta"] == 0
    assert fit_eta({"G": ([[2, -1]], [[0, -1]])}, minimum_genes=1)["status"] == "fallback"


def test_label_fit_and_test_value_invariance(cfg):
    rows = labels()
    held = {
        **rows[0],
        "label_id": "held",
        "element_id": "heldE",
        "gene_id": "heldG",
        "group_id": "held",
        "split": "test",
    }
    first = fit_from_labels(edges(), rows + [held], cfg)
    changed = {**held, "label_status": "powered_negative", "effect_direction": "none"}
    second = fit_from_labels(edges(), rows + [changed], cfg)
    # Both genes favor eta=1 in the training loss, but shared elements leave no
    # independent validation groups. Automatic deployment must retain eta=0.
    assert first["eta"] == second["eta"] == 0
    assert first["candidate_eta"] == second["candidate_eta"] == 1
    assert first["fit_label_ids"] == second["fit_label_ids"]
    assert first["excluded_counts"]["held_out_test"] == 1


@pytest.mark.parametrize("field", ["gene_id", "element_id", "group_id"])
def test_split_leakage(cfg, field):
    rows = labels()
    held = {
        **rows[0],
        "label_id": "held",
        "element_id": "heldE",
        "gene_id": "heldG",
        "group_id": "held",
        "split": "test",
    }
    held[field] = rows[0][field]
    with pytest.raises(PaceError, match="leak"):
        fit_from_labels(edges(), rows + [held], cfg)


def test_inapplicable_data_falls_back(cfg):
    data = labels()
    for row in data:
        row["context_id"] = "different_tissue"
    result = fit_from_labels(edges(), data, cfg)
    assert result["eta"] == 0
    assert result["excluded_counts"] == {"context_mismatch": 4}
    cfg["allocation"]["minimum_genes"] = 3
    assert fit_from_labels(edges(), labels(), cfg)["status"] == "fallback"


def test_zero_evidence_and_duplicate_labels(cfg):
    data = edges()
    data[0]["A_used"] = 0
    result = fit_from_labels(data, labels(), cfg)
    assert result["eta"] == 0
    assert result["excluded_counts"]["zero_or_unresolved_support"] == 1
    with pytest.raises(PaceError, match="Duplicate"):
        fit_from_labels(edges(), labels() + labels(), cfg)


def test_automatic_fallback_frozen_reuse_and_scope(cfg, tmp_path):
    no_labels = compute(cfg)
    fixed = copy.deepcopy(cfg)
    fixed["allocation"]["eta"] = 0
    assert no_labels["scores"] == compute(fixed)["scores"]
    assert no_labels["eta_calibration"]["reason"] == "no_functional_labels"
    table = tmp_path / "eta.tsv"
    write_table(table, labels())
    cfg["allocation"]["labels_path"] = str(table)
    result = run(cfg, tmp_path / "fit")
    artifact = tmp_path / "fit/eta_calibration.json"
    assert read_json(artifact)["eta"] == result["manifest"]["comparison_contract"]["eta"]
    cfg["allocation"].update(labels_path=None, calibrator_path=str(artifact))
    frozen = compute(cfg)
    assert frozen["scores"] == result["scores"]
    assert frozen["eta_calibration"]["reuse"] is True
    cfg["contact"]["near_diagonal_bp"] = 1
    with pytest.raises(PaceError, match="scope"):
        compute(cfg)


def test_manual_continuous_config_and_conflicts(cfg):
    overrides = {"context": cfg["context"], "allocation": {"eta": 0.35}}
    assert load_config(overrides=overrides)["allocation"]["eta"] == 0.35
    overrides["allocation"]["labels_path"] = "labels.tsv"
    with pytest.raises(PaceError, match="auto"):
        load_config(overrides=overrides)


def test_invalid_calibration_artifact(cfg, tmp_path):
    path = tmp_path / "not_an_artifact.json"
    path.write_text("[]")
    cfg["allocation"]["calibrator_path"] = str(path)
    with pytest.raises(PaceError, match="schema/kind"):
        compute(cfg)


def test_benchmark_excludes_fitting_labels_and_keeps_baselines(cfg, tmp_path):
    calibration_path = tmp_path / "calibration.tsv"
    write_table(calibration_path, labels())
    cfg["allocation"]["labels_path"] = str(calibration_path)
    run_path = tmp_path / "run.yaml"
    run_path.write_text(yaml.safe_dump(cfg))
    benchmark_labels = [{**r, "assayed_region_id": r["element_id"]} for r in labels()]
    table = tmp_path / "benchmark.tsv"
    write_table(table, benchmark_labels)
    membership = tmp_path / "membership.tsv"
    write_table(membership, [dict(region_id=e, element_id=e) for e in ("E1", "E2", "E3", "E4")])
    config_path = tmp_path / "benchmark.yaml"
    config_path.write_text(
        yaml.safe_dump(
            dict(run_config=str(run_path), labels=str(table), region_membership=str(membership))
        )
    )
    with pytest.raises(PaceError, match="overlap eta fitting"):
        benchmark_command(config_path, tmp_path / "leak")
    # Independent assayed edges outside the candidate catalog are retained as missed predictions.
    held = {
        **benchmark_labels[0],
        "label_id": "held",
        "gene_id": "G4",
        "element_id": "E4",
        "assayed_region_id": "E4",
        "group_id": "held",
    }
    write_table(table, [held])
    report = benchmark_command(config_path, tmp_path / "independent")
    assert {r["method"] for r in report} == {
        "PACE_calibrated",
        "PACE_eta0",
        "PACE_eta1",
        "ABC_style_single_TSS",
        "negative_distance",
    }
    assert all(r["n_scored"] == 0 for r in report)
