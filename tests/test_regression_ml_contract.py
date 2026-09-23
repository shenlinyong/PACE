"""Regressions for independent ML validation, scientific scope and feature quality."""

import copy
import json
import math
from pathlib import Path

import pytest
import yaml

from pace_livestock.config import load_config
from pace_livestock.demo import create_example
from pace_livestock.errors import PaceError
from pace_livestock.evaluation.benchmark import predict_ml_command
from pace_livestock.io.tables import read_table, write_table
from pace_livestock.learning.model import (
    canonical_feature_contract,
    independent_training_components,
    predict_score_rows,
    train_classifier,
)
from pace_livestock.pipeline import compute, run

CONTEXT = dict(species="synthetic", assembly="toy", context_id="liver")
CONTRACT = {
    "target_level": "individual",
    "estimand": "bulk_proxy",
    "activity_panel": ["ATAC", "H3K27ac"],
    "activity_scales": [
        ["ATAC", "normalized", "training_norm", "w500"],
        ["H3K27ac", "normalized", "training_norm", "w500"],
    ],
    "contact_definition": {
        "scale": "depth_normalized_contact",
        "resolution_bp": 500,
        "normalization_id": "depth_norm",
    },
    "candidate_construction": {
        "profile": "canonical_grid",
        "window_bp": 500,
        "radius_bp": 5000000,
        "include_promoter_units": True,
    },
}


def fixture_data(tmp_path):
    rows = []
    for split, n in (("train", 6), ("calibration", 2), ("test", 2)):
        for group in range(n):
            for label in (0, 1):
                element = f"{split}_{group}_{label}"
                rows.append(
                    dict(
                        element_id=element,
                        gene_id="G",
                        assayed_region_id=element,
                        mapping_count=1,
                        group_id=f"{split}_{group}",
                        split=split,
                        label_status="enhancing_positive" if label else "powered_negative",
                        effect_direction="down" if label else "none",
                        A_used=1,
                        Cbar=1,
                        distance_bp=100,
                        pace_score=0.5,
                        regime="measured",
                        activity_sources="observed",
                        contact_sources="observed",
                        **{"RNA:TPM": label * 10},
                    )
                )
    config = dict(
        data="labels.tsv",
        model_id="contract_test",
        is_synthetic=True,
        context=CONTEXT,
        feature_contract=CONTRACT,
        extra_features=["RNA:TPM"],
        penalties=[[0.01, 0.01], [0.02, 0.01]],
        folds=3,
        seed=17,
        calibrate=True,
    )
    write_table(tmp_path / "labels.tsv", rows)
    path = tmp_path / "train.yaml"
    path.write_text(yaml.safe_dump(config))
    return path, config, rows


@pytest.fixture(scope="module")
def trained(tmp_path_factory):
    folder = tmp_path_factory.mktemp("ml_contract")
    path, config, rows = fixture_data(folder)
    train_classifier(path, folder / "model")
    return folder / "model", config, rows


def predict(model, rows, contract=CONTRACT, *, status="resolved", value=10):
    rows = copy.deepcopy(rows)
    features = [
        dict(entity_type="gene", entity_id="G", feature_name="RNA:TPM", value=value, status=status)
    ]
    return predict_score_rows(
        rows,
        features,
        model,
        execution_profile="demonstration",
        context=CONTEXT,
        feature_contract=contract,
    )


def test_shared_regions_collapse_apparent_inner_cv_groups(tmp_path):
    path, _, rows = fixture_data(tmp_path)
    for row in rows:
        if row["split"] == "train":
            row["assayed_region_id"] = "same_physical_region"
    write_table(tmp_path / "labels.tsv", rows)
    with pytest.raises(PaceError, match="Insufficient independent groups"):
        train_classifier(path, tmp_path / "model")


def test_shared_edges_collapse_apparent_inner_cv_groups(tmp_path):
    path, _, rows = fixture_data(tmp_path)
    for row in rows:
        if row["split"] == "train":
            row["element_id"] = "same_enhancer"
    write_table(tmp_path / "labels.tsv", rows)
    with pytest.raises(PaceError, match="Insufficient independent groups"):
        train_classifier(path, tmp_path / "model")


def test_fold_manifest_keeps_connected_groups_together(tmp_path):
    path, _, rows = fixture_data(tmp_path)
    # Sharing one assayed region joins both declared groups without removing records.
    rows[0]["assayed_region_id"] = rows[2]["assayed_region_id"]
    write_table(tmp_path / "labels.tsv", rows)
    report = train_classifier(path, tmp_path / "model")
    assignment = {r["group_id"]: r for r in report["fold_assignments"]}
    assert report["n_independent_training_groups"] == 5
    assert assignment["train_0"]["fold"] == assignment["train_1"]["fold"]
    for fold in range(3):
        training = [
            r for r in rows if r["split"] == "train" and assignment[r["group_id"]]["fold"] != fold
        ]
        validation = [
            r for r in rows if r["split"] == "train" and assignment[r["group_id"]]["fold"] == fold
        ]
        assert not {r["assayed_region_id"] for r in training} & {
            r["assayed_region_id"] for r in validation
        }
        assert not {(r["element_id"], r["gene_id"]) for r in training} & {
            (r["element_id"], r["gene_id"]) for r in validation
        }


def test_transitive_group_connections():
    rows = [
        dict(group_id="a", element_id="e1", gene_id="g", assayed_region_id="r1"),
        dict(group_id="b", element_id="e1", gene_id="g", assayed_region_id="r2"),
        dict(group_id="c", element_id="e3", gene_id="g", assayed_region_id="r2"),
    ]
    components = independent_training_components(rows, range(3))
    assert len(set(components.values())) == 1


@pytest.mark.parametrize(
    "field",
    [
        "target_level",
        "estimand",
        "activity_panel",
        "activity_scales",
        "contact_definition",
        "candidate_construction",
    ],
)
def test_different_scientific_contract_suppresses_ml(trained, field):
    model, _, rows = trained
    changed = copy.deepcopy(CONTRACT)
    changes = {
        "target_level": "population_mean",
        "estimand": "another_estimand",
        "activity_panel": ["H3K27ac"],
        "activity_scales": [["H3K27ac", "normalized", "different_norm", "w500"]],
        "contact_definition": {**CONTRACT["contact_definition"], "resolution_bp": 1000},
        "candidate_construction": {**CONTRACT["candidate_construction"], "radius_bp": 1000000},
    }
    changed[field] = changes[field]
    prediction = predict(model, rows[-1:], changed)[0]
    assert prediction["ml_status"] == "out_of_scope"
    assert math.isnan(prediction["pace_ml_score"])
    assert math.isnan(prediction["pace_ml_probability"])


def test_new_loci_with_same_contract_are_allowed(trained):
    model, _, rows = trained
    row = {
        **rows[-1],
        "element_id": "unseen_enhancer",
        "gene_id": "unseen_gene",
        "candidate_universe_id": "new_cohort",
    }
    prediction = predict(model, [row])[0]
    assert prediction["ml_status"] == "resolved"
    assert math.isfinite(prediction["pace_ml_score"])
    assert math.isfinite(prediction["pace_ml_probability"])
    assert prediction["ml_probability_scope"] == "synthetic_demonstration_only"


def test_contract_panel_order_is_irrelevant():
    reverse = copy.deepcopy(CONTRACT)
    reverse["activity_panel"].reverse()
    reverse["activity_scales"].reverse()
    assert canonical_feature_contract(reverse) == canonical_feature_contract(CONTRACT)


def test_missing_prediction_contract_is_explicit_and_uncalibrated(trained):
    model, _, rows = trained
    prediction = predict(model, rows[-1:], None)[0]
    assert prediction["ml_status"] == "demonstration_unverified_contract"
    assert math.isfinite(prediction["pace_ml_score"])
    assert math.isnan(prediction["pace_ml_probability"])


def test_legacy_research_model_cannot_silently_claim_scope(trained, tmp_path):
    model, _, rows = trained
    asset = json.loads((model / "model.json").read_text())
    asset.pop("feature_contract")
    asset["is_synthetic"] = False
    legacy = tmp_path / "legacy.json"
    legacy.write_text(json.dumps(asset))
    predictions = predict_score_rows(
        copy.deepcopy(rows[-1:]),
        [],
        legacy,
        execution_profile="research",
        context=CONTEXT,
        feature_contract=CONTRACT,
    )
    assert predictions[0]["ml_status"] == "unverified_contract"
    assert math.isnan(predictions[0]["pace_ml_score"])
    assert math.isnan(predictions[0]["pace_ml_probability"])


def test_research_training_requires_declared_contract(tmp_path):
    path, config, _ = fixture_data(tmp_path)
    config.pop("feature_contract")
    config["is_synthetic"] = False
    path.write_text(yaml.safe_dump(config))
    with pytest.raises(PaceError, match="explicit feature_contract"):
        train_classifier(path, tmp_path / "model")


@pytest.mark.parametrize("bad_status", ["invalid", "unresolved", "missing", "low_quality"])
def test_invalid_feature_value_has_no_influence(trained, bad_status):
    model, _, rows = trained
    low = predict(model, rows[-1:], status=bad_status, value=0)[0]
    high = predict(model, rows[-1:], status=bad_status, value=1e6)[0]
    missing = predict(model, rows[-1:], status="resolved", value=None)[0]
    assert low["pace_ml_score"] == high["pace_ml_score"] == missing["pace_ml_score"]
    valid = predict(model, rows[-1:], status="resolved", value=10)[0]
    assert valid["pace_ml_score"] != missing["pace_ml_score"]


def test_resolution_status_overrides_residual_valid_status(trained):
    model, _, rows = trained
    features = [
        dict(
            entity_type="gene",
            entity_id="G",
            feature_name="RNA:TPM",
            value=1e6,
            status="observed",
            resolution_status="unresolved",
        )
    ]
    actual = predict_score_rows(
        copy.deepcopy(rows[-1:]),
        features,
        model,
        execution_profile="demonstration",
        context=CONTEXT,
        feature_contract=CONTRACT,
    )[0]
    expected = predict(model, rows[-1:], value=None)[0]
    assert actual["pace_ml_score"] == expected["pace_ml_score"]


@pytest.mark.parametrize(
    "field,value",
    [
        ("target_level", "population_mean"),
        ("normalization_id", "unrelated_normalization"),
        ("contact_scale", "raw_counts"),
        ("contact_resolution_bp", 1000),
    ],
)
def test_conflicting_row_metadata_cannot_bypass_contract(trained, field, value):
    model, _, rows = trained
    prediction = predict(model, [{**rows[-1], field: value}])[0]
    assert prediction["ml_status"] == "out_of_scope"
    assert math.isnan(prediction["pace_ml_probability"])


def test_mixed_training_measurement_contracts_are_rejected(tmp_path):
    path, _, rows = fixture_data(tmp_path)
    rows[-1]["normalization_id"] = "different_scale"
    write_table(tmp_path / "labels.tsv", rows)
    with pytest.raises(PaceError, match="scientific metadata"):
        train_classifier(path, tmp_path / "model")


def test_calibration_cannot_change_the_classifier_evidence_domain(tmp_path):
    path, _, rows = fixture_data(tmp_path)
    for row in rows:
        if row["split"] == "calibration":
            row["activity_sources"] = "aggregate"
    write_table(tmp_path / "labels.tsv", rows)
    with pytest.raises(PaceError, match="Calibration evidence sources/regime"):
        train_classifier(path, tmp_path / "model")


def test_published_training_example_and_both_inference_paths(tmp_path):
    """The shipped example is matched to the manifest consumed by both public paths."""
    repo = Path(__file__).resolve().parents[1]
    example = repo / "examples/training/learning.yaml"
    config_path = create_example(tmp_path / "inputs", "measured")
    cfg = load_config(config_path)
    samples = read_table(cfg["inputs"]["samples"])
    samples.append({**samples[0], "sample_id": "S_H3K4me1", "assay": "H3K4me1"})
    write_table(cfg["inputs"]["samples"], samples)
    observations = read_table(cfg["inputs"]["observed_activity"])
    for element, value in (("E1", 0), ("E2", 5), ("E3", 10)):
        observations.append(
            {
                **observations[0],
                "element_id": element,
                "sample_id": "S_H3K4me1",
                "assay": "H3K4me1",
                "signal": value,
            }
        )
    write_table(cfg["inputs"]["observed_activity"], observations)
    baseline = run(cfg, tmp_path / "baseline")
    training_config = yaml.safe_load(example.read_text())
    actual_contract = baseline["manifest"]["ml_feature_contract"]
    assert canonical_feature_contract(
        training_config["feature_contract"]
    ) == canonical_feature_contract(actual_contract)
    train_classifier(example, tmp_path / "model")
    cfg["multiomics"] = {"mode": "ml", "model_path": str(tmp_path / "model")}
    integrated = compute(cfg)
    assert all(row["ml_status"] == "resolved" for row in integrated["scores"])
    assert all(
        row["ml_probability_scope"] == "synthetic_demonstration_only"
        for row in integrated["scores"]
    )
    probabilities = {r["element_id"]: r["pace_ml_probability"] for r in integrated["scores"]}
    assert probabilities["E1"] < probabilities["E2"] < probabilities["E3"]
    assert [r["pace_score"] for r in integrated["scores"]] == [
        r["pace_score"] for r in baseline["scores"]
    ]
    assert integrated["qc"]["multiomics_roles"]["H3K4me1"] == "active_ml"
    prediction_config = tmp_path / "predict.yaml"
    prediction_config.write_text(
        yaml.safe_dump(
            {
                "run": "baseline",
                "model": "model",
                "features": "baseline/multiomics_features.tsv.gz",
            }
        )
    )
    standalone = predict_ml_command(prediction_config, tmp_path / "predictions")
    assert [r["pace_ml_probability"] for r in standalone] == [
        r["pace_ml_probability"] for r in integrated["scores"]
    ]
    # PACE itself is a learned feature: changing allocation changes its estimand.
    cfg["allocation"]["eta"] = 1.0
    mismatch = compute(cfg)
    assert all(row["ml_status"] == "out_of_scope" for row in mismatch["scores"])
    assert all(math.isnan(row["pace_ml_probability"]) for row in mismatch["scores"])
    cfg["allocation"]["eta"] = 0.0
    cfg["methylation"]["minimum_coverage"] += 1
    mismatch = compute(cfg)
    assert all(row["ml_status"] == "out_of_scope" for row in mismatch["scores"])
