"""Independent objective checks, leakage guards, and trained-feature inference."""

import json
import math

import numpy as np
import pytest
import yaml

from pace_livestock.errors import PaceError
from pace_livestock.io.tables import write_table
from pace_livestock.learning.model import (
    fit_elastic_net,
    fit_preprocessor,
    predict_score_rows,
    train_classifier,
    training_labels,
    transform,
)


def test_intercept_is_not_penalized():
    model = fit_elastic_net(np.zeros((4, 1)), [0, 1, 1, 1], lambda1=10, lambda2=10)
    # Intercept-only Bernoulli maximum likelihood: log(3/1), even under heavy slope penalties.
    assert model["intercept"] == pytest.approx(math.log(3), abs=2e-6)
    assert model["coefficients"] == [0.0]
    assert model["converged"]


def test_training_only_imputation():
    x = np.array([[1.0, np.nan, 4.0], [3.0, np.nan, np.nan], [5.0, np.nan, 8.0]])
    p = fit_preprocessor(x)
    # Train medians 3 and 6; all-missing feature is removed, missing indicator retained for last.
    assert p["keep"] == [True, False, True]
    assert p["medians"] == [3.0, 6.0]
    before = json.dumps(p, sort_keys=True)
    transform(np.array([[1e100, 8.0, np.nan]]), p)
    assert json.dumps(p, sort_keys=True) == before


def test_label_exclusions():
    rows = [
        dict(label_status=s, effect_direction=d, mapping_count=m)
        for s, d, m in [
            ("enhancing_positive", "down", 1),
            ("enhancing_positive", "up", 1),
            ("low_power", "none", 1),
            ("not_tested", "none", 1),
            ("enhancing_positive", "down", 2),
        ]
    ]
    usable, rejected = training_labels(rows)
    assert len(usable) == 1 and usable[0]["label"] == 1
    assert len(rejected) == 4


def learning_fixture(tmp_path, *, calibrate=True):
    rows = []
    for split, groups in [("train", 6), ("calibration", 2), ("test", 2)]:
        for group in range(groups):
            for label in (0, 1):
                key = f"{split}_{group}_{label}"
                rows.append(
                    dict(
                        element_id=key,
                        gene_id="G",
                        assayed_region_id=key,
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
                        H3K4me1=10 * label,
                    )
                )
    write_table(tmp_path / "training.tsv", rows)
    context = dict(species="synthetic", assembly="toy", context_id="tissue")
    config = dict(
        data="training.tsv",
        extra_features=["H3K4me1"],
        penalties=[[0.01, 0.01], [0.02, 0.01]],
        folds=3,
        seed=17,
        model_id="test_classifier",
        calibrate=calibrate,
        is_synthetic=True,
        context=context,
        feature_contract={
            "target_level": "individual",
            "estimand": "bulk_proxy",
            "activity_panel": ["ATAC", "H3K27ac"],
            "activity_scales": [
                ["ATAC", "normalized", "toy", "w500"],
                ["H3K27ac", "normalized", "toy", "w500"],
            ],
            "contact_definition": {"scale": "normalized", "resolution_bp": 500},
            "candidate_construction": {"profile": "canonical_grid", "radius_bp": 5000000},
        },
    )
    path = tmp_path / "training.yaml"
    path.write_text(yaml.safe_dump(config))
    return path, rows, context


def test_train_calibrate_save_predict(tmp_path):
    path, rows, context = learning_fixture(tmp_path)
    report = train_classifier(path, tmp_path / "model")
    assert report["n_train"] == 12 and report["n_calibration"] == 4 and report["n_test"] == 4
    asset = json.loads((tmp_path / "model/model.json").read_text())
    assert asset["calibrator"] is not None
    assert any(abs(c) > 0 for c in asset["model"]["coefficients"])
    r = rows[-2:]
    features = [
        dict(
            entity_type="element",
            entity_id=x["element_id"],
            feature_name="H3K4me1",
            value=x["H3K4me1"],
        )
        for x in r
    ]
    predictions = predict_score_rows(
        r,
        features,
        tmp_path / "model",
        execution_profile="demonstration",
        context=context,
        feature_contract=yaml.safe_load(path.read_text())["feature_contract"],
    )
    # Identical core features; the actually fitted extra assay must affect inference.
    assert predictions[1]["pace_ml_score"] > predictions[0]["pace_ml_score"]
    assert all(math.isfinite(x["pace_ml_probability"]) for x in predictions)
    with pytest.raises(PaceError, match="Synthetic"):
        predict_score_rows(
            r, features, tmp_path / "model", execution_profile="research", context=context
        )


def test_test_extremes_cannot_change_fit(tmp_path):
    path, rows, _ = learning_fixture(tmp_path, calibrate=False)
    train_classifier(path, tmp_path / "a")
    for row in rows:
        if row["split"] == "test":
            row["H3K4me1"] = 1e100
    write_table(tmp_path / "training.tsv", rows)
    train_classifier(path, tmp_path / "b")
    a, b = [json.loads((tmp_path / d / "model.json").read_text()) for d in ("a", "b")]
    assert a["preprocessing"] == b["preprocessing"]
    assert a["model"] == b["model"]
    assert a["calibrator"] is None


def test_group_leakage_rejected(tmp_path):
    path, rows, _ = learning_fixture(tmp_path)
    rows[-1]["group_id"] = rows[0]["group_id"]
    write_table(tmp_path / "training.tsv", rows)
    with pytest.raises(PaceError, match="leaks"):
        train_classifier(path, tmp_path / "model")
