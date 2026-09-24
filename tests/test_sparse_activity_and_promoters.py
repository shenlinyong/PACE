"""Sparse-assay policies checked against independently calculated scores."""

import copy
import json
import math

import pytest

from pace_livestock.cli import main
from pace_livestock.config import load_config
from pace_livestock.core import activity, score
from pace_livestock.demo import create_example
from pace_livestock.errors import PaceError
from pace_livestock.evaluation.compare import compare_runs
from pace_livestock.evidence.assets import load_asset
from pace_livestock.io.tables import read_table, write_table
from pace_livestock.pipeline import compute, run


def setup_run(tmp_path):
    return load_config(create_example(tmp_path / "input", "measured"))


def test_activity_offsets_preserve_na_and_require_scaled_offsets():
    assert activity([0, 9]) == 0
    assert activity([0, 9], [1, 0]) == pytest.approx(3)
    assert math.isnan(activity([math.nan, 9], [1, 0]))
    for offsets in ([math.nan, 0], [-1, 0], [1]):
        with pytest.raises(PaceError):
            activity([0, 9], offsets)
    original = [activity(x, [1, 0.5]) for x in ([0, 2], [4, 8])]
    scaled = [activity(x, [10, 0.5]) for x in ([0, 2], [40, 8])]
    assert original[0] / sum(original) == pytest.approx(scaled[0] / sum(scaled))


def test_activity_offsets_are_applied_once_after_aggregation(tmp_path):
    cfg = setup_run(tmp_path)
    data = read_table(cfg["inputs"]["observed_activity"])
    for row in data:
        if row["element_id"] == "E1" and row["assay"] == "ATAC":
            row["signal"] = 0
    write_table(cfg["inputs"]["observed_activity"], data)
    assert all(r["A_used"] == 0 for r in compute(cfg)["scores"] if r["element_id"] == "E1")
    cfg["activity"]["pseudocounts"] = {"ATAC": 1}
    result = run(cfg, tmp_path / "smoothed")
    raw = next(
        r for r in result["resolved_activity"] if r["element_id"] == "E1" and r["assay"] == "ATAC"
    )
    assert raw["resolved_value"] == 0 and raw["activity_pseudocount"] == 1
    other = next(
        r["resolved_value"]
        for r in result["resolved_activity"]
        if r["element_id"] == "E1" and r["assay"] == "H3K27ac"
    )
    assert next(r["A_used"] for r in result["scores"] if r["element_id"] == "E1") == pytest.approx(
        math.sqrt(other)
    )
    for name in ("resolved_activity", "evidence", "sources"):
        cfg["inputs"][name] = str(tmp_path / "smoothed" / f"{name}.tsv")
    cfg["inputs"]["observed_activity"] = None
    assert [r["pace_score"] for r in compute(cfg)["scores"]] == pytest.approx(
        [r["pace_score"] for r in result["scores"]]
    )


@pytest.mark.parametrize(
    "section,value",
    [
        ("activity", {"pseudocounts": {"RNA": 1}}),
        ("activity", {"pseudocounts": {"ATAC": -1}}),
        ("activity", {"pseudocounts": {"ATAC": True}}),
        ("promoters", {"minimum_retained_weight": 0}),
        ("promoters", {"minimum_weight": 1.1}),
        ("promoters", {"missing_policy": "per_pair"}),
        ("contact", {"allow_cross_context_prior": "yes"}),
    ],
)
def test_reject_invalid_sparse_settings(tmp_path, section, value):
    cfg_path = create_example(tmp_path / "input", "measured")
    with pytest.raises(PaceError):
        load_config(cfg_path, overrides={section: value})


def add_alternative_tss(cfg, weight=0.05):
    promoters = read_table(cfg["inputs"]["promoters"])
    original = next(p for p in promoters if p["gene_id"] == "G1")
    original["pi"] = 1 - weight
    promoters.append({**original, "promoter_id": "P1alt", "tss0": 11000, "pi": weight})
    write_table(cfg["inputs"]["promoters"], promoters)
    contacts = read_table(cfg["inputs"]["observed_contacts"])
    for e, value in [("E2", 50), ("E3", 5)]:
        base = next(r for r in contacts if r["element_id"] == e and r["promoter_id"] == "P1")
        contacts.append(
            {**base, "promoter_id": "P1alt", "contact_value": value, "bin_pair_id": e + ":P1alt"}
        )
    write_table(cfg["inputs"]["observed_contacts"], contacts)


def test_missing_tss_filter_is_gene_wide_and_records_weight_coverage(tmp_path):
    cfg = setup_run(tmp_path)
    expected = compute(cfg)
    add_alternative_tss(cfg)
    strict = compute(cfg)
    assert all(math.isnan(r["pace_score"]) for r in strict["scores"] if r["gene_id"] == "G1")
    cfg["promoters"]["missing_policy"] = "drop_missing"
    result = run(cfg, tmp_path / "filtered")
    assert [r["pace_score"] for r in result["scores"]] == pytest.approx(
        [r["pace_score"] for r in expected["scores"]]
    )
    for row in result["scores"]:
        if row["gene_id"] == "G1":
            assert row["tss_retained_weight"] == pytest.approx(0.95)
            assert row["tss_dropped_ids"] == "P1alt"
            assert row["tss_contact_scope"] == "selected_tss_set"
    table = read_table(tmp_path / "filtered/promoter_weights.tsv")
    alt = next(r for r in table if r["promoter_id"] == "P1alt")
    assert float(alt["pi_original"]) == 0.05 and float(alt["pi"]) == 0
    assert int(alt["n_missing_candidate_contacts"]) == 1
    assert alt["tss_selection_reason"] == "missing_contact"
    assert result["qc"]["promoter_processing"]["n_filtered_genes"] == 1
    # Filtering low original weights uses the same global set even with strict missingness.
    cfg["promoters"].update(missing_policy="strict", minimum_weight=0.1)
    low_weight = compute(cfg)
    assert [r["pace_score"] for r in low_weight["scores"]] == pytest.approx(
        [r["pace_score"] for r in expected["scores"]]
    )


def test_tss_filter_fails_when_too_much_original_weight_is_lost(tmp_path):
    cfg = setup_run(tmp_path)
    add_alternative_tss(cfg, weight=0.5)
    cfg["promoters"]["missing_policy"] = "drop_missing"
    result = compute(cfg)
    rows = [r for r in result["scores"] if r["gene_id"] == "G1"]
    assert all(math.isnan(r["Cbar"]) and math.isnan(r["pace_score"]) for r in rows)
    assert all(r["reason"] == "insufficient_retained_tss_weight" for r in rows)
    assert result["qc"]["promoter_processing"]["n_failed_genes"] == 1


def test_changed_promoter_selection_or_activity_offsets_cannot_be_compared_as_full_effect(tmp_path):
    cfg = setup_run(tmp_path)
    run(cfg, tmp_path / "baseline")
    shifted = copy.deepcopy(cfg)
    shifted["activity"]["pseudocounts"] = {"ATAC": 1}
    run(shifted, tmp_path / "offset")
    with pytest.raises(PaceError, match="activity_pseudocounts"):
        compare_runs(tmp_path / "baseline", tmp_path / "offset", allow_evidence_difference=True)
    add_alternative_tss(cfg)
    cfg["promoters"]["missing_policy"] = "drop_missing"
    run(cfg, tmp_path / "filtered")
    with pytest.raises(PaceError, match="promoter"):
        compare_runs(tmp_path / "baseline", tmp_path / "filtered", allow_evidence_difference=True)


def transferred_prior(cfg):
    return dict(
        model_id="other_tissue",
        kind="contact_prior",
        is_synthetic=True,
        **{**cfg["context"], "context_id": "other_tissue"},
        target_level=cfg["target_level"],
        training_sources=["toy"],
        calibration_sources=[],
        test_sources=[],
        validation={},
        a=1,
        gamma=1,
        d_ref=1000,
        d_min=500,
        scale="toy_contact",
        resolution=500,
        normalization_id="toy_mean",
    )


def test_other_tissue_prior_requires_opt_in_and_retains_source(tmp_path):
    cfg = setup_run(tmp_path)
    prior = tmp_path / "prior.json"
    prior.write_text(json.dumps(transferred_prior(cfg)))
    cfg["contact"]["prior_path"] = str(prior)
    with pytest.raises(PaceError, match="context_id mismatch"):
        compute(cfg)
    cfg["contact"]["allow_cross_context_prior"] = True
    result = compute(cfg)
    asset = result["manifest"]["asset_manifests"]["contact_prior"]
    assert asset["source_context_id"] == asset["context_id"] == "other_tissue"
    assert asset["target_context_id"] == cfg["context"]["context_id"]
    assert asset["transfer_status"] == "cross_context_unvalidated" and asset["validation"] == {}
    assert all(
        r["prior_transfer_status"] == "cross_context_unvalidated"
        for r in result["resolved_contacts"]
    )
    for field in ("species", "assembly", "target_level"):
        altered = transferred_prior(cfg)
        altered[field] = "different"
        prior.write_text(json.dumps(altered))
        with pytest.raises(PaceError, match=field + " mismatch"):
            load_asset(prior, cfg, kind="contact_prior")
    altered = transferred_prior(cfg)
    altered["is_synthetic"] = False
    prior.write_text(json.dumps(altered))
    cfg["execution_profile"] = "validated"
    with pytest.raises(PaceError, match="context_id mismatch"):
        load_asset(prior, cfg, kind="contact_prior")


def test_zero_observed_subset_is_still_partial_when_other_support_is_missing():
    rows, summary = score(
        [
            dict(element_id="a", gene_id="G", A_used=1, Cbar=0),
            dict(element_id="b", gene_id="G", A_used=1, Cbar=math.nan),
        ]
    )
    assert summary[0]["normalization_status"] == "partial"
    assert all(math.isnan(r["pace_score"]) for r in rows)


def test_coarse_cooler_prior_default_starts_at_its_resolution(tmp_path):
    cooler = pytest.importorskip("cooler")
    pd = pytest.importorskip("pandas")
    bins = pd.DataFrame(
        {
            "chrom": ["chr1"] * 12,
            "start": [25000 * i for i in range(12)],
            "end": [25000 * (i + 1) for i in range(12)],
        }
    )
    pixels = pd.DataFrame(
        [(i, j, 120 // (j - i)) for i in range(12) for j in range(i + 1, 12)],
        columns=["bin1_id", "bin2_id", "count"],
    )
    path = tmp_path / "coarse.cool"
    cooler.create_cooler(str(path), bins, pixels)
    out = tmp_path / "fit"
    assert (
        main(
            [
                "fit-prior",
                "--cooler",
                str(path),
                "--no-balanced",
                "--species",
                "toy",
                "--assembly",
                "toy",
                "--tissue",
                "liver",
                "--max-distance",
                "300000",
                "--out",
                str(out),
                "--synthetic",
            ]
        )
        == 0
    )
    manifest = json.loads((out / "manifest.json").read_text())
    assert manifest["fit_range"][0] == manifest["d_min"] == 25000
    assert manifest["gamma"] == pytest.approx(1, abs=0.12)
