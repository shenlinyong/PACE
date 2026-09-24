"""Independent numerical and file-based regressions for sparse experimental data."""

import json
import math
from pathlib import Path

import numpy as np
import pytest

from pace_livestock.cli import main
from pace_livestock.config import load_config
from pace_livestock.core import score
from pace_livestock.demo import create_example
from pace_livestock.errors import PaceError
from pace_livestock.evidence.resolve import resolve_contacts
from pace_livestock.io.tables import read_table, write_table
from pace_livestock.pipeline import compute, run
from pace_livestock.provenance import output_directory, output_policy
from pace_livestock.reporting import regional_scores
from pace_livestock.schemas import load_tables


def edges():
    return [
        dict(element_id=e, gene_id="G", A_used=1, Cbar=c)
        for e, c in [("E1", 2), ("E2", 1), ("E3", math.nan)]
    ]


def test_missing_denominator_is_not_a_complete_score():
    rows, _ = score(edges())
    assert all(math.isnan(r["pace_score"]) for r in rows)
    assert rows[0]["pace_score_conditional"] == pytest.approx(2 / 3)
    assert (rows[0]["pace_score_lo"], rows[0]["pace_score_hi"]) == pytest.approx((0, 2 / 3))
    assert (rows[2]["pace_score_lo"], rows[2]["pace_score_hi"]) == (0, 1)
    # A hidden support of 97 would make the first share 0.02, not 2/3.
    full = edges()
    full[2]["Cbar"] = 97
    actual, _ = score(full)
    assert actual[0]["pace_score"] == pytest.approx(0.02)
    assert rows[0]["pace_score_lo"] <= 0.02 <= rows[0]["pace_score_hi"]
    compatibility, _ = score(edges(), partial_policy="conditional")
    assert compatibility[0]["pace_score"] == pytest.approx(2 / 3)
    assert compatibility[0]["score_scope"] == "conditional_only"


def test_supplied_bounds_cover_all_feasible_completions():
    missing = edges()
    missing[2].update(
        support_lower=2, support_upper=4, bound_source="independent upper/lower assessment"
    )
    rows, _ = score(missing)
    assert (rows[0]["pace_score_lo"], rows[0]["pace_score_hi"]) == pytest.approx((2 / 7, 2 / 5))
    assert (rows[2]["pace_score_lo"], rows[2]["pace_score_hi"]) == pytest.approx((2 / 5, 4 / 7))
    for completion in [2, 2.5, 3, 4]:
        truth = [2 / (3 + completion), 1 / (3 + completion), completion / (3 + completion)]
        for row, value in zip(rows, truth, strict=True):
            assert row["pace_score_lo"] - 1e-12 <= value <= row["pace_score_hi"] + 1e-12


def test_bounds_zero_total_single_unknown_and_overflow():
    rows, _ = score([dict(element_id="E", gene_id="G", A_used=1, Cbar=math.nan)])
    assert rows[0]["pace_score_lo"] == rows[0]["pace_score_hi"] == 1
    rows, _ = score([dict(element_id="E", gene_id="G", A_used=0, Cbar=2)])
    assert math.isnan(rows[0]["pace_score_lo"])
    rows, _ = score([dict(element_id=e, gene_id="G", A_used=1e300, Cbar=1e300) for e in ("a", "b")])
    for r in rows:
        assert r["pace_score_lo"] == pytest.approx(0.5)
        assert r["pace_score_hi"] == pytest.approx(0.5)


def prior_asset():
    return dict(
        model_id="toy_fitted_prior",
        kind="contact_prior",
        is_synthetic=True,
        species="synthetic",
        assembly="toy_assembly",
        context_id="toy_tissue",
        target_level="individual",
        training_sources=["synthetic"],
        calibration_sources=[],
        test_sources=[],
        validation={},
        a=1.0,
        gamma=1.0,
        d_ref=1000.0,
        d_min=1000.0,
        scale="toy_contact",
        resolution=500,
        normalization_id="toy_mean",
    )


def configured(tmp_path):
    path = create_example(tmp_path / "input", "measured")
    cfg = load_config(path)
    p = tmp_path / "prior.json"
    p.write_text(json.dumps(prior_asset()))
    cfg["contact"]["prior_path"] = str(p)
    return cfg


def test_pseudocount_rescues_sampling_zero_but_not_missing_activity(tmp_path):
    cfg = configured(tmp_path)
    t = load_tables(cfg)
    t["observed_contacts"][0]["contact_value"] = 0
    rows = resolve_contacts(t, cfg, prior_asset=prior_asset())
    r = next(r for r in rows if r["element_id"] == "E1" and r["promoter_id"] == "P1")
    assert r["observed_value"] == 0
    assert r["pseudocount_value"] == pytest.approx(0.2)
    assert r["resolved_value"] == pytest.approx(0.2)
    assert r["evidence_type"] == "regularized"
    cfg["contact"]["pseudocount"] = "none"
    assert (
        next(
            r
            for r in resolve_contacts(t, cfg, prior_asset=prior_asset())
            if r["element_id"] == "E1" and r["promoter_id"] == "P1"
        )["resolved_value"]
        == 0
    )


def test_prior_and_neighbor_diagonal_correction(tmp_path):
    cfg = configured(tmp_path)
    t = load_tables(cfg)
    t["promoters"][0]["tss0"] = 5249
    t["observed_contacts"][0].update(
        contact_value=1000.0, near_diagonal_value=4.0, near_diagonal_method="neighbor_max"
    )
    # A supplied compatible prior applies near the diagonal even if general fallback is false.
    rows = resolve_contacts(t, cfg, prior_asset=prior_asset())
    first = next(r for r in rows if r["element_id"] == "E1" and r["promoter_id"] == "P1")
    assert first["resolved_value"] == 1 and first["reason"] == "near_diagonal_prior"
    rows = resolve_contacts(t, cfg)
    first = next(r for r in rows if r["element_id"] == "E1" and r["promoter_id"] == "P1")
    assert first["resolved_value"] == 4 and first["reason"] == "near_diagonal_neighbor_max"
    cfg["contact"]["near_diagonal_policy"] = "unresolved"
    rows = resolve_contacts(t, cfg, prior_asset=prior_asset())
    assert math.isnan(
        next(r for r in rows if r["element_id"] == "E1" and r["promoter_id"] == "P1")[
            "resolved_value"
        ]
    )


def test_same_bin_tss_reuses_one_measurement_without_fake_replicates(tmp_path):
    cfg = load_config(create_example(tmp_path / "input", "measured"))
    t = load_tables(cfg)
    p = t["promoters"][0]
    p["pi"] = 0.25
    t["promoters"].append({**p, "promoter_id": "P1b", "tss0": p["tss0"] + 1, "pi": 0.75})
    rows = resolve_contacts(t, cfg)
    a = next(r for r in rows if r["element_id"] == "E1" and r["promoter_id"] == "P1")
    b = next(r for r in rows if r["element_id"] == "E1" and r["promoter_id"] == "P1b")
    assert b["resolved_value"] == a["resolved_value"] == 3
    assert b["shared_tss_bin"] and "contact:E1:P1:S_Hi-C" in b["parent_evidence_ids"]


def test_regularized_contacts_roundtrip_and_tampering(tmp_path):
    cfg = configured(tmp_path)
    original = run(cfg, tmp_path / "result")
    cfg["inputs"]["observed_contacts"] = None
    for name in ("resolved_contacts", "evidence", "sources"):
        cfg["inputs"][name] = str(tmp_path / "result" / f"{name}.tsv")
    restored = compute(cfg)
    assert [r["pace_score"] for r in restored["scores"]] == pytest.approx(
        [r["pace_score"] for r in original["scores"]]
    )
    rows = read_table(cfg["inputs"]["resolved_contacts"])
    rows[0]["pseudocount_value"] = 999
    write_table(cfg["inputs"]["resolved_contacts"], rows)
    with pytest.raises(PaceError, match="pseudocount"):
        compute(cfg)


def test_region_aggregation_counts_each_unit_once():
    full = edges()
    full[2]["Cbar"] = 1
    rows, _ = score(full)
    membership = [
        dict(region_id="peak", source_id="experiment", element_id=e) for e in ("E1", "E1", "E2")
    ]
    out = regional_scores(rows, membership)
    assert out[0]["n_units"] == 2 and out[0]["pace_score"] == pytest.approx(0.75)


def test_human_prior_is_explicit_research_only(tmp_path):
    overrides = dict(
        context=dict(species="Gallus gallus", assembly="GRCg7w", context_id="liver"),
        contact=dict(mode="prior_only", prior_preset="abc_human"),
    )
    cfg = load_config(overrides=overrides)
    assert cfg["contact"]["scale"] == "relative_distance_contact"
    from pace_livestock.evidence.assets import builtin_contact_prior

    assert builtin_contact_prior(cfg)["transfer_status"] == "unvalidated_for_target_context"
    for profile in ("validated", "demonstration"):
        with pytest.raises(PaceError, match="research"):
            load_config(overrides={**overrides, "execution_profile": profile})


def test_force_preserves_results_and_failed_replacements(tmp_path):
    out = tmp_path / "run"
    with output_directory(out) as d:
        (d / "data.txt").write_text("old")
    with output_policy(force=True):
        with pytest.raises(RuntimeError):
            with output_directory(out) as d:
                (d / "data.txt").write_text("bad")
                raise RuntimeError("computation failed")
        assert (out / "data.txt").read_text() == "old"
        with output_directory(out) as d:
            (d / "data.txt").write_text("new")
    assert (out / "data.txt").read_text() == "new"
    assert (next(tmp_path.glob("run.backup-*")) / "data.txt").read_text() == "old"
    unsafe = tmp_path / "inputs"
    unsafe.mkdir()
    with output_policy(force=True), pytest.raises(PaceError, match="identifiable"):
        with output_directory(unsafe):
            pass


def test_human_prior_run_exports_transfer_source_and_retains_measured_activity(tmp_path):
    fixture = create_example(tmp_path / "fixture", "measured")
    cfg = load_config(
        fixture,
        overrides={
            "execution_profile": "research",
            "contact": {"mode": "prior_only", "prior_preset": "abc_human"},
            "inputs": {"observed_contacts": None},
        },
    )
    # Synthetic input measurements exercise the research code path, not biological accuracy.
    result = compute(cfg)
    prior = result["manifest"]["asset_manifests"]["contact_prior"]
    assert prior["prior_source"] == "abc_human_default"
    assert prior["transfer_status"] == "unvalidated_for_target_context"
    assert prior["source_species"] == "Homo sapiens"
    assert all(r["evidence_type"] == "contact_prior" for r in result["resolved_contacts"])
    assert all(math.isfinite(r["pace_score"]) for r in result["scores"])
    observations = read_table(cfg["inputs"]["observed_activity"])
    observations[0].update(signal=None, measurement_status="unmeasured")
    write_table(cfg["inputs"]["observed_activity"], observations)
    assert all(math.isnan(r["pace_score"]) for r in compute(cfg)["scores"])


def test_preparation_commands_and_count_normalization(tmp_path):
    cfg = create_example(tmp_path / "input", "measured")
    assert (
        main(["prepare-pairs", "--catalog-dir", str(cfg.parent), "--out", str(tmp_path / "pairs")])
        == 0
    )
    assert len(read_table(tmp_path / "pairs/pairs.tsv")) == 6
    write_table(
        tmp_path / "counts.tsv", [dict(element_id="E1", sample_id="S", assay="ATAC", count=100)]
    )
    write_table(tmp_path / "libs.tsv", [dict(sample_id="S", library_size=1_000_000)])
    assert (
        main(
            [
                "normalize-activity",
                "--counts",
                str(tmp_path / "counts.tsv"),
                "--library-sizes",
                str(tmp_path / "libs.tsv"),
                "--units",
                str(cfg.parent / "units.tsv"),
                "--out",
                str(tmp_path / "normalized"),
            ]
        )
        == 0
    )
    row = read_table(tmp_path / "normalized/observed_activity.tsv")[0]
    assert float(row["signal"]) == pytest.approx(0.2)  # 100 CPM over 500 bp
    assert (
        main(
            [
                "merge-tables",
                "--table",
                "observed_activity",
                "--inputs",
                str(tmp_path / "normalized/observed_activity.tsv"),
                str(tmp_path / "normalized/observed_activity.tsv"),
                "--out",
                str(tmp_path / "bad"),
            ]
        )
        == 2
    )
    assert not (tmp_path / "bad").exists()


def test_init_creates_short_config_without_inventing_activity(tmp_path):
    assert (
        main(
            [
                "init",
                "--species",
                "Gallus gallus",
                "--assembly",
                "GRCg7w",
                "--tissue",
                "liver",
                "--panel",
                "H3K27ac",
                "--prior-preset",
                "abc_human",
                "--out",
                str(tmp_path / "project"),
            ]
        )
        == 0
    )
    import os
    import shlex

    from pace_livestock.cli import build_parser
    from pace_livestock.run_options import config_from_args

    script = (tmp_path / "project/run.sh").read_text()
    assert "--config" not in script and not (tmp_path / "project/config.yaml").exists()
    command = shlex.split(script.split("pace run", 1)[1].replace("\\\n", " "))
    cwd = os.getcwd()
    os.chdir(tmp_path / "project")
    try:
        cfg = config_from_args(build_parser().parse_args(["run", *command]))
    finally:
        os.chdir(cwd)
    assert cfg["contact"]["mode"] == "prior_only"
    assert read_table(cfg["inputs"]["observed_activity"]) == []
    assert cfg["activity"]["panel"] == ["H3K27ac"]


def test_promoter_weight_command(tmp_path):
    promoters = [
        dict(
            gene_id="G",
            promoter_id=p,
            chrom="chr1",
            tss0=i + 1,
            strand="+",
            pi=0.5,
            pi_source="equal",
        )
        for i, p in enumerate(["p1", "p2"])
    ]
    write_table(tmp_path / "promoters.tsv", promoters)
    write_table(
        tmp_path / "signals.tsv",
        [dict(promoter_id="p1", signal=2), dict(promoter_id="p2", signal=6)],
    )
    assert (
        main(
            [
                "prepare-promoter-weights",
                "--promoters",
                str(tmp_path / "promoters.tsv"),
                "--signals",
                str(tmp_path / "signals.tsv"),
                "--assay",
                "H3K4me3",
                "--normalization-id",
                "CPM",
                "--out",
                str(tmp_path / "weights"),
            ]
        )
        == 0
    )
    assert [float(r["pi"]) for r in read_table(tmp_path / "weights/promoters.tsv")] == [0.25, 0.75]


def tiny_cooler(tmp_path):
    cooler = pytest.importorskip("cooler")
    pd = pytest.importorskip("pandas")
    path = tmp_path / "tiny.cool"
    bins = pd.DataFrame(
        {"chrom": ["chr1"] * 8, "start": np.arange(8) * 1000, "end": np.arange(1, 9) * 1000}
    )
    pixels = pd.DataFrame(
        {"bin1_id": [0, 0, 0, 0], "bin2_id": [0, 1, 2, 3], "count": [1000.0, 14.0, 6.0, 10 / 3]}
    )
    cooler.create_cooler(str(path), bins, pixels, dtypes={"count": "float64"})
    return path


def test_sparse_cooler_fit_includes_all_zero_opportunities(tmp_path):
    path = tiny_cooler(tmp_path)
    from pace_livestock.io.prior_fit import cooler_distance_bins

    bins, _ = cooler_distance_bins(
        path, balanced=False, min_distance=1000, max_distance=4000, n_bins=3
    )
    assert [r["n_pairs"] for r in bins] == [7, 6, 5]
    assert [r["mean_contact"] for r in bins] == pytest.approx([2, 1, 2 / 3])
    assert (
        main(
            [
                "fit-prior",
                "--cooler",
                str(path),
                "--no-balanced",
                "--species",
                "synthetic",
                "--assembly",
                "toy",
                "--tissue",
                "toy",
                "--synthetic",
                "--min-distance",
                "1000",
                "--max-distance",
                "4000",
                "--reference-distance",
                "1000",
                "--minimum-distance",
                "1000",
                "--distance-bins",
                "3",
                "--out",
                str(tmp_path / "fit"),
            ]
        )
        == 0
    )
    manifest = json.loads((tmp_path / "fit/manifest.json").read_text())
    assert manifest["a"] == pytest.approx(2)
    assert manifest["gamma"] == pytest.approx(1)
    assert manifest["validation"] == {}


def test_cooler_diagonal_keeps_raw_value_and_reports_neighbor_max(tmp_path):
    from pace_livestock.io.cooler import query_contacts

    path = tiny_cooler(tmp_path)
    rows = query_contacts(
        path,
        [dict(element_id="E", promoter_id="P", chrom="chr1", anchor0=100, tss0=200)],
        resolution=1000,
        balanced=False,
        missing_pixels_are_zero=True,
    )
    assert rows[0]["contact_value"] == 1000
    assert rows[0]["near_diagonal_value"] == 14


def test_balanced_prior_fit_excludes_invalid_bin_opportunities(tmp_path):
    from pace_livestock.io.prior_fit import cooler_distance_bins

    path = tiny_cooler(tmp_path)
    h5py = pytest.importorskip("h5py")
    with h5py.File(path, "r+") as handle:
        handle["bins"].create_dataset("weight", data=[1, 2, np.nan, 1, 0, 1, 2, 1])
    bins, _ = cooler_distance_bins(
        path, balanced=True, min_distance=1000, max_distance=4000, n_bins=3
    )
    # Valid pairs at offsets 1, 2, 3 are (01,56,67), (13,35,57), (03,36).
    assert [r["n_pairs"] for r in bins] == [3, 3, 2]
    assert [r["mean_contact"] for r in bins] == pytest.approx([28 / 3, 0, 5 / 3])


def test_quoted_false_cannot_silently_turn_missing_pixels_into_zero():
    from pace_livestock.io.bigwig import quantify_bigwig
    from pace_livestock.io.cooler import query_contacts

    with pytest.raises(PaceError, match="boolean"):
        quantify_bigwig("unused", [], missing_is_measured_zero="false")
    for options in (
        dict(balanced=False, missing_pixels_are_zero="false"),
        dict(balanced="false", missing_pixels_are_zero=False),
    ):
        with pytest.raises(PaceError, match="booleans"):
            query_contacts("unused", [], resolution=1000, **options)


def test_console_entries_do_not_collide_on_case_insensitive_filesystems():
    import tomllib

    metadata = tomllib.loads((Path(__file__).parents[1] / "pyproject.toml").read_text())
    scripts = metadata["project"]["scripts"]
    assert scripts == {
        "pace": "pace_livestock.cli:main",
        "pace-livestock": "pace_livestock.cli:main",
    }
    assert len({name.casefold() for name in scripts}) == len(scripts)
