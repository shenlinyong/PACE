"""Counterexamples from the September 2026 source review."""

import copy
import math

import pytest

from pace_livestock.boundary_prior import BoundaryIndex
from pace_livestock.chunked import run_by_chromosome
from pace_livestock.config import load_config
from pace_livestock.demo import create_example
from pace_livestock.errors import PaceError
from pace_livestock.evaluation.metrics import binary_metrics
from pace_livestock.evidence.resolve import resolve_contacts
from pace_livestock.io.tables import read_table, write_table
from pace_livestock.pipeline import compute
from pace_livestock.schemas import load_tables


@pytest.mark.parametrize("labels", [[0.2, 1], [-0.2, 1], [0, 1.9]])
def test_binary_metrics_reject_fractional_labels(labels):
    with pytest.raises(PaceError, match="binary"):
        binary_metrics(labels, [0.1, 0.9])


def two_chromosomes(tmp_path):
    cfg = load_config(create_example(tmp_path / "inputs", "measured"))
    for name in ("units", "promoters", "candidates", "observed_activity", "observed_contacts"):
        path = cfg["inputs"][name]
        original = read_table(path)
        duplicate = copy.deepcopy(original)
        for row in duplicate:
            for key in ("element_id", "gene_id", "promoter_id", "bin_pair_id"):
                if row.get(key):
                    row[key] = "second_" + row[key]
            if "chrom" in row:
                row["chrom"] = "chrSecond"
        write_table(path, original + duplicate)
    return cfg


@pytest.mark.parametrize("problem", ["activity_scale", "contact_resolution"])
def test_chunks_reject_the_same_incompatible_measurements_as_whole_run(tmp_path, problem):
    cfg = two_chromosomes(tmp_path)
    name = "observed_activity" if problem == "activity_scale" else "observed_contacts"
    rows = read_table(cfg["inputs"][name])
    for row in rows:
        if row["element_id"].startswith("second_"):
            if problem == "activity_scale":
                row["unit"] = "different_signal_unit"
            else:
                row["resolution"] = 1000
    write_table(cfg["inputs"][name], rows)
    with pytest.raises(PaceError):
        compute(cfg)
    output = tmp_path / "chunked"
    with pytest.raises(PaceError):
        run_by_chromosome(cfg, output, max_pairs=1)
    assert not output.exists()


def test_fitted_bin_prior_uses_same_coordinates_in_all_contact_modes(tmp_path):
    cfg = load_config(create_example(tmp_path / "inputs", "measured"))
    tables = load_tables(cfg)
    unit, promoter = tables["units"][0], tables["promoters"][0]
    left = unit["anchor0"] // 500 * 500 + 250
    right = promoter["tss0"] // 500 * 500 + 250
    # Boundary lies between the genomic TSS and the center of its measured bin.
    boundary = BoundaryIndex(
        [dict(chrom=unit["chrom"], position0=promoter["tss0"] + 100, strength=1)]
    )
    prior = dict(
        a=1.0,
        gamma=1.0,
        beta=1.0,
        d_ref=500.0,
        d_min=500.0,
        model_id="fit",
        scale=cfg["contact"]["scale"],
        resolution=500,
        fitting_coordinate_policy="bin_centers",
        normalization_id="toy_mean",
    )
    from unittest.mock import patch

    expected = 500 / abs(left - right) * math.exp(-1)
    for mode, reliability in (("prior_only", None), ("shrinkage", 0.5), ("observed", None)):
        cfg["contact"].update(mode=mode, reliability=reliability, reliability_source="test")
        with patch("pace_livestock.evidence.resolve.load_boundaries", return_value=boundary):
            rows = resolve_contacts(tables, cfg, prior_asset=prior)
        found = next(
            r
            for r in rows
            if r["element_id"] == unit["element_id"] and r["promoter_id"] == promoter["promoter_id"]
        )
        assert found["prior_value"] == pytest.approx(expected)
        assert found["prior_coordinate_policy"] == "bin_centers"


def test_identities_do_not_depend_on_chunk_size(tmp_path):
    """Chunk size is an execution detail; chunked and single runs stay comparable."""
    import json

    from pace_livestock.evaluation.compare import compare_runs
    from pace_livestock.pipeline import run

    cfg = two_chromosomes(tmp_path)
    run(copy.deepcopy(cfg), tmp_path / "single")
    run_by_chromosome(copy.deepcopy(cfg), tmp_path / "small", max_pairs=1)
    run_by_chromosome(copy.deepcopy(cfg), tmp_path / "large", max_pairs=1000)
    ids = [
        json.loads((tmp_path / name / "run_manifest.json").read_text())["universe_ids"]
        for name in ("single", "small", "large")
    ]
    assert ids[0] == ids[1] == ids[2]
    for left, right in (("small", "large"), ("single", "small")):
        rows, _, _ = compare_runs(tmp_path / left, tmp_path / right)
        assert all(r["full_delta_pace"] == 0 for r in rows if r["reason"] == "complete")


def test_renamed_copy_of_a_pixel_is_counted_once():
    from pace_livestock.boundary_prior import unique_count_pairs

    row = dict(
        measurement_status="observed",
        resolution=500,
        anchor0=5249,
        tss0=10000,
        raw_count=3,
        count_to_contact=1.0,
        contact_value=3,
        sample_id="S",
        bin_pair_id="a",
        chrom="c1",
    )
    other = dict(row, anchor0=20249, bin_pair_id="b")
    renamed = dict(row, bin_pair_id="renamed")
    assert (
        len(unique_count_pairs([row, other, renamed], resolution=500, boundaries=BoundaryIndex([])))
        == 2
    )
    conflicting = dict(renamed, raw_count=4, contact_value=4)
    with pytest.raises(PaceError, match="conflicting"):
        unique_count_pairs([row, conflicting], resolution=500, boundaries=BoundaryIndex([]))


@pytest.mark.parametrize("module", ["benchmark", "training", "eta"])
def test_one_label_rule_everywhere(module):
    from pace_livestock.catalog import map_labels
    from pace_livestock.labels import classify_label
    from pace_livestock.learning.model import training_labels

    negative_down = dict(label_status="powered_negative", effect_direction="down")
    typo = dict(label_status="powered_negative", effect_direction="donw")
    assert classify_label(negative_down) == (None, "negative_requires_no_effect")
    if module == "benchmark":
        membership = [dict(region_id="R", element_id="E")]
        usable, rejected = map_labels([dict(negative_down, assayed_region_id="R")], membership)
        assert not usable and rejected[0]["reason"] == "negative_requires_no_effect"
        with pytest.raises(PaceError, match="effect_direction"):
            map_labels([dict(typo, assayed_region_id="R")], membership)
    elif module == "training":
        usable, rejected = training_labels([negative_down])
        assert not usable and rejected[0]["reason"] == "negative_requires_no_effect"
        with pytest.raises(PaceError, match="effect_direction"):
            training_labels([typo])
    else:
        with pytest.raises(PaceError, match="effect_direction"):
            classify_label(typo)


def test_anchor_prior_is_identical_in_prior_only_and_per_pair_modes(tmp_path):
    cfg = load_config(create_example(tmp_path / "inputs", "measured"))
    contacts = read_table(cfg["inputs"]["observed_contacts"])
    for row in contacts:
        row.update(raw_count=row["contact_value"], count_to_contact=1.0)
    write_table(cfg["inputs"]["observed_contacts"], contacts)
    tables = load_tables(cfg)
    prior = dict(
        a=1.0,
        gamma=1.0,
        beta=0.0,
        d_ref=500.0,
        d_min=500.0,
        kappa=2.0,
        model_id="anchors",
        scale=cfg["contact"]["scale"],
        resolution=500,
        normalization_id="toy_mean",
        manifest_sha256="x",
    )
    values = {}
    for mode, reliability in (("prior_only", None), ("shrinkage", "per_pair")):
        cfg["contact"].update(mode=mode, reliability=reliability, reliability_source="test")
        rows = resolve_contacts(tables, cfg, prior_asset=prior)
        values[mode] = {(r["element_id"], r["promoter_id"]): r["prior_value"] for r in rows}
        assert {r["prior_coordinate_policy"] for r in rows} == {"genomic_anchors"}
    assert values["prior_only"] == values["shrinkage"]
