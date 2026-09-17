"""Regression cases for contact contracts, boundary catalogs and log-space reloads."""

import json
import math
from pathlib import Path

import pytest
import yaml

from pace_livestock.config import load_config
from pace_livestock.demo import create_example
from pace_livestock.errors import PaceError
from pace_livestock.evaluation.compare import compare_runs, load_run
from pace_livestock.evidence.assets import load_asset
from pace_livestock.evidence.resolve import contact_measurement_contract, resolve_contacts
from pace_livestock.io.tables import read_table, write_table
from pace_livestock.operations import prepare_command
from pace_livestock.pipeline import compute, run
from pace_livestock.schemas import load_tables, universe_ids


def test_contact_replicates_require_one_resolution(tmp_path):
    path = create_example(tmp_path / "inputs", "measured")
    samples = read_table(path.parent / "samples.tsv")
    hic = next(row for row in samples if row["assay"] == "Hi-C")
    samples.append({**hic, "sample_id": "HIC_second", "technical_replicate": "2"})
    write_table(path.parent / "samples.tsv", samples)
    original = read_table(path.parent / "observed_contacts.tsv")
    second = [
        {**row, "sample_id": "HIC_second", "contact_value": 3 * float(row["contact_value"])}
        for row in original
    ]
    write_table(path.parent / "observed_contacts.tsv", original + second)
    compatible = compute(load_config(path))
    assert compatible["resolved_contacts"][0]["observed_value"] == pytest.approx(6)
    for row in second:
        row["resolution"] = 1000
    write_table(path.parent / "observed_contacts.tsv", original + second)
    with pytest.raises(PaceError, match="resolution"):
        compute(load_config(path))


@pytest.mark.parametrize("field", ["normalization_id", "balancing", "window_id"])
def test_contact_measurement_metadata_cannot_change_between_pairs(tmp_path, field):
    path = create_example(tmp_path / "inputs", "measured")
    rows = read_table(path.parent / "observed_contacts.tsv")
    rows[0][field] = "different_measurement"
    write_table(path.parent / "observed_contacts.tsv", rows)
    with pytest.raises(PaceError, match=field.replace("_id", "")):
        compute(load_config(path))


def test_prior_resolution_must_match_observed_contacts(tmp_path):
    path = create_example(tmp_path / "inputs", "measured")
    cfg = load_config(path)
    genome = load_config(create_example(tmp_path / "genome", "genome_only"))
    manifest = Path(genome["contact"]["prior_path"])
    if manifest.is_dir():
        manifest /= "manifest.json"
    data = json.loads(manifest.read_text())
    data["resolution"] = 10000
    manifest.write_text(json.dumps(data))
    cfg["contact"].update(
        prior_path=str(manifest),
        mode="shrinkage",
        reliability=0.5,
        reliability_source="independent_calibration",
    )
    with pytest.raises(PaceError, match="prior resolution"):
        compute(cfg)
    data.update(resolution=500, normalization_id="toy_mean")
    manifest.write_text(json.dumps(data))
    result = compute(cfg)
    assert all(math.isfinite(row["pace_score"]) for row in result["scores"])
    assert (
        result["manifest"]["comparison_contract"]["contact_measurement_contract"][
            "normalization_id"
        ]
        == "toy_mean"
    )


def test_source_normalization_is_checked_even_without_contact_row_metadata(tmp_path):
    path = create_example(tmp_path / "inputs", "measured")
    sources = read_table(path.parent / "sources.tsv")
    sources.append({**sources[0], "source_id": "different", "normalization_id": "counts_raw"})
    write_table(path.parent / "sources.tsv", sources)
    rows = read_table(path.parent / "observed_contacts.tsv")
    rows[0]["source_id"] = "different"
    write_table(path.parent / "observed_contacts.tsv", rows)
    with pytest.raises(PaceError, match="normalization"):
        compute(load_config(path))


def test_imported_contacts_need_a_measurement_contract(tmp_path):
    path = create_example(tmp_path / "inputs", "measured")
    cfg = load_config(path)
    tables = load_tables(cfg)
    tables["observed_contacts"] = []
    tables["resolved_contacts"] = [{"scale": cfg["contact"]["scale"]}]
    with pytest.raises(PaceError, match="declared resolution"):
        contact_measurement_contract(tables, cfg)
    tables["resolved_contacts"][0]["resolution"] = 500
    assert contact_measurement_contract(tables, cfg)["resolution"] == 500
    tables["resolved_contacts"].append({"scale": cfg["contact"]["scale"]})
    # Missing metadata cannot be borrowed from the first imported row: sorting
    # input must neither create nor remove knowledge about its measurement.
    for rows in (tables["resolved_contacts"], tables["resolved_contacts"][::-1]):
        tables["resolved_contacts"] = rows
        with pytest.raises(PaceError, match="declared resolution"):
            contact_measurement_contract(tables, cfg)


def test_different_resolutions_cannot_produce_a_biological_delta(tmp_path):
    path = create_example(tmp_path / "inputs", "measured")
    run(path, tmp_path / "run500")
    rows = read_table(path.parent / "observed_contacts.tsv")
    for row in rows:
        row["resolution"] = 1000
        if row["element_id"] == "E1":
            row["contact_value"] = 2 * float(row["contact_value"])
    write_table(path.parent / "observed_contacts.tsv", rows)
    run(path, tmp_path / "run1000")
    with pytest.raises(PaceError, match="contact_measurement_contract"):
        compare_runs(tmp_path / "run500", tmp_path / "run1000")
    # The evidence-policy override must not bypass a measurement mismatch.
    with pytest.raises(PaceError, match="contact_measurement_contract"):
        compare_runs(tmp_path / "run500", tmp_path / "run1000", allow_evidence_difference=True)


def test_underflow_and_true_zero_survive_result_reload(tmp_path):
    path = create_example(tmp_path / "inputs", "measured")
    activity = read_table(path.parent / "observed_activity.tsv")
    for row in activity:
        row["signal"] = 0 if row["element_id"] == "E1" else 1e-200
    write_table(path.parent / "observed_activity.tsv", activity)
    contacts = read_table(path.parent / "observed_contacts.tsv")
    for row in contacts:
        row["contact_value"] = 1e-200
    write_table(path.parent / "observed_contacts.tsv", contacts)
    result = run(path, tmp_path / "tiny")
    assert all(row["support"] == 0 for row in result["scores"])
    restored, _ = load_run(tmp_path / "tiny")
    for row in restored:
        assert (row["log_support"] == -math.inf) == (row["element_id"] == "E1")
        if row["element_id"] != "E1":
            assert math.isfinite(row["log_support"])
    compared, _, _ = compare_runs(tmp_path / "tiny", tmp_path / "tiny")
    assert all(row["reason"] == "complete" for row in compared)
    assert all(row["full_delta_pace"] == pytest.approx(0) for row in compared)


def test_evidence_policy_override_does_not_report_absolute_biological_deltas(tmp_path):
    path = create_example(tmp_path / "inputs", "measured")
    run(path, tmp_path / "original")
    cfg = yaml.safe_load(path.read_text())
    cfg["activity"]["minimum_callable_fraction"] = 0.5
    path.write_text(yaml.safe_dump(cfg))
    run(path, tmp_path / "stricter")
    compared, _, _ = compare_runs(
        tmp_path / "original", tmp_path / "stricter", allow_evidence_difference=True
    )
    assert all(row["conditional_delta_pace"] == pytest.approx(0) for row in compared)
    assert all(math.isnan(row["full_delta_pace"]) for row in compared)
    assert all(math.isnan(row["delta_support"]) for row in compared)


def test_imported_observation_cannot_bypass_near_diagonal_policy(tmp_path):
    path = create_example(tmp_path / "inputs", "measured")
    cfg = load_config(path)
    tables = load_tables(cfg)
    tables["resolved_contacts"] = resolve_contacts(tables, cfg)
    cfg["contact"]["near_diagonal_bp"] = 100000
    result = resolve_contacts(tables, cfg)
    assert all(math.isnan(row["resolved_value"]) for row in result)
    assert all(row["resolution_status"] == "unresolved" for row in result)


def test_imported_contacts_cannot_override_configured_evidence_policy(tmp_path):
    cfg = load_config(create_example(tmp_path / "inputs", "measured"))
    genome = load_config(create_example(tmp_path / "genome", "genome_only"))
    prior = load_asset(genome["contact"]["prior_path"], cfg, kind="contact_prior")
    prior["normalization_id"] = "toy_mean"
    tables = load_tables(cfg)
    imported = resolve_contacts(tables, cfg)
    for row in imported:
        row.update(
            evidence_type="contact_prior",
            prior_id=prior["model_id"],
            observation_sample_id=None,
            reliability=0,
            resolved_mode="prior_only",
        )
    tables["resolved_contacts"] = imported
    with pytest.raises(PaceError, match="fallback policy"):
        resolve_contacts(tables, cfg, prior_asset=prior)
    cfg["contact"]["allow_prior_fallback"] = True
    assert all(
        row["evidence_type"] == "contact_prior"
        for row in resolve_contacts(tables, cfg, prior_asset=prior)
    )
    tables["resolved_contacts"][0]["reliability"] = 1
    with pytest.raises(PaceError, match="reliability conflicts"):
        resolve_contacts(tables, cfg, prior_asset=prior)


@pytest.mark.parametrize("offset,tss0", [(0, 1750), (100, 50)])
def test_terminal_tss_catalog_loads_with_verified_reference_bounds(tmp_path, offset, tss0):
    source = tmp_path / "source"
    source.mkdir()
    (source / "regions.bed").write_text("chr1\t600\t1000\tE\n")
    (source / "genes.gtf").write_text(
        f'chr1\ttoy\ttranscript\t{tss0 + 1}\t1800\t.\t+\t.\tgene_id "G"; transcript_id "T";\n'
    )
    write_table(source / "sizes.tsv", [{"chrom": "chr1", "length": 1800}])
    prepare = source / "prepare.yaml"
    prepare.write_text(
        yaml.safe_dump(
            dict(
                kind="catalog",
                bed="regions.bed",
                gtf="genes.gtf",
                chrom_sizes="sizes.tsv",
                source_id="S",
                offset=offset,
            )
        )
    )
    prepare_command(prepare, tmp_path / "prepared")
    cfg = load_config(create_example(tmp_path / "template", "measured"))
    cfg["inputs"] = {key: None for key in cfg["inputs"]}
    for key in ("units", "promoters", "candidates"):
        cfg["inputs"][key] = str(tmp_path / "prepared" / f"{key}.tsv")
    cfg["catalog"].update(
        include_promoter_units=True,
        offset_bp=offset,
        chrom_sizes_path=str(tmp_path / "prepared" / "chrom_sizes.tsv"),
    )
    tables = load_tables(cfg)
    assert len(tables["promoters"]) == 1
    assert all(row["end"] - row["start"] == 500 for row in tables["units"])
    cfg["catalog"]["chrom_sizes_path"] = None
    with pytest.raises(PaceError, match="Promoter unit absent"):
        load_tables(cfg)


def test_interior_missing_promoter_is_not_excused_by_reference_sizes(tmp_path):
    path = create_example(tmp_path / "inputs", "measured")
    cfg = load_config(path)
    write_table(tmp_path / "sizes.tsv", [{"chrom": "chrToy", "length": 50000}])
    cfg["catalog"].update(include_promoter_units=True, chrom_sizes_path=str(tmp_path / "sizes.tsv"))
    with pytest.raises(PaceError, match="Promoter unit absent"):
        load_tables(cfg)


def test_catalog_identity_uses_reference_content_not_file_location(tmp_path):
    path = create_example(tmp_path / "inputs", "measured")
    cfg = load_config(path)
    ids = []
    for name in ("first.tsv", "relocated.tsv"):
        sizes = tmp_path / name
        write_table(sizes, [{"chrom": "chrToy", "length": 50000}])
        cfg["catalog"]["chrom_sizes_path"] = str(sizes)
        ids.append(universe_ids(load_tables(cfg), cfg))
    assert ids[0] == ids[1]
