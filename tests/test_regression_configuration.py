"""Public configuration and prepared-catalog integration regressions."""

import argparse
import copy

import pytest
import yaml

from pace_livestock.config import load_config
from pace_livestock.demo import create_example
from pace_livestock.errors import PaceError
from pace_livestock.io.tables import read_table, write_table
from pace_livestock.pipeline import compute
from pace_livestock.run_options import add_run_options, config_from_args


@pytest.mark.parametrize("offset", [500, 1000])
def test_grid_offset_is_checked_even_when_promoters_are_excluded(tmp_path, offset):
    path = create_example(tmp_path / "data", "measured")
    with pytest.raises(PaceError, match="offset_bp"):
        load_config(
            path, overrides={"catalog": {"offset_bp": offset, "include_promoter_units": False}}
        )


@pytest.mark.parametrize("value", [[], {}, 12, ""])
def test_contact_definition_rejects_nonidentifiers(tmp_path, value):
    path = create_example(tmp_path / "data", "measured")
    with pytest.raises(PaceError, match="normalization_id"):
        load_config(path, overrides={"contact": {"normalization_id": value}})


def test_catalog_directory_reads_portable_metadata_and_honors_explicit_override(tmp_path):
    root = tmp_path / "data"
    path = create_example(root, "measured")
    write_table(root / "chrom_sizes.tsv", [{"chrom": "chrToy", "length": 20000}])
    (root / "run_catalog_config.yaml").write_text(
        yaml.safe_dump(
            {
                "catalog": {
                    "chrom_sizes_path": "chrom_sizes.tsv",
                    "include_promoter_units": False,
                    "candidate_radius_bp": 20000,
                }
            }
        )
    )
    parser = argparse.ArgumentParser()
    add_run_options(parser)
    args = parser.parse_args(
        [
            "--config",
            str(path),
            "--catalog-dir",
            str(root),
            "--candidate-radius",
            "10000",
            "--out",
            str(tmp_path / "out"),
        ]
    )
    cfg = config_from_args(args)
    assert cfg["catalog"]["chrom_sizes_path"] == str(root / "chrom_sizes.tsv")
    assert cfg["catalog"]["candidate_radius_bp"] == 10000
    assert compute(cfg)["qc"]["n_scoreable"] == 6


def test_contact_cannot_impersonate_activity_sample(tmp_path):
    path = create_example(tmp_path / "data", "measured")
    cfg = load_config(path)
    rows = read_table(cfg["inputs"]["observed_contacts"])
    samples = read_table(cfg["inputs"]["samples"])
    rows_sample = next(row["sample_id"] for row in samples if row["assay"] == "ATAC")
    for row in rows:
        row["sample_id"] = rows_sample
    write_table(cfg["inputs"]["observed_contacts"], rows)
    with pytest.raises(PaceError, match="contact assay"):
        compute(cfg)


def test_model_feature_contract_changes_with_scoring_definition(tmp_path):
    path = create_example(tmp_path / "data", "measured")
    cfg = load_config(path)
    first = compute(cfg)["manifest"]["ml_feature_contract"]
    changed = copy.deepcopy(cfg)
    changed["allocation"]["eta"] = 1
    second = compute(changed)["manifest"]["ml_feature_contract"]
    assert first != second
    assert first["allocation_eta"] == 0 and second["allocation_eta"] == 1
