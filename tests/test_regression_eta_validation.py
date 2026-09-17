"""Deployment gates for automatic allocation calibration on independent labels."""

import copy
import json

import pytest

from pace_livestock.config import load_config
from pace_livestock.demo import create_example
from pace_livestock.errors import PaceError
from pace_livestock.learning.allocation import (
    _independent_components,
    fit_from_labels,
    resolve_eta,
)
from pace_livestock.provenance import digest


@pytest.fixture
def cfg(tmp_path):
    return load_config(create_example(tmp_path / "inputs", "measured"))


def population(cfg, *, n=6, prefix="train", split="train", saturated=False, harmful=False):
    edges, labels = [], []
    for i in range(n):
        gene = f"{prefix}_G{i}"
        for positive in (True, False):
            element = f"{prefix}_E{i}_{int(positive)}"
            activity = (4 if saturated else 0.6) if positive else 1
            other_contact = 0.1 if positive else 9
            if harmful:
                activity = 1 if positive else 0.9
                other_contact = 9 if positive else 0.1
            edges.extend(
                [
                    dict(element_id=element, gene_id=gene, A_used=activity, Cbar=1),
                    dict(
                        element_id=element,
                        gene_id=f"{prefix}_other{i}",
                        A_used=activity,
                        Cbar=other_contact,
                    ),
                ]
            )
            labels.append(
                dict(
                    label_id=element,
                    element_id=element,
                    gene_id=gene,
                    species=cfg["context"]["species"],
                    assembly=cfg["context"]["assembly"],
                    context_id=cfg["context"]["context_id"],
                    target_level=cfg["target_level"],
                    perturbation_type="synthetic_inhibition",
                    effect_direction="down" if positive else "none",
                    label_status="enhancing_positive" if positive else "powered_negative",
                    split=split,
                    group_id=f"{prefix}_independent_locus{i}",
                    assay_id="toy_assay",
                    source_id="toy_source",
                )
            )
    return edges, labels


def test_single_group_does_not_deploy_nonzero_eta(cfg):
    edges, rows = population(cfg)
    for row in rows:
        row["group_id"] = "one_experiment"
    fitted = fit_from_labels(edges, rows, cfg)
    assert fitted["candidate_eta"] == 1
    assert fitted["eta"] == 0
    assert fitted["reason"] == "insufficient_independent_groups"
    assert fitted["deployment_validated"] is False


def test_no_ap_gain_keeps_zero_even_when_ranking_loss_improves(cfg):
    edges, rows = population(cfg, saturated=True)
    fitted = fit_from_labels(edges, rows, cfg)
    assert fitted["candidate_eta"] == 1
    assert fitted["eta"] == 0
    assert fitted["reason"] == "no_stable_independent_ap_gain"
    assert fitted["cross_validation"]["selected_shrinkage"] == 0


def test_independent_gain_selects_smallest_validated_nonzero_eta(cfg):
    edges, rows = population(cfg)
    fitted = fit_from_labels(edges, rows, cfg)
    assert fitted["eta"] == pytest.approx(0.25)
    assert fitted["deployment_validated"] is True
    cv = fitted["cross_validation"]
    assert cv["mean_ap_at_zero"] == 0.5
    assert cv["mean_ap_selected"] == 1
    assert cv["positive_group_fraction"] == 1
    for fold in cv["folds"]:
        assert not set(fold["training_group_ids"]) & set(fold["validation_group_ids"])
    assert fitted == fit_from_labels(list(reversed(edges)), list(reversed(rows)), cfg)


def test_test_labels_never_tune_eta(cfg):
    edges, rows = population(cfg)
    held_edges, held_labels = population(cfg, prefix="test", split="test", harmful=True)
    first = fit_from_labels(edges + held_edges, rows + held_labels, cfg)
    for row in held_labels:
        row["label_status"] = "powered_negative"
        row["effect_direction"] = "none"
    second = fit_from_labels(edges + held_edges, rows + held_labels, cfg)
    assert first == second
    assert first["excluded_counts"]["held_out_test"] == len(held_labels)


def test_calibration_confirms_but_does_not_fit_or_tune(cfg):
    edges, rows = population(cfg)
    held_edges, held_labels = population(cfg, prefix="cal", split="calibration", harmful=True)
    base = fit_from_labels(edges, rows, cfg)
    result = fit_from_labels(edges + held_edges, rows + held_labels, cfg)
    assert result["candidate_eta"] == base["candidate_eta"]
    assert result["cross_validation"] == base["cross_validation"]
    assert result["eta"] == 0
    assert result["reason"] == "calibration_gain_not_confirmed"
    assert result["estimation_label_ids"] == sorted(r["label_id"] for r in rows)
    assert result["confirmation_label_ids"] == sorted(r["label_id"] for r in held_labels)


def test_calibration_only_uses_grouped_cv(cfg):
    edges, rows = population(cfg, split="calibration")
    result = fit_from_labels(edges, rows, cfg)
    assert result["eta"] > 0
    assert result["estimation_split"] == "calibration_cross_validation"
    assert not result["confirmation_label_ids"]


def test_connected_identities_cannot_cross_folds():
    rows = [
        dict(label_id="a", gene_id="G1", element_id="E1", group_id="one"),
        dict(label_id="b", gene_id="G1", element_id="E2", group_id="two"),
        dict(label_id="c", gene_id="G2", element_id="E2", group_id="three"),
    ]
    assert len(_independent_components(rows)) == 1


def test_frozen_unvalidated_nonzero_eta_is_rejected(cfg, tmp_path):
    path = tmp_path / "calibrator.json"
    scope = {"context": cfg["context"]}
    artifact = dict(
        schema_version="pace-eta-2",
        kind="allocation_calibrator",
        scope=scope,
        scope_sha256=digest(scope),
        is_synthetic=True,
        status="fixed",
        eta=0.5,
    )
    path.write_text(json.dumps(artifact))
    cfg = copy.deepcopy(cfg)
    cfg["allocation"].update(eta="auto", calibrator_path=str(path))
    with pytest.raises(PaceError, match="validated grouped"):
        resolve_eta([], cfg, scope)


def test_declared_unpaired_confirmation_cannot_be_silently_dropped(cfg):
    edges, rows = population(cfg)
    held_edges, held_labels = population(cfg, prefix="cal", split="calibration")
    # A calibration set with positives only cannot validate a ranking change.
    held_labels = [row for row in held_labels if row["label_status"] == "enhancing_positive"]
    result = fit_from_labels(edges + held_edges, rows + held_labels, cfg)
    assert result["eta"] == 0
    assert result["reason"] == "insufficient_calibration_groups"


def test_unstable_group_gains_do_not_deploy(cfg):
    edges, rows = population(cfg, n=4)
    harmful_edges, harmful_rows = population(cfg, n=2, prefix="harmful", harmful=True)
    result = fit_from_labels(edges + harmful_edges, rows + harmful_rows, cfg)
    assert result["eta"] == 0
    assert result["reason"] == "no_stable_independent_ap_gain"


def test_legacy_calibrator_must_be_retrained(cfg, tmp_path):
    path = tmp_path / "legacy.json"
    path.write_text(
        json.dumps(
            dict(schema_version="pace-eta-1", kind="allocation_calibrator", eta=1, status="fitted")
        )
    )
    cfg["allocation"].update(eta="auto", calibrator_path=str(path))
    with pytest.raises(PaceError, match="schema/kind"):
        resolve_eta([], cfg, {})


def test_validated_frozen_nonzero_eta_can_be_reused(cfg, tmp_path):
    from pace_livestock.io.tables import write_table

    edges, rows = population(cfg)
    labels_path = tmp_path / "labels.tsv"
    write_table(labels_path, rows)
    cfg["allocation"].update(eta="auto", labels_path=str(labels_path))
    scope = {"context": cfg["context"]}
    fitted = resolve_eta(edges, cfg, scope)
    assert fitted["eta"] > 0
    path = tmp_path / "validated.json"
    path.write_text(json.dumps(fitted))
    cfg["allocation"].update(labels_path=None, calibrator_path=str(path))
    reloaded = resolve_eta(edges, cfg, scope)
    assert reloaded["eta"] == fitted["eta"]
    assert reloaded["reuse"] is True


def test_eta_scope_binds_contact_measurement_resolution(cfg):
    from pace_livestock.learning.allocation import calibration_scope
    from pace_livestock.schemas import load_tables

    tables = load_tables(cfg)
    first = calibration_scope(cfg, tables, {}, [])
    changed = copy.deepcopy(tables)
    for row in changed["observed_contacts"]:
        row["resolution"] *= 2
    second = calibration_scope(cfg, changed, {}, [])
    assert digest(first) != digest(second)
    assert (
        first["contact_measurement_contract"]["resolution"] * 2
        == second["contact_measurement_contract"]["resolution"]
    )


def test_eta_scope_is_portable_but_binds_reference_size_contents(cfg, tmp_path):
    from pace_livestock.learning.allocation import calibration_scope
    from pace_livestock.schemas import load_tables

    tables = load_tables(cfg)
    first_path = tmp_path / "first" / "chrom_sizes.tsv"
    moved_path = tmp_path / "relocated" / "chrom_sizes.tsv"
    for path in (first_path, moved_path):
        path.parent.mkdir()
        path.write_text("chrom\tlength\nchrToy\t1000000\n")
    original = copy.deepcopy(cfg)
    original["catalog"]["chrom_sizes_path"] = str(first_path)
    relocated = copy.deepcopy(original)
    relocated["catalog"]["chrom_sizes_path"] = str(moved_path)
    # Prior locations are likewise excluded; their manifest digests are bound
    # separately in asset_hashes whenever an actual prior is supplied.
    original["contact"]["prior_path"] = str(tmp_path / "first" / "prior.json")
    relocated["contact"]["prior_path"] = str(tmp_path / "relocated" / "prior.json")
    first = calibration_scope(original, tables, {}, [])
    moved = calibration_scope(relocated, tables, {}, [])
    assert first == moved
    assert str(first_path) not in json.dumps(first)
    moved_path.write_text("chrom\tlength\nchrToy\t1000001\n")
    changed = calibration_scope(relocated, tables, {}, [])
    assert digest(changed) != digest(first)
