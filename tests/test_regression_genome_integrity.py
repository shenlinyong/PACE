"""Individual imports and multiallelic/SV regressions use explicit toy genomes."""

from __future__ import annotations

import copy
import math
from pathlib import Path

import pytest

from pace_livestock.config import load_config
from pace_livestock.errors import PaceError
from pace_livestock.evidence.assets import load_asset
from pace_livestock.evidence.resolve import resolve_activity
from pace_livestock.io.tables import write_table
from pace_livestock.io.variants import Reference, VariantIndex, relationship_affected, selected_alts
from pace_livestock.operations import prepare_genome_command
from pace_livestock.pipeline import compute, run
from pace_livestock.schemas import load_tables
from pace_livestock.sequence.binding import bind_predictions, genome_binding_id
from pace_livestock.sequence.genome import build_window, prepare_windows
from pace_livestock.sequence.model import predict_windows


@pytest.fixture
def cfg():
    return load_config(Path(__file__).parents[1] / "examples/genome_only/config.yaml")


def variant(**kwargs):
    return (
        dict(
            chrom="chrToy",
            pos0=5100,
            end=5101,
            ref="A",
            alts=("C",),
            gt=(1, 1),
            phased=False,
            phase_set=None,
            filter="PASS",
        )
        | kwargs
    )


def write_vcf(path, *, pos=7001, alt="C,AC", gt="1/1", info="."):
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chrToy,length=20000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttoy_animal\n"
        f"chrToy\t{pos}\t.\tA\t{alt}\t60\tPASS\t{info}\tGT\t{gt}\n"
    )
    return str(path)


def predictions_for(cfg):
    tables = load_tables(cfg)
    asset = load_asset(cfg["sequence"]["model_path"], cfg, kind="sequence")
    windows, _ = prepare_windows(cfg, tables["units"], input_length=asset["input_length"])
    binding = genome_binding_id(cfg, tables["units"], asset)
    return bind_predictions(predict_windows(windows, asset), binding), tables, asset, windows


def test_unselected_structural_alt_does_not_disable_snv_window(cfg):
    tables = load_tables(cfg)
    v = variant(alts=("C", "<DEL>"), end=18000)
    with Reference(cfg["genome"]["reference_path"]) as ref:
        result = build_window(
            ref, tables["units"][0], [v], input_length=8192, ploidy=2, callable_fraction=1
        )
    assert result["status"] == "resolved"
    assert result["mapping_status"] == "identity_or_snv"
    assert not relationship_affected(v)
    # END belongs to the unused deletion, not to the carried SNP.
    assert not VariantIndex([v]).query("chrToy", 12000, 12500)


def test_unselected_short_indel_preserves_contacts(cfg, tmp_path):
    baseline = compute(cfg)
    cfg["genome"]["variant_path"] = write_vcf(tmp_path / "multi.vcf")
    actual = compute(cfg)
    assert actual["qc"]["n_scoreable"] == baseline["qc"]["n_scoreable"] == 6
    assert all(r["resolution_status"] == "resolved" for r in actual["resolved_contacts"])


@pytest.mark.parametrize(
    "alt", ["A]chrToy:10001]", "A[chrToy:10001[", "]chrToy:10001]A", "[chrToy:10001[A"]
)
def test_both_breakend_endpoints_indexed(alt):
    v = variant(pos0=2000, end=2001, alts=(alt,), gt=(0, 1), phased=True)
    index = VariantIndex([v])
    assert index.query("chrToy", 2000, 2001)
    remote = index.query("chrToy", 10000, 10001)
    assert len(remote) == 1 and remote[0]["remote_breakend"]
    assert relationship_affected(remote[0])
    assert not VariantIndex([{**v, "gt": (0, 0)}]).query("chrToy", 10000, 10001)


def test_remote_breakend_at_promoter_disables_contact(cfg, tmp_path):
    cfg["genome"]["variant_path"] = write_vcf(
        tmp_path / "bnd.vcf", pos=2001, alt="A]chrToy:10001]", gt="0|1", info="SVTYPE=BND"
    )
    actual = compute(cfg)
    contact = next(
        r
        for r in actual["resolved_contacts"]
        if (r["element_id"], r["promoter_id"]) == ("E3", "P1")
    )
    assert math.isnan(contact["resolved_value"])
    assert contact["resolution_status"] == "unresolved"
    assert contact["structural_status"] == "reported_variant_affects_relationship"


def test_missing_genotype_is_uncertainty_not_an_alt_choice():
    v = variant(alts=("C", "A]chrToy:10001]"), gt=(None, None))
    assert selected_alts(v) == ()
    assert relationship_affected(v)
    assert VariantIndex([v]).query("chrToy", 10000, 10001)
    with pytest.raises(PaceError, match="allele index"):
        selected_alts(variant(gt=(3, 3)))


def test_remote_breakend_outside_reference_rejected(cfg, tmp_path):
    cfg["genome"]["variant_path"] = write_vcf(
        tmp_path / "bnd.vcf", pos=2001, alt="A]missing:10001]", gt="0|1"
    )
    with pytest.raises(PaceError, match="remote breakend"):
        compute(cfg)


def test_individual_import_requires_reference_and_binding(cfg, tmp_path):
    predictions, _, _, _ = predictions_for(cfg)
    path = tmp_path / "predictions.tsv"
    write_table(path, predictions)
    cfg["inputs"]["predictions"] = str(path)
    missing_reference = copy.deepcopy(cfg)
    missing_reference["genome"]["reference_path"] = None
    with pytest.raises(PaceError, match="reference"):
        compute(missing_reference)
    unbound = [{k: v for k, v in r.items() if k != "genome_binding_id"} for r in predictions]
    write_table(path, unbound)
    with pytest.raises(PaceError, match="genome_binding_id"):
        compute(cfg)


def test_matching_import_reuses_values_but_changed_vcf_rejected(cfg, tmp_path):
    baseline = compute(cfg)
    predictions, _, _, _ = predictions_for(cfg)
    path = tmp_path / "predictions.tsv"
    write_table(path, predictions)
    cfg["inputs"]["predictions"] = str(path)
    imported = compute(cfg)
    assert [r["pace_score"] for r in imported["scores"]] == [
        r["pace_score"] for r in baseline["scores"]
    ]
    cfg["genome"]["variant_path"] = write_vcf(
        tmp_path / "changed.vcf", pos=5001, alt="<DEL>", gt="0|1", info="END=11000"
    )
    with pytest.raises(PaceError, match="genome_binding_id"):
        compute(cfg)


def test_import_cannot_override_failed_window(cfg):
    predictions, tables, asset, windows = predictions_for(cfg)
    baseline = compute(cfg)
    tables["resolved_activity"] = baseline["resolved_activity"]
    windows[0].update(
        status="unresolved",
        reason="reported_sv_unsupported",
        structural_status="reported_sv_unsupported",
    )
    rows = resolve_activity(
        tables, cfg, predictions=predictions, sequence_asset=asset, windows=windows
    )
    affected = [r for r in rows if r["element_id"] == windows[0]["element_id"]]
    assert affected and all(math.isnan(r["resolved_value"]) for r in affected)
    assert all(r["resolution_status"] == "unresolved" for r in affected)
    assert all(r["reason"] == "reported_sv_unsupported" for r in affected)


def test_prediction_cannot_override_failed_window(cfg):
    predictions, tables, asset, windows = predictions_for(cfg)
    windows[0].update(status="unresolved", reason="uncallable_input_window")
    rows = resolve_activity(
        tables, cfg, predictions=predictions, sequence_asset=asset, windows=windows
    )
    affected = [r for r in rows if r["element_id"] == windows[0]["element_id"]]
    assert affected and all(math.isnan(r["resolved_value"]) for r in affected)
    assert all(r["reason"] == "uncallable_input_window" for r in affected)


def test_run_exports_reusable_bound_predictions(cfg, tmp_path):
    original = run(cfg, tmp_path / "original")
    cfg["inputs"]["predictions"] = str(tmp_path / "original/predictions.tsv")
    imported = compute(cfg)
    assert [r["pace_score"] for r in imported["scores"]] == [
        r["pace_score"] for r in original["scores"]
    ]


def test_predict_sequence_export_reimports(cfg, tmp_path):
    config_path = Path(__file__).parents[1] / "examples/genome_only/config.yaml"
    prepare_genome_command(config_path, tmp_path / "prepared", predict=True)
    cfg["inputs"]["predictions"] = str(tmp_path / "prepared/predictions.tsv")
    imported = compute(cfg)
    assert imported["qc"]["n_scoreable"] == 6


def test_empty_vcf_without_header_is_not_reference_call(tmp_path):
    from pace_livestock.io.variants import read_variants

    vcf = tmp_path / "empty.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n")
    with pytest.raises(PaceError, match="#CHROM"):
        read_variants(vcf)
