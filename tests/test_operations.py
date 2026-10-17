"""Actual command paths for calibration, preparation, benchmarks and measured evidence."""

import math

import pytest
import yaml
from fixture import create_example

from pace_livestock.cli import main
from pace_livestock.errors import PaceError
from pace_livestock.evaluation.benchmark import benchmark_command
from pace_livestock.io.tables import read_table, write_table
from pace_livestock.operations import fit_contact_command, prepare_command


def metadata():
    return dict(
        model_id="test",
        species="synthetic",
        assembly="toy",
        context_id="tissue",
        target_level="individual",
        is_synthetic=True,
    )


def test_fit_contact_includes_zeros(tmp_path):
    # Two pairs at each distance, one observed zero: means 4,2,1; a=4 at d_ref=1000, gamma=1.
    rows = [
        dict(
            bin_pair_id=f"{d}_{i}", distance_bp=d, contact_value=v, split="train", region_id="train"
        )
        for d, c in [(1000, 8), (2000, 4), (4000, 2)]
        for i, v in enumerate([c, 0])
    ]
    rows.append(
        dict(
            bin_pair_id="heldout",
            distance_bp=3000,
            contact_value=4 / 3,
            split="test",
            region_id="test",
        )
    )
    write_table(tmp_path / "contact.tsv", rows)
    cfg = dict(
        **metadata(),
        data="contact.tsv",
        scale="counts",
        resolution=500,
        bin_edges=[500, 1500, 3000, 5000],
        d_ref=1000,
        d_min=500,
    )
    p = tmp_path / "fit.yaml"
    p.write_text(yaml.safe_dump(cfg))
    asset = fit_contact_command(p, tmp_path / "model")
    assert asset["a"] == pytest.approx(4, abs=1e-10)
    assert asset["gamma"] == pytest.approx(1, abs=1e-10)
    residual = read_table(tmp_path / "model/held_out_residuals.tsv")[0]
    assert float(residual["residual"]) == pytest.approx(0, abs=1e-10)


def test_catalog_prepare_command(tmp_path):
    (tmp_path / "a.bed").write_text("chr1\t500\t1000\tE\n")
    (tmp_path / "a.gtf").write_text(
        'chr1\ttoy\ttranscript\t1501\t1900\t.\t+\t.\tgene_id "G"; transcript_id "T";\n'
    )
    write_table(tmp_path / "sizes.tsv", [dict(chrom="chr1", length=3000)])
    p = tmp_path / "prepare.yaml"
    p.write_text(
        yaml.safe_dump(
            dict(kind="catalog", bed="a.bed", gtf="a.gtf", chrom_sizes="sizes.tsv", source_id="S")
        )
    )
    result = prepare_command(p, tmp_path / "prepared")
    assert result == {"units": 2, "candidates": 2}


def test_benchmark_four_explicit_baselines(tmp_path):
    run_config = create_example(tmp_path / "inputs", "measured")
    write_table(
        tmp_path / "membership.tsv", [dict(region_id=e, element_id=e) for e in ("E1", "E2", "E3")]
    )
    labels = [
        dict(
            label_id=f"L{i}",
            assayed_region_id=e,
            gene_id=g,
            species="synthetic",
            assembly="toy_assembly",
            context_id="toy_tissue",
            effect_direction="down" if i % 2 else "none",
            label_status="enhancing_positive" if i % 2 else "powered_negative",
        )
        for i, (e, g) in enumerate([(e, g) for e in ("E1", "E2", "E3") for g in ("G1", "G2")])
    ]
    write_table(tmp_path / "labels.tsv", labels)
    p = tmp_path / "benchmark.yaml"
    p.write_text(
        yaml.safe_dump(
            dict(
                run_config=str(run_config), labels="labels.tsv", region_membership="membership.tsv"
            )
        )
    )
    report = benchmark_command(p, tmp_path / "benchmark")
    assert {r["method"] for r in report} == {
        "PACE_eta0",
        "PACE_eta1",
        "ABC_style_single_TSS",
        "negative_distance",
        "contact_power_2",
    }
    assert all(r["n_scored"] == 6 for r in report)
    assert all(
        math.isnan(r["precision"]) for r in report
    )  # No frozen deployment threshold supplied.


def test_cli_help_and_errors(tmp_path, capsys):
    cfg = create_example(tmp_path / "inputs", "measured")
    assert main(["validate", "--config", str(cfg)]) == 0
    assert main(["run", "--config", str(cfg), "--out", str(tmp_path / "run")]) == 0
    assert main(["run", "--config", str(cfg), "--out", str(tmp_path / "run")]) == 2
    assert "already exists" in capsys.readouterr().err


def benchmark_inputs(tmp_path, **label_fields):
    run_config = create_example(tmp_path / "inputs", "measured")
    write_table(
        tmp_path / "membership.tsv", [dict(region_id=e, element_id=e) for e in ("E1", "E2", "E3")]
    )
    labels = [
        dict(
            label_id=f"L{i}",
            assayed_region_id=e,
            gene_id=g,
            species="synthetic",
            assembly="toy_assembly",
            context_id="toy_tissue",
            effect_direction="down" if i % 2 else "none",
            label_status="enhancing_positive" if i % 2 else "powered_negative",
            **label_fields,
        )
        for i, (e, g) in enumerate([(e, g) for e in ("E1", "E2", "E3") for g in ("G1", "G2")])
    ]
    write_table(tmp_path / "labels.tsv", labels)
    p = tmp_path / "benchmark.yaml"
    p.write_text(
        yaml.safe_dump(
            dict(
                run_config=str(run_config), labels="labels.tsv", region_membership="membership.tsv"
            )
        )
    )
    return p, run_config


def test_benchmark_baselines_share_the_missing_value_policy(tmp_path):
    p, run_config = benchmark_inputs(tmp_path)
    cfg = yaml.safe_load(run_config.read_text())
    cfg.setdefault("scoring", {})["partial_policy"] = "conditional"
    run_config.write_text(yaml.safe_dump(cfg))
    contacts = run_config.parent / "observed_contacts.tsv"
    rows = read_table(contacts)
    write_table(contacts, [r for r in rows if r["bin_pair_id"] != "E1:P1"])
    report = benchmark_command(p, tmp_path / "benchmark")
    scored = {r["method"]: r["n_scored"] for r in report if r["scope"] == "all_tested_labels"}
    assert scored["PACE_eta0"] == scored["ABC_style_single_TSS"] == scored["contact_power_2"] == 5
    common = {r["method"]: r["n_tested"] for r in report if r["scope"] == "common_labels"}
    # eta = 1 also needs E1's contact with G1 to score E1-G2, so four labels are common.
    assert set(common.values()) == {4}


@pytest.mark.parametrize("field", ["species", "assembly"])
def test_benchmark_rejects_labels_of_another_genome(tmp_path, field):
    p, _ = benchmark_inputs(tmp_path)
    labels = read_table(tmp_path / "labels.tsv")
    write_table(tmp_path / "labels.tsv", [{**r, field: "other"} for r in labels])
    with pytest.raises(PaceError, match=field):
        benchmark_command(p, tmp_path / "benchmark")


def test_benchmark_rejects_membership_of_another_catalog(tmp_path):
    p, _ = benchmark_inputs(tmp_path)
    write_table(
        tmp_path / "membership.tsv", [dict(region_id=e, element_id=f"x{e}") for e in ("E1", "E2")]
    )
    with pytest.raises(PaceError, match="catalog"):
        benchmark_command(p, tmp_path / "benchmark")


def test_single_class_metrics_keep_the_defined_average_precision():
    from pace_livestock.evaluation.metrics import binary_metrics

    positives = binary_metrics([1, 1], [0.2, 0.5])
    assert positives["average_precision"] == 1 and math.isnan(positives["auroc"])
    negatives = binary_metrics([0, 0], [0.2, 0.5])
    assert math.isnan(negatives["average_precision"]) and math.isnan(negatives["auroc"])


@pytest.mark.parametrize("scale", [1e-200, 1.0, 1e200])
def test_contact_power_is_computed_in_log_space(scale):
    from pace_livestock.core.scoring import score

    edges = [
        dict(element_id=e, gene_id="G", A_used=a, Cbar=c * scale)
        for e, a, c in (("E1", 4.0, 3.0), ("E2", 2.0, 2.0), ("E3", 1.0, 1.0))
    ]
    rows, _ = score(edges, contact_power=2)
    support = [4 * 9, 2 * 4, 1]
    assert [r["pace_score"] for r in rows] == pytest.approx([s / sum(support) for s in support])


def test_benchmark_compares_external_methods_on_shared_labels(tmp_path):
    p, run_config = benchmark_inputs(tmp_path)
    external = [
        dict(element_id=e, gene_id=g, score=s)
        for (e, g), s in [(("E1", "G1"), 0.9), (("E2", "G1"), 0.1), (("E3", "G2"), 0.7)]
    ]
    write_table(tmp_path / "external.tsv", external)
    cfg = yaml.safe_load(p.read_text())
    cfg["external_methods"] = [
        dict(name="other", path="external.tsv", version="1", configuration="default")
    ]
    p.write_text(yaml.safe_dump(cfg))
    report = benchmark_command(p, tmp_path / "benchmark")
    shared = {r["method"]: r for r in report if r["scope"] == "common_labels_all_methods"}
    assert set(shared) >= {"other", "PACE_eta0", "ABC_style_single_TSS"}
    assert {r["n_tested"] for r in shared.values()} == {3}
    assert {r["n_scored"] for r in shared.values()} == {3}
