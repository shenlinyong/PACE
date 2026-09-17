"""Actual command paths for calibration, preparation, benchmarks and variant scenarios."""

import math

import pytest
import yaml

from pace_livestock.cli import main
from pace_livestock.demo import create_example
from pace_livestock.errors import PaceError
from pace_livestock.evaluation.benchmark import benchmark_command
from pace_livestock.io.tables import read_table, write_table
from pace_livestock.operations import fit_contact_command, fit_fusion_command, prepare_command
from pace_livestock.sequence.effects import variant_effects_command


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


def test_fit_fusion_independent_target(tmp_path):
    rows = [
        dict(
            assay="ATAC",
            quality_stratum="default",
            observed=9,
            predicted=3,
            target=math.sqrt(40) - 1,
            split="calibration" if i < 4 else "test",
            group_id=str(i),
            input_donor=f"D{i}",
            target_donor=f"D{i}",
            input_measurement_id=f"low{i}",
            target_measurement_id=f"independent{i}",
        )
        for i in range(6)
    ]
    write_table(tmp_path / "fusion.tsv", rows)
    cfg = dict(
        **metadata(),
        data="fusion.tsv",
        calibration_target="individual_state",
        signal_unit="signal",
        normalization_id="norm",
        output_window=500,
        scales={"ATAC": 1},
        minimum_samples=4,
        measurement_design="independent_technical_replicates_same_donor",
    )
    p = tmp_path / "fit.yaml"
    p.write_text(yaml.safe_dump(cfg))
    asset = fit_fusion_command(p, tmp_path / "model")
    assert asset["strata"]["default"]["ATAC"]["weight"] == pytest.approx(0.5, abs=1e-10)
    rows[0]["target_donor"] = "other_animal"
    write_table(tmp_path / "fusion.tsv", rows)
    with pytest.raises(PaceError, match="matching donors"):
        fit_fusion_command(p, tmp_path / "bad")


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
    }
    assert all(r["n_scored"] == 6 for r in report)
    assert all(
        math.isnan(r["precision"]) for r in report
    )  # No frozen deployment threshold supplied.


def test_variant_effects_are_explicit_scenarios(tmp_path):
    cfg = create_example(tmp_path / "inputs", "genome_only")
    p = tmp_path / "effects.yaml"
    p.write_text(yaml.safe_dump(dict(run_config=str(cfg), variants=str(cfg.parent / "sample.vcf"))))
    rows = variant_effects_command(p, tmp_path / "effects")
    assert len(rows) == 6  # Variant lies inside three model input windows, two assay heads each.
    assert all(r["scenario"] == "single_variant_on_reference_context" for r in rows)
    assert all(math.isnan(r["full_delta_pace"]) for r in rows)
    assert any(abs(r["delta_signal"]) > 0 for r in rows)


def test_cli_help_and_errors(tmp_path, capsys):
    cfg = create_example(tmp_path / "inputs", "measured")
    assert main(["validate", "--config", str(cfg)]) == 0
    assert main(["run", "--config", str(cfg), "--out", str(tmp_path / "run")]) == 0
    assert main(["run", "--config", str(cfg), "--out", str(tmp_path / "run")]) == 2
    assert "already exists" in capsys.readouterr().err
