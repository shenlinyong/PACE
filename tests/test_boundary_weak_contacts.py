"""Numerical identities, raw-count provenance and end-to-end weak calibration."""

import copy
import json
import math

import numpy as np
import pytest
import yaml

from pace_livestock.boundary_prior import (
    BoundaryIndex,
    contact_prior,
    fit_contact_grid,
    fit_kappa,
    load_boundaries,
    motif_boundaries,
    posterior_contact,
    scan_motifs,
    unique_count_pairs,
)
from pace_livestock.cli import main
from pace_livestock.config import load_config
from pace_livestock.core import score
from pace_livestock.demo import create_example
from pace_livestock.errors import PaceError
from pace_livestock.evidence.assets import load_asset
from pace_livestock.evidence.contact import distance_prior
from pace_livestock.evidence.resolve import resolve_contacts
from pace_livestock.io.cooler import query_contacts
from pace_livestock.io.tables import read_table, write_table
from pace_livestock.pipeline import compute, run
from pace_livestock.provenance import file_hash, write_json
from pace_livestock.schemas import load_tables
from pace_livestock.weak_labels import (
    calibrate,
    fit_chromosome_grid,
    map_weak_labels,
    soft_average_precision,
)


def asset():
    return dict(
        model_id="toy_boundary",
        kind="contact_prior",
        is_synthetic=True,
        species="synthetic",
        assembly="toy_assembly",
        context_id="toy_tissue",
        target_level="individual",
        training_sources=[],
        calibration_sources=[],
        test_sources=[],
        validation={},
        a=1.0,
        gamma=1.0,
        beta=0.0,
        d_ref=1000.0,
        d_min=1000.0,
        scale="toy_contact",
        resolution=500,
        normalization_id="toy_mean",
        balancing="unbalanced",
        window_id="bin_pair",
        kappa=2.0,
    )


@pytest.fixture
def configured(tmp_path):
    cfg = load_config(create_example(tmp_path / "input", "measured"))
    prior = tmp_path / "prior.json"
    write_json(prior, asset())
    cfg["contact"].update(
        prior_path=str(prior),
        mode="shrinkage",
        reliability="per_pair",
        reliability_source="gamma_poisson",
        normalization_id="toy_mean",
        balancing="unbalanced",
        window_id="bin_pair",
    )
    rows = read_table(cfg["inputs"]["observed_contacts"])
    for row in rows:
        row.update(raw_count=int(row["contact_value"]), count_to_contact=1.0)
    write_table(cfg["inputs"]["observed_contacts"], rows)
    return cfg


def write_config(path, cfg):
    path.write_text(yaml.safe_dump(cfg))
    return str(path)


def test_boundary_sum_endpoints_and_symmetry():
    index = BoundaryIndex(
        [dict(chrom="1", position0=x, strength=s) for x, s in [(10, 0.5), (20, 1), (30, 2)]]
    )
    assert index.strength("1", 10, 30) == 1
    assert index.strength("1", 31, 9) == 3.5
    assert index.strength("2", 0, 100) == 0
    assert index.strength("1", 20, 20) == 0


def test_beta_zero_is_exact_powerlaw():
    prior = asset()
    for d in (0, 25, 1000, 30000):
        assert contact_prior("1", 0, d, prior) == distance_prior(d, prior)
    index = BoundaryIndex([dict(chrom="1", position0=500, strength=0.5)])
    assert contact_prior("1", 0, 2000, {**prior, "beta": 2}, index) == pytest.approx(0.5 / math.e)


@pytest.mark.parametrize("strength", [-1, float("nan"), True])
def test_boundary_rejects_invalid_strength(strength):
    with pytest.raises(PaceError):
        BoundaryIndex([dict(chrom="1", position0=0, strength=strength)])


def test_boundary_rejects_duplicate_position():
    row = dict(chrom="1", position0=10, strength=1)
    with pytest.raises(PaceError, match="Duplicate"):
        BoundaryIndex([row, row])


def test_boundary_file_is_checksum_bound(tmp_path):
    path = tmp_path / "b.tsv"
    write_table(path, [dict(chrom="1", position0=500, strength=1)])
    prior = {
        **asset(),
        "beta": 1,
        "asset_directory": str(tmp_path),
        "boundary_table": "b.tsv",
        "boundary_sha256": file_hash(path),
    }
    assert load_boundaries(prior).strength("1", 0, 1000) == 1
    path.write_text(path.read_text().replace("500", "600"))
    with pytest.raises(PaceError, match="checksum"):
        load_boundaries(prior)


def test_divergent_motifs_and_chip_occupancy():
    motifs = [
        dict(chrom="1", start=x, end=x + 10, strand=s, strength=0.8)
        for x, s in [(10, "+"), (30, "-"), (60, "+"), (90, "-")]
    ]
    result = motif_boundaries(motifs)
    assert len(result) == 1 and result[0]["position0"] == 50
    assert result[0]["boundary_status"] == "candidate_unvalidated"
    assert motif_boundaries(motifs, chip=[dict(chrom="1", start=20, end=50)]) == []
    assert (
        motif_boundaries(motifs, chip=[dict(chrom="1", start=0, end=100)])[0]["evidence_type"]
        == "chip_supported_motif"
    )
    assert motif_boundaries(motifs, max_gap_bp=10) == []


def test_pwm_both_strands_and_chunk_overlap(tmp_path):
    path = tmp_path / "ref.fa"
    path.write_text(">1\nAACGNN\nCGTAAA\n")
    pwm = [dict(zip("ACGT", r, strict=True)) for r in [(10, 0, 0, 0), (0, 10, 0, 0), (0, 0, 10, 0)]]
    small = scan_motifs(path, pwm, threshold=5, chunk_size=2)
    assert small == scan_motifs(path, pwm, threshold=5, chunk_size=100)
    assert {(r["start"], r["strand"]) for r in small} == {(1, "+"), (6, "-")}


def test_gamma_poisson_matches_posterior_rate_identity():
    # Prior Gamma(2, rate=.5), exposure=2 => posterior Gamma(8, rate=2.5).
    result = posterior_contact(6, 0.5, 4, 2)
    assert result["resolved_value"] == pytest.approx(8 / 2.5)
    assert result["posterior_variance"] == pytest.approx(8 / 2.5**2)
    assert result["reliability"] == 0.8
    assert posterior_contact(0, 0.5, 4, 2)["resolved_value"] == pytest.approx(0.8)
    assert posterior_contact(6, 0.05, 4, 2)["reliability"] > result["reliability"]


@pytest.mark.parametrize(
    "args",
    [
        (1.2, 1, 1, 2),
        (-1, 1, 1, 2),
        (1, 0, 1, 2),
        (1, 1, 0, 2),
        (1, 1, 1, 0),
        (1, 1, 1, float("inf")),
    ],
)
def test_invalid_poisson_inputs(args):
    with pytest.raises(PaceError):
        posterior_contact(*args)


def test_kappa_moment_equation():
    # Numerator= (0-2)^2 + ((2-2)^2-2) + ((6-2)^2-6) = 12; denominator=12.
    assert fit_kappa([0, 2, 6], [2, 2, 2])["kappa"] == pytest.approx(1)
    assert fit_kappa([2, 2, 2], [2, 2, 2])["at_bound"]
    with pytest.raises(PaceError):
        fit_kappa([2], [2])


def count_rows():
    return [
        dict(
            chrom="1",
            anchor0=0,
            tss0=d,
            sample_id="s",
            bin_pair_id=f"0:{d // 1000}",
            raw_count=y,
            count_to_contact=0.01,
            contact_value=y * 0.01,
            resolution=1000,
            measurement_status="observed",
        )
        for d, y in [(1000, 100), (2000, 50), (4000, 25), (8000, 13)]
    ]


def test_count_pairs_deduplicate_and_reject_conflicts():
    rows = count_rows()
    index = BoundaryIndex([])
    unique = unique_count_pairs(rows + rows, resolution=1000, boundaries=index)
    assert len(unique) == 4
    with pytest.raises(PaceError, match="reproduce"):
        unique_count_pairs([{**rows[0], "raw_count": 2}], resolution=1000, boundaries=index)
    with pytest.raises(PaceError, match="disagree"):
        unique_count_pairs(
            rows + [{**rows[0], "raw_count": 2, "contact_value": 0.02}],
            resolution=1000,
            boundaries=index,
        )


def test_profile_poisson_fit_recovers_decay():
    rows = unique_count_pairs(count_rows(), resolution=1000, boundaries=BoundaryIndex([]))
    fit = fit_contact_grid(rows, gamma_grid=[0.5, 1, 1.5], beta_grid=[0], d_ref=1000, d_min=1000)
    assert fit["gamma"] == 1 and fit["beta"] == 0
    assert fit["a"] == pytest.approx(1, abs=0.01)
    assert fit["n_unique_pairs"] == 4


def test_pipeline_per_pair_with_masked_bins_and_true_zero(configured):
    rows = read_table(configured["inputs"]["observed_contacts"])
    rows[0].update(raw_count=0, contact_value=0)
    rows[1].update(contact_value=None, measurement_status="unmappable", count_to_contact=None)
    write_table(configured["inputs"]["observed_contacts"], rows)
    result = compute(configured)
    contacts = {(r["element_id"], r["promoter_id"]): r for r in result["resolved_contacts"]}
    zero = contacts["E1", "P1"]
    assert (
        zero["observed_value"] == 0
        and zero["resolved_value"] > 0
        and zero["evidence_type"] == "fused"
    )
    masked = contacts["E2", "P1"]
    assert masked["reliability"] == 0 and masked["reason"] == "unavailable_bin_prior"
    assert masked["observation_sample_id"] is None
    assert all(math.isfinite(r["pace_score"]) for r in result["scores"])


def test_pipeline_raw_count_required(configured):
    rows = read_table(configured["inputs"]["observed_contacts"])
    rows[0].pop("raw_count")
    write_table(configured["inputs"]["observed_contacts"], rows)
    with pytest.raises(PaceError, match="raw_count"):
        compute(configured)


def test_auto_kappa_uses_physical_pairs_once(configured):
    prior = asset()
    prior.pop("kappa")
    write_json(configured["contact"]["prior_path"], prior)
    result = compute(configured)
    assert {r["kappa_source"] for r in result["resolved_contacts"]} == {"run_unique_bin_pairs"}
    tables = load_tables(configured)
    loaded = load_asset(configured["contact"]["prior_path"], configured, kind="contact_prior")
    base = resolve_contacts(tables, configured, prior_asset=loaded)
    # A duplicate alternative TSS queries the same bins but is not an extra replicate.
    tables["promoters"].append({**tables["promoters"][0], "promoter_id": "P3"})
    tables["observed_contacts"] += [
        {**r, "promoter_id": "P3"} for r in tables["observed_contacts"] if r["promoter_id"] == "P1"
    ]
    duplicate = resolve_contacts(tables, configured, prior_asset=loaded)
    assert base[0]["kappa"] == duplicate[0]["kappa"]


def test_cooler_preserves_raw_counts_and_weights(tmp_path):
    cooler = pytest.importorskip("cooler")
    pd = pytest.importorskip("pandas")
    path = tmp_path / "raw.cool"
    cooler.create_cooler(
        str(path),
        pd.DataFrame(
            dict(
                chrom=["1"] * 3,
                start=[0, 1000, 2000],
                end=[1000, 2000, 3000],
                weight=[0.5, 0.2, np.nan],
            )
        ),
        pd.DataFrame(dict(bin1_id=[0, 0], bin2_id=[1, 2], count=[20, 7])),
    )
    pairs = [
        dict(element_id="E", promoter_id=str(t), chrom="1", anchor0=10, tss0=t)
        for t in (1100, 2100)
    ]
    rows = query_contacts(path, pairs, resolution=1000, balanced=True, missing_pixels_are_zero=True)
    assert rows[0]["raw_count"] == 20 and rows[0]["count_to_contact"] == pytest.approx(0.1)
    assert rows[0]["contact_value"] == pytest.approx(2)
    assert rows[1]["raw_count"] == 7 and rows[1]["measurement_status"] == "unmappable"
    assert math.isnan(rows[1]["count_to_contact"])


def label_inputs():
    units = [
        dict(element_id="E", chrom="1", start=10, end=20),
        dict(element_id="F", chrom="1", start=20, end=30),
    ]
    candidates = [dict(element_id=e, gene_id="G") for e in ("E", "F")]
    variants = [
        dict(variant_id=str(i), chrom="1", pos0=p, gene_id="G", pip=q, signal_id="S")
        for i, p, q in [(1, 10, 0.2), (2, 19, 0.3)]
    ]
    return variants, units, candidates


def test_pip_independence_and_signal_aware_sum():
    variants, units, candidates = label_inputs()
    independent = map_weak_labels(variants, units, candidates)
    assert independent[0]["weak_label"] == pytest.approx(0.44)
    assert independent[1]["weak_label"] == 0
    assert independent[1]["label_status"] == "unlabelled_background"
    signals = map_weak_labels(variants, units, candidates, aggregation="independent_signals")
    assert signals[0]["weak_label"] == 0.5


@pytest.mark.parametrize("problem", ["duplicate", "pip", "signal_mass", "chromosome"])
def test_bad_fine_mapping_rejected(problem):
    variants, units, candidates = label_inputs()
    if problem == "duplicate":
        variants.append(variants[0])
    elif problem == "pip":
        variants[0]["pip"] = 1.1
    elif problem == "signal_mass":
        variants[0]["pip"] = 0.9
    else:
        variants[0]["chrom"] = "2"
    with pytest.raises(PaceError):
        map_weak_labels(variants, units, candidates, aggregation="independent_signals")


def test_soft_ap_is_tie_aware_and_reduces_to_binary_ap():
    assert soft_average_precision([1, 0, 1], [3, 2, 1]) == pytest.approx(5 / 6)
    assert soft_average_precision([1, 0], [1, 1]) == 0.5
    assert soft_average_precision([0.2, 0.8], [1, 1]) == 0.5
    assert soft_average_precision([0.8, 0.2], [1, 1]) == 0.5


def test_nested_selection_never_uses_final_test_labels():
    labels = [dict(chrom=c, weak_label=q) for c in ("1", "2", "3", "4") for q in (1, 0)]
    predictions = {(1, 0, 0): [0.2, 0.8] * 4, (1, 0, 1): [0.8, 0.2] * 4}
    a = fit_chromosome_grid(labels, predictions, test_chromosomes=["4"])
    changed = copy.deepcopy(labels)
    changed[-2]["weak_label"], changed[-1]["weak_label"] = 0, 1
    b = fit_chromosome_grid(changed, predictions, test_chromosomes=["4"])
    assert a["selected"] == b["selected"] == dict(gamma=1, beta=0, eta=1)
    assert a["nested_folds"] == b["nested_folds"]
    assert a["test_soft_ap"] != b["test_soft_ap"]
    assert all(f["test_chromosome"] not in f["training_chromosomes"] for f in a["nested_folds"])


def test_ties_leave_B_off_and_too_few_chromosomes_fail():
    labels = [dict(chrom=c, weak_label=q) for c in ("1", "2", "3", "4") for q in (1, 0)]
    result = fit_chromosome_grid(
        labels, {(1, 0, 0): [1, 0] * 4, (1, 1, 1): [1, 0] * 4}, test_chromosomes=["4"]
    )
    assert result["selected"]["eta"] == 0 and result["selected"]["beta"] == 0
    with pytest.raises(PaceError, match="three"):
        fit_chromosome_grid(labels[:6], {(1, 0, 0): [1, 0] * 3}, test_chromosomes=["3"])


def four_chromosomes(cfg):
    units = read_table(cfg["inputs"]["units"])
    promoters = read_table(cfg["inputs"]["promoters"])
    candidates = read_table(cfg["inputs"]["candidates"])
    activity = read_table(cfg["inputs"]["observed_activity"])
    batches = {k: [] for k in ("units", "promoters", "candidates", "observed_activity")}
    variants = []
    for chrom in ("c1", "c2", "c3", "c4"):
        for row, start in zip(units, [5000, 20000, 30000], strict=True):
            batches["units"].append(
                {
                    **row,
                    "element_id": chrom + row["element_id"],
                    "chrom": chrom,
                    "start": start,
                    "end": start + 500,
                    "anchor0": start + 249,
                }
            )
        for row, tss in zip(promoters, [10000, 40000], strict=True):
            batches["promoters"].append(
                {
                    **row,
                    "gene_id": chrom + row["gene_id"],
                    "promoter_id": chrom + row["promoter_id"],
                    "chrom": chrom,
                    "tss0": tss,
                }
            )
        batches["candidates"] += [
            {**r, "element_id": chrom + r["element_id"], "gene_id": chrom + r["gene_id"]}
            for r in candidates
        ]
        batches["observed_activity"] += [
            {
                **r,
                "element_id": chrom + r["element_id"],
                "signal": {"E1": 1.0, "E2": 2.8, "E3": 1.0}[r["element_id"]],
            }
            for r in activity
        ]
        variants += [
            dict(
                variant_id=chrom + gene,
                chrom=chrom,
                pos0=pos,
                gene_id=chrom + gene,
                pip=0.9,
                signal_id="signal",
            )
            for gene, pos in [("G1", 5000), ("G2", 30000)]
        ]
    for key, rows in batches.items():
        write_table(cfg["inputs"][key], rows)
    cfg["inputs"]["observed_contacts"] = None
    cfg["contact"].update(mode="prior_only", reliability=None, reliability_source=None)
    return variants


def test_learned_eta_reaches_main_scores_and_reusable_artifact(configured, tmp_path):
    variants = four_chromosomes(configured)
    model, labels = calibrate(
        configured,
        variants,
        gamma_grid=[1],
        beta_grid=[0],
        eta_grid=[0, 1],
        test_chromosomes=["c4"],
        aggregation="independent_signals",
    )
    assert model["selected"]["eta"] == 1
    assert len(labels) == 24 and model["functional_validation"] is False
    base = compute(configured)
    path = tmp_path / "weak_model.json"
    write_json(path, model)
    configured["allocation"]["weak_model_path"] = str(path)
    result = run(configured, tmp_path / "weak_run")
    expected, _ = score(base["scores"], eta=1)
    assert [r["pace_score"] for r in result["scores"]] == pytest.approx(
        [r["pace_score"] for r in expected]
    )
    assert all(
        math.isfinite(r["B"]) and r["allocation_evidence"] == "eqtl_weak" for r in result["scores"]
    )
    assert result["eta_calibration"]["reuse"] is True
    configured["activity"]["pseudocounts"] = {"ATAC": 1}
    with pytest.raises(PaceError, match="differs"):
        compute(configured)


def test_all_five_cli_commands_and_main_run(configured, tmp_path):
    motifs = tmp_path / "motifs.tsv"
    write_table(
        motifs,
        [
            dict(chrom="chrToy", start=x, end=x + 10, strand=s, strength=0.8)
            for x, s in [(6900, "-"), (7100, "+")]
        ],
    )
    config = tmp_path / "command.yaml"
    assert (
        main(
            [
                "boundaries",
                "--config",
                write_config(config, dict(motifs=str(motifs))),
                "--out",
                str(tmp_path / "boundaries"),
            ]
        )
        == 0
    )
    prior_cfg = {
        k: v
        for k, v in asset().items()
        if k
        in {
            "model_id",
            "is_synthetic",
            "species",
            "assembly",
            "context_id",
            "target_level",
            "scale",
            "resolution",
            "normalization_id",
            "balancing",
            "window_id",
            "a",
            "gamma",
            "beta",
            "d_ref",
            "d_min",
            "kappa",
        }
    }
    prior_cfg.update(beta=0.5, boundaries=str(tmp_path / "boundaries/boundaries.tsv"))
    assert (
        main(
            [
                "prior",
                "--config",
                write_config(config, prior_cfg),
                "--out",
                str(tmp_path / "prior_asset"),
            ]
        )
        == 0
    )
    configured["contact"]["prior_path"] = str(tmp_path / "prior_asset")
    assert (
        main(
            ["fuse", "--config", write_config(config, configured), "--out", str(tmp_path / "fused")]
        )
        == 0
    )
    # Fit existing raw-contact output with measurement metadata and fixed-bin coordinates.
    rows = count_rows()
    for row in rows:
        row.update(
            scale="toy_contact",
            normalization_id="toy_mean",
            balancing="unbalanced",
            window_id="bin_pair",
        )
    raw = tmp_path / "raw.tsv"
    write_table(raw, rows)
    fitting = {k: v for k, v in prior_cfg.items() if k not in ("a", "gamma", "beta", "kappa")}
    fitting.update(data=str(raw), resolution=1000, gamma_grid=[0.5, 1, 1.5], beta_grid=[0])
    assert (
        main(
            [
                "fit-hic",
                "--config",
                write_config(config, fitting),
                "--out",
                str(tmp_path / "fitted"),
            ]
        )
        == 0
    )
    configured["contact"]["prior_path"] = str(tmp_path / "prior.json")
    variants = four_chromosomes(configured)
    labels = tmp_path / "eqtl.tsv"
    write_table(labels, variants)
    run_cfg = write_config(tmp_path / "run.yaml", configured)
    weak = dict(
        run_config=run_cfg,
        labels=str(labels),
        **configured["context"],
        gamma_grid=[1],
        beta_grid=[0],
        eta_grid=[0, 1],
        aggregation="independent_signals",
        test_chromosomes=["c4"],
    )
    assert (
        main(
            ["fit-labels", "--config", write_config(config, weak), "--out", str(tmp_path / "weak")]
        )
        == 0
    )
    configured["allocation"]["weak_model_path"] = str(tmp_path / "weak/weak_model.json")
    assert (
        main(
            ["run", "--config", write_config(config, configured), "--out", str(tmp_path / "final")]
        )
        == 0
    )
    assert json.loads((tmp_path / "final/eta_calibration.json").read_text())["eta"] == 1


def test_bad_configuration_is_rejected(configured, tmp_path):
    configured["contact"]["mode"] = "observed"
    with pytest.raises(PaceError, match="shrinkage"):
        load_config(write_config(tmp_path / "bad.yaml", configured))


def test_joint_decay_and_boundary_fit():
    rows = [
        dict(raw_count=y, count_to_contact=0.01, distance_bp=d, boundary_strength=b)
        for b, counts in [(0, [100, 50, 25]), (1, [50, 25, 12])]
        for d, y in zip([1000, 2000, 4000], counts, strict=True)
    ]
    fitted = fit_contact_grid(
        rows, gamma_grid=[0.5, 1, 1.5], beta_grid=[0, math.log(2), 1.4], d_ref=1000, d_min=1000
    )
    assert fitted["gamma"] == 1
    assert fitted["beta"] == math.log(2)


def test_resolved_per_pair_reuse_verifies_raw_provenance(configured, tmp_path):
    original = run(configured, tmp_path / "original")
    configured["inputs"]["resolved_contacts"] = str(tmp_path / "original/resolved_contacts.tsv")
    configured["inputs"]["evidence"] = str(tmp_path / "original/evidence.tsv")
    configured["inputs"]["sources"] = str(tmp_path / "original/sources.tsv")
    repeated = compute(configured)
    assert [r["pace_score"] for r in repeated["scores"]] == [
        r["pace_score"] for r in original["scores"]
    ]
    rows = read_table(configured["inputs"]["resolved_contacts"])
    rows[0]["kappa"] = 7
    write_table(configured["inputs"]["resolved_contacts"], rows)
    with pytest.raises(PaceError, match="recomputation"):
        compute(configured)


def test_benchmark_keeps_frozen_contact_and_varies_eta(configured, tmp_path):
    from pace_livestock.evaluation.benchmark import benchmark_command

    variants = four_chromosomes(configured)
    model, _ = calibrate(
        configured,
        variants,
        gamma_grid=[1],
        beta_grid=[0],
        eta_grid=[0, 1],
        test_chromosomes=["c4"],
        aggregation="independent_signals",
    )
    path = tmp_path / "weak.json"
    write_json(path, model)
    configured["allocation"]["weak_model_path"] = str(path)
    run_path = write_config(tmp_path / "run.yaml", configured)
    labels = [
        dict(
            label_id=f"L{i}",
            assayed_region_id=e,
            gene_id=g,
            context_id="toy_tissue",
            effect_direction="down" if positive else "none",
            label_status="enhancing_positive" if positive else "powered_negative",
        )
        for i, (e, g, positive) in enumerate(
            [
                (e, g, (e, g) in (("c4E1", "c4G1"), ("c4E3", "c4G2")))
                for e in ("c4E1", "c4E2", "c4E3")
                for g in ("c4G1", "c4G2")
            ]
        )
    ]
    write_table(tmp_path / "labels.tsv", labels)
    write_table(
        tmp_path / "membership.tsv",
        [dict(region_id=e, element_id=e) for e in ("c4E1", "c4E2", "c4E3")],
    )
    command = write_config(
        tmp_path / "benchmark.yaml",
        dict(run_config=run_path, labels="labels.tsv", region_membership="membership.tsv"),
    )
    report = benchmark_command(command, tmp_path / "benchmark")
    metrics = {r["method"]: r for r in report if r["scope"] == "all_tested_labels"}
    assert (
        metrics["PACE_eta1"]["average_precision"] == metrics["PACE_eqtl_weak"]["average_precision"]
    )
    assert metrics["PACE_eta0"]["average_precision"] < metrics["PACE_eta1"]["average_precision"]


def test_published_contact_workflow(tmp_path):
    import shutil
    from pathlib import Path

    example = Path(__file__).resolve().parents[1] / "examples/contact"
    root = tmp_path / "demo"
    shutil.copytree(example, root)
    for command, config, output in [
        ("boundaries", "boundaries.yaml", "boundary_asset"),
        ("prior", "prior.yaml", "prior_asset"),
        ("fit-hic", "fit-hic.yaml", "fitted_prior"),
        ("fuse", "fuse.yaml", "fused"),
        ("fit-labels", "fit-labels.yaml", "weak_fit"),
        ("run", "calibrated.yaml", "calibrated"),
    ]:
        assert main([command, "--config", str(root / config), "--out", str(root / output)]) == 0
    manifest = json.loads((root / "calibrated/eta_calibration.json").read_text())
    assert manifest["status"] == "weak_fitted" and not manifest["functional_validation"]
