"""Regression checks for invalid annotations, replicate pooling and evidence closure."""

import math

import pytest
import yaml

from pace_livestock.config import load_config
from pace_livestock.demo import create_example
from pace_livestock.errors import PaceError
from pace_livestock.io.bed_gtf import read_bed
from pace_livestock.io.intervals import interval_features
from pace_livestock.io.methylation import (
    load_reference_cpg,
    promoter_methylation_regions,
    summarize_methylation,
)
from pace_livestock.io.rna import transcript_tpm_to_gene
from pace_livestock.io.tables import read_table, write_table
from pace_livestock.operations import prepare_command
from pace_livestock.pipeline import compute


def omics_config(tmp_path):
    cfg = load_config(create_example(tmp_path / "inputs", "measured"))
    samples = read_table(cfg["inputs"]["samples"])
    for assay in ("RNA", "H3K4me1", "H3K9me3", "WGBS"):
        samples.append({**samples[0], "sample_id": f"S_{assay}", "assay": assay})
    write_table(cfg["inputs"]["samples"], samples)
    observed = read_table(cfg["inputs"]["observed_activity"])
    for assay in ("H3K4me1", "H3K9me3"):
        observed.append({**observed[0], "sample_id": f"S_{assay}", "assay": assay})
    write_table(cfg["inputs"]["observed_activity"], observed)
    expression = tmp_path / "expression.tsv"
    write_table(
        expression, [{"gene_id": "G1", "sample_id": "S_RNA", "tpm": 123, "status": "invalid"}]
    )
    cfg["inputs"]["expression"] = str(expression)
    methylation = tmp_path / "methylation.tsv"
    write_table(
        methylation,
        [
            dict(
                chrom="chrToy",
                dyad_start0=5200,
                methylated_count=1,
                total_count=2,
                sample_id="S_WGBS",
                assay="WGBS",
            )
        ],
    )
    cfg["inputs"]["methylation"] = str(methylation)
    return cfg


def test_failed_rna_is_missing_and_all_feature_evidence_resolves(tmp_path):
    result = compute(omics_config(tmp_path))
    rna = next(r for r in result["features"] if r["feature_name"] == "RNA:TPM")
    assert math.isnan(rna["value"])
    assert rna["status"] == "invalid"
    evidence = {r["evidence_id"]: r for r in result["evidence"]}
    sources = {r["source_id"] for r in result["sources"]}
    assert all(r["evidence_id"] in evidence for r in result["features"])
    assert all(r["source_id"] in sources for r in evidence.values())
    assert all(
        parent in evidence
        for r in evidence.values()
        for parent in (r.get("parent_evidence_ids") or "").split(";")
        if parent
    )
    assert result["qc"]["multiomics_roles"]["H3K9me3"] == "annotation_only"


def test_rna_donor_aggregation_is_order_independent_and_not_sample_count_weighted(tmp_path):
    cfg = omics_config(tmp_path)
    cfg["target_level"] = "population_mean"
    samples = read_table(cfg["inputs"]["samples"])
    template = next(r for r in samples if r["assay"] == "RNA")
    rows = []
    for i, (donor, value, status) in enumerate(
        [
            ("donor_A", 0, "observed"),
            ("donor_A", 2, "observed"),
            ("donor_B", 9, "observed"),
            ("donor_B", 1000, "invalid"),
        ]
    ):
        sample = f"RNA_{i}"
        samples.append(
            {**template, "sample_id": sample, "donor_id": donor, "technical_replicate": str(i)}
        )
        rows.append(dict(gene_id="G1", sample_id=sample, tpm=value, status=status))
    write_table(cfg["inputs"]["samples"], samples)
    actual = []
    for ordered in (rows, rows[::-1]):
        write_table(cfg["inputs"]["expression"], ordered)
        result = compute(cfg)
        features = [r for r in result["features"] if r["feature_name"] == "RNA:TPM"]
        assert len(features) == 1
        feature = features[0]
        assert feature["n_observed_donors"] == 2
        assert feature["n_observed_samples"] == 3
        actual.append(feature["value"])
    # mean(mean(0,2),9) = 5; neither mean(0,2,9) nor the final input row.
    assert actual == [5, 5]


def test_methylation_main_pipeline_uses_coverage_and_promoter_windows(tmp_path):
    cfg = omics_config(tmp_path)
    promoters = read_table(cfg["inputs"]["promoters"])
    p = promoters[0]
    write_table(
        cfg["inputs"]["methylation"],
        [
            dict(
                chrom=p["chrom"],
                dyad_start0=int(p["tss0"]),
                methylated_count=0,
                total_count=2,
                sample_id="S_WGBS",
                assay="WGBS",
            ),
            dict(
                chrom=p["chrom"],
                dyad_start0=int(p["tss0"]) + 1,
                methylated_count=8,
                total_count=8,
                sample_id="S_WGBS",
                assay="WGBS",
            ),
        ],
    )
    reference = tmp_path / "reference_cpg.tsv"
    write_table(reference, [dict(entity_type="promoter", entity_id=p["promoter_id"], n_cpg=4)])
    cfg["methylation"] = dict(
        minimum_coverage=5,
        promoter_upstream_bp=2,
        promoter_downstream_bp=2,
        reference_cpg_path=str(reference),
    )
    result = compute(cfg)
    features = {
        r["feature_name"]: r
        for r in result["features"]
        if r["entity_type"] == "promoter" and r["entity_id"] == p["promoter_id"]
    }
    assert features["DNA_methylation:M_site"]["value"] == 1
    assert features["DNA_methylation:cpg_coverage_fraction"]["value"] == 0.25
    assert any(r["feature_name"] == "promoter:DNA_methylation:M_site" for r in result["features"])


def test_unequal_technical_replicates_do_not_reweight_biological_replicates(tmp_path):
    cfg = omics_config(tmp_path)
    samples = read_table(cfg["inputs"]["samples"])
    template = next(r for r in samples if r["assay"] == "RNA")
    rows = []
    for i, (biological_replicate, value) in enumerate(
        [
            ("bio_A", 0),
            ("bio_A", 0),
            ("bio_A", 0),
            ("bio_B", 100),
        ]
    ):
        sample_id = f"RNA_bio_{i}"
        samples.append(
            {
                **template,
                "sample_id": sample_id,
                "biological_replicate": biological_replicate,
                "technical_replicate": str(i),
            }
        )
        rows.append(dict(gene_id="G1", sample_id=sample_id, tpm=value, status="observed"))
    write_table(cfg["inputs"]["samples"], samples)
    write_table(cfg["inputs"]["expression"], rows)
    result = compute(cfg)
    feature = next(r for r in result["features"] if r["feature_name"] == "RNA:TPM")
    # The three technical repeats estimate one biological replicate's mean.
    assert feature["value"] == 50
    assert feature["n_observed_samples"] == 4
    assert feature["n_observed_biological_replicates"] == 2
    assert feature["n_observed_donors"] == 1


def test_promoter_windows_are_strand_aware_and_reference_denominators_explicit(tmp_path):
    rows = [
        dict(promoter_id="P+", chrom="1", tss0=100, strand="+"),
        dict(promoter_id="P-", chrom="1", tss0=100, strand="-"),
    ]
    windows = promoter_methylation_regions(rows, upstream_bp=20, downstream_bp=5)
    assert [(r["start"], r["end"]) for r in windows] == [(80, 106), (95, 121)]
    reference = tmp_path / "cpg.tsv"
    write_table(reference, [dict(entity_type="promoter", entity_id="P+", n_cpg=0)])
    assert load_reference_cpg(reference) == {("promoter", "P+"): 0}
    with pytest.raises(PaceError, match="minimum_coverage"):
        summarize_methylation([], [], minimum_coverage=0)


def test_invalid_transcript_cannot_create_a_partial_gene_tpm():
    rows = [
        dict(transcript_id="T1", sample_id="S", tpm=2, status="observed"),
        dict(transcript_id="T2", sample_id="S", tpm=999, status="invalid"),
    ]
    mapping = [dict(transcript_id=t, gene_id="G", promoter_id="P") for t in ("T1", "T2")]
    result = transcript_tpm_to_gene(rows, mapping)[0]
    assert math.isnan(result["tpm"]) and result["status"] == "unmeasured"


def test_bed6_motif_orientation_survives_file_adapter(tmp_path):
    path = tmp_path / "ctcf.bed"
    path.write_text("1\t10\t20\tmotif1\t1\t+\n1\t40\t50\tmotif2\t1\t-\n")
    peaks = read_bed(path, source_id="S")
    rows = interval_features(
        [dict(element_id="E", chrom="1", start=0, end=100)],
        peaks,
        evidence_id="ev",
        feature_prefix="CTCF",
        motif_strands=True,
    )
    values = {r["feature_name"]: r["value"] for r in rows}
    assert values["CTCF:plus_motifs"] == values["CTCF:minus_motifs"] == 1


def test_ctcf_motif_configuration_is_used_by_prepare_command(tmp_path):
    (tmp_path / "ctcf.bed").write_text("1\t10\t20\tm1\t1\t+\n1\t40\t50\tm2\t1\t-\n")
    write_table(tmp_path / "units.tsv", [dict(element_id="E", chrom="1", start=0, end=100)])
    path = tmp_path / "prepare.yaml"
    path.write_text(
        yaml.safe_dump(
            dict(
                kind="bed_features",
                bed="ctcf.bed",
                units="units.tsv",
                source_id="S",
                evidence_id="ev",
                feature_prefix="CTCF",
                motif_strands=True,
            )
        )
    )
    prepare_command(path, tmp_path / "prepared")
    values = {
        r["feature_name"]: float(r["value"])
        for r in read_table(tmp_path / "prepared" / "features.tsv")
    }
    assert values["CTCF:plus_motifs"] == values["CTCF:minus_motifs"] == 1


@pytest.mark.parametrize(
    "change,pattern",
    [
        ({"status": "wrong"}, "invalid status"),
        ({"status": "observed", "tpm": None}, "finite TPM"),
        ({"gene_id": "missing_gene"}, "unknown gene_id"),
    ],
)
def test_expression_schema_rejects_inconsistent_annotations(tmp_path, change, pattern):
    cfg = omics_config(tmp_path)
    rows = read_table(cfg["inputs"]["expression"])
    rows[0].update(change)
    write_table(cfg["inputs"]["expression"], rows)
    with pytest.raises(PaceError, match=pattern):
        compute(cfg)
