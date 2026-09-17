"""Catalog and real small-file format adapters; heavy IO tests require the io extra."""

import math

import numpy as np
import pytest

from pace_livestock.catalog import candidate_edges, canonical_units, map_labels
from pace_livestock.errors import PaceError
from pace_livestock.io.bed_gtf import read_gtf
from pace_livestock.io.bigwig import quantify_bigwig
from pace_livestock.io.cooler import query_contacts
from pace_livestock.io.methylation import merge_stranded_cpg
from pace_livestock.io.rna import transcript_tpm_to_gene
from pace_livestock.io.variants import Reference


def test_catalog_duplicate_regions_and_bounds():
    regions = [dict(chrom="1", start=450, end=1100, region_id="R", source_id="S")]
    u, m, x = canonical_units(regions * 2, {"1": 1200}, include_promoters=False)
    assert [(r["start"], r["end"]) for r in u] == [(0, 500), (500, 1000)]
    assert len(m) == 2
    assert x[0]["reason"] == "incomplete_unit"
    reverse = canonical_units(regions[::-1], {"1": 1200}, include_promoters=False)[0]
    assert u == reverse


def test_sparse_candidates_match_independent_brute_force():
    rng = np.random.default_rng(17)
    units = [
        dict(element_id=f"E{i}", chrom="1", anchor0=int(x))
        for i, x in enumerate(rng.integers(0, 10000, 100))
    ]
    promoters = [
        dict(gene_id=f"G{i // 2}", chrom="1", tss0=int(x))
        for i, x in enumerate(rng.integers(0, 10000, 20))
    ]
    actual = {
        (r["element_id"], r["gene_id"]) for r in candidate_edges(units, promoters, radius=200)
    }
    # Independent direct predicate on this small fixture; production uses sparse interval search.
    expected = {
        (u["element_id"], p["gene_id"])
        for u in units
        for p in promoters
        if abs(u["anchor0"] - p["tss0"]) <= 200
    }
    assert actual == expected


def test_gtf_tss_strands_and_physical_dedup(tmp_path):
    p = tmp_path / "a.gtf"
    p.write_text(
        '1\ttoy\ttranscript\t101\t500\t.\t+\t.\tgene_id "G.1"; transcript_id "T1";\n'
        '1\ttoy\ttranscript\t101\t600\t.\t+\t.\tgene_id "G.1"; transcript_id "T2";\n'
        '1\ttoy\ttranscript\t201\t800\t.\t-\t.\tgene_id "H"; transcript_id "T3";\n'
    )
    promoters, mapping = read_gtf(p, aliases={"1": "chr1"})
    # GTF is one-based closed: + uses start-1=100; - uses end-1=799.
    assert [r["tss0"] for r in promoters] == [100, 799]
    assert len(mapping) == 3
    assert promoters[0]["gene_id"] == "G.1"


def test_label_regions_not_duplicated():
    membership = [{"region_id": "R", "element_id": e} for e in ("E1", "E2")]
    labels = [
        dict(
            assayed_region_id="R",
            gene_id="G",
            label_status="enhancing_positive",
            effect_direction="down",
        )
    ]
    usable, rejected = map_labels(labels, membership)
    assert not usable and rejected[0]["reason"] == "ambiguous_or_unmapped_region"


def test_real_bigwig_missing_and_negative(tmp_path):
    bw = pytest.importorskip(
        "pyBigWig", reason="install pace-livestock[io] to test real bigWig files"
    )
    p = tmp_path / "a.bw"
    with bw.open(str(p), "w") as w:
        w.addHeader([("chr1", 1000)])
        w.addEntries(["chr1"], [0], ends=[250], values=[4.0])
    units = [dict(element_id="E", chrom="chr1", start=0, end=500)]
    r = quantify_bigwig(p, units)[0]
    # Stored half-window has mean 4 and coverage .5; declaring unstored measured zero gives mean 2.
    assert r["signal"] == 4 and r["callable_fraction"] == 0.5
    assert quantify_bigwig(p, units, missing_is_measured_zero=True)[0]["signal"] == 2
    assert math.isnan(quantify_bigwig(p, units, minimum_callable_fraction=0.9)[0]["signal"])
    q = tmp_path / "negative.bw"
    with bw.open(str(q), "w") as w:
        w.addHeader([("chr1", 1000)])
        w.addEntries(["chr1"], [0], ends=[500], values=[-1.0])
    with pytest.raises(PaceError, match="Negative"):
        quantify_bigwig(q, units)


def test_real_cooler_sparse_shared_bin(tmp_path):
    cooler = pytest.importorskip("cooler", reason="install pace-livestock[io] to test cooler")
    pd = pytest.importorskip("pandas", reason="cooler dataframe dependency")
    p = tmp_path / "a.cool"
    bins = pd.DataFrame(
        {"chrom": ["chr1"] * 4, "start": [0, 500, 1000, 1500], "end": [500, 1000, 1500, 2000]}
    )
    pixels = pd.DataFrame({"bin1_id": [0, 0], "bin2_id": [1, 2], "count": [7, 3]})
    cooler.create_cooler(str(p), bins, pixels)
    pairs = [
        dict(element_id="E", promoter_id=f"P{i}", chrom="chr1", anchor0=100, tss0=t)
        for i, t in enumerate([550, 600, 1600])
    ]
    rows = query_contacts(p, pairs, resolution=500, balanced=False, missing_pixels_are_zero=True)
    assert rows[0]["contact_value"] == rows[1]["contact_value"] == 7
    assert rows[0]["bin_pair_id"] == rows[1]["bin_pair_id"]
    assert rows[2]["contact_value"] == 0
    assert math.isnan(
        query_contacts(p, pairs, resolution=500, balanced=False, missing_pixels_are_zero=False)[2][
            "contact_value"
        ]
    )
    with pytest.raises(PaceError, match="resolution"):
        query_contacts(p, pairs, resolution=1000, balanced=False, missing_pixels_are_zero=True)


def test_stranded_cpg_and_transcript_rna(tmp_path):
    p = tmp_path / "cpg.fa"
    p.write_text(">chr1\nAACGAA\n")
    rows = [
        dict(
            chrom="chr1",
            pos0=pos,
            strand=strand,
            methylated_count=1,
            total_count=2,
            sample_id="S",
            assay="WGBS",
        )
        for pos, strand in [(2, "+"), (3, "-")]
    ]
    with Reference(p) as reference:
        merged = merge_stranded_cpg(rows, reference)
    assert len(merged) == 1 and merged[0]["dyad_start0"] == 2
    assert merged[0]["total_count"] == 4
    rna = transcript_tpm_to_gene(
        [
            dict(transcript_id="T1", sample_id="S", tpm=2),
            dict(transcript_id="T2", sample_id="S", tpm=3),
        ],
        [dict(transcript_id=t, gene_id="G", promoter_id="P") for t in ("T1", "T2")],
    )
    assert rna[0]["tpm"] == 5


def test_interval_union_and_motif_orientation():
    from pace_livestock.io.intervals import interval_features

    region = dict(element_id="E", chrom="chr1", start=0, end=100)
    peaks = [
        dict(chrom="chr1", start=10, end=40, strand="+"),
        dict(chrom="chr1", start=30, end=60, strand="-"),
    ]
    rows = interval_features(
        [region], peaks, evidence_id="ev", feature_prefix="CTCF", motif_strands=True
    )
    values = {r["feature_name"]: r["value"] for r in rows}
    # Union [10,60) covers 50/100 bases, not 60/100 from summing two overlaps.
    assert values["CTCF:overlap_fraction"] == 0.5
    assert values["CTCF:peak_count"] == 2
    assert values["CTCF:plus_motifs"] == values["CTCF:minus_motifs"] == 1
    with pytest.raises(PaceError, match="strand"):
        interval_features(
            [region],
            [{k: v for k, v in peaks[0].items() if k != "strand"}],
            evidence_id="ev",
            feature_prefix="CTCF",
            motif_strands=True,
        )
