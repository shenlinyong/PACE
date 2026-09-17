"""Genotype missingness, callability, phase, REF and fixed-target mapping contracts."""

import pytest

from pace_livestock.errors import PaceError
from pace_livestock.io.variants import Reference, VariantIndex, read_variants
from pace_livestock.sequence.genome import Callability, build_window


@pytest.fixture
def reference(tmp_path):
    path = tmp_path / "reference.fa"
    path.write_text(">chr1\n" + "\n".join(["A" * 80] * 200) + "\n")
    with Reference(path) as fasta:
        yield fasta


def unit():
    return {"element_id": "E", "chrom": "chr1", "start": 5000, "end": 5500}


def variant(**overrides):
    return {
        "chrom": "chr1",
        "pos0": 5100,
        "end": 5101,
        "ref": "A",
        "alts": ("C",),
        "gt": (0, 1),
        "phased": True,
        "phase_set": "1",
        "filter": "PASS",
        **overrides,
    }


def build(reference, variants, **kwargs):
    return build_window(
        reference,
        unit(),
        variants,
        input_length=8192,
        ploidy=kwargs.pop("ploidy", 2),
        callable_fraction=kwargs.pop("callable_fraction", 1),
        **kwargs,
    )


def test_snv_changes_only_declared_copy(reference):
    result = build(reference, [variant()])
    assert result["status"] == "resolved"
    # Target starts at offset (8192-500)/2=3846; genomic 5100 is target offset 100.
    assert result["sequences"][0][3946] == "A"
    assert result["sequences"][1][3946] == "C"
    assert all(len(s) == 8192 for s in result["sequences"])


@pytest.mark.parametrize(
    "record,reason",
    [
        (variant(gt=(None, None)), "missing_genotype"),
        (variant(phased=False), "unphased_heterozygote"),
        (variant(gt=(1,)), "ploidy_mismatch"),
        (variant(alts=("<DEL>",), end=6000), "reported_sv_unsupported"),
        (variant(alts=("AC",)), "indel_changes_output_target"),
        (variant(filter="LowQual"), "filtered_variant"),
    ],
)
def test_unsupported_states(reference, record, reason):
    result = build(reference, [record])
    assert result["status"] == "unresolved"
    assert result["reason"] == reason
    assert result["sequences"] == []


def test_ref_and_allele_validation(reference):
    with pytest.raises(PaceError, match="REF mismatch"):
        build(reference, [variant(ref="T")])
    with pytest.raises(PaceError, match="allele index"):
        build(reference, [variant(gt=(0, 3))])


def test_multiallelic_haploid(reference):
    result = build(reference, [variant(alts=("C", "G"), gt=(2,))], ploidy=1)
    assert len(result["sequences"]) == 1
    assert result["sequences"][0][3946] == "G"


def test_phase_blocks_and_overlap(reference):
    result = build(reference, [variant(), variant(pos0=5200, end=5201, phase_set="2")])
    assert result["reason"] == "unlinked_phase_blocks"
    result = build(reference, [variant(), variant()])
    assert result["reason"] == "overlapping_variants"


def test_uncallable_is_not_reference(reference):
    assert build(reference, [], callable_fraction=0.9)["reason"] == "uncallable_input_window"
    assumed = build(reference, [], callable_fraction=0.9, unrecorded_policy="assume_reference")
    assert assumed["status"] == "resolved"
    assert assumed["assumed_reference_fraction"] == pytest.approx(0.1)
    assert (
        build(reference, [variant(gt=(None, None))], unrecorded_policy="assume_reference")["reason"]
        == "missing_genotype"
    )


def test_flanking_insertion_preserves_target(reference):
    result = build(reference, [variant(pos0=4000, end=4001, alts=("AC",), gt=(1, 1))])
    assert result["status"] == "resolved"
    assert result["mapping_status"] == "fixed_target_flanking_indel"
    # The original central A-only target remains 500 bases, after deterministic recentering.
    assert result["sequences"][0][3846:4346] == "A" * 500


def test_flanking_deletion_with_verified_margin(reference):
    deletion = variant(pos0=4000, end=4002, ref="AA", alts=("A",), gt=(1, 1))
    assert build(reference, [deletion])["reason"] == "indel_requires_additional_context"
    result = build(reference, [deletion], context_margin=1)
    assert result["status"] == "resolved"
    assert len(result["sequences"][0]) == 8192
    assert result["sequences"][0][3846:4346] == "A" * 500


def test_spanning_sv_masks_window_and_reference_gt_does_not(reference):
    structural = variant(pos0=0, end=8000, alts=("<DEL>",), gt=(0, 1))
    assert build(reference, [structural])["structural_status"] == "reported_sv_unsupported"
    assert build(reference, [variant(alts=("<DEL>",), gt=(0, 0))])["status"] == "resolved"


def test_callability_interval_union_and_sv_index():
    c = Callability([dict(chrom="chr1", start=0, end=8), dict(chrom="chr1", start=5, end=10)])
    assert c.fraction("chr1", 0, 20) == 0.5
    index = VariantIndex([variant(pos0=0, end=9000, alts=("<DEL>",)), variant(pos0=200, end=201)])
    assert len(index.query("chr1", 5000, 5500)) == 1


def test_vcf_missing_gt_is_not_reference(tmp_path):
    p = tmp_path / "a.vcf"
    p.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\nchr1\t2\t.\tA\tC,G\t.\tPASS\t.\tGT:PS\t./.:.\n"
    )
    r = read_variants(p)[0]
    assert r["pos0"] == 1
    assert r["gt"] == (None, None)
    assert r["alts"] == ("C", "G")


def test_reference_index_across_lines(reference):
    assert reference.fetch("chr1", 79, 83) == "AAAA"
    assert reference.fetch("chr1", 15998, 16000) == "AA"
    with pytest.raises(PaceError):
        reference.fetch("chr1", 0, 16001)
