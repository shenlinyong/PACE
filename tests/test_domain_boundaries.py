"""CTCF domain boundaries and the boundary effect fitted on Hi-C."""

import json
import math

import numpy as np
import pytest

from pace_livestock.boundary_prior import CTCF_PFM, BoundaryIndex, contact_prior, ctcf_boundaries
from pace_livestock.cli import main
from pace_livestock.errors import PaceError
from pace_livestock.io.tables import read_table

RESOLUTION = 10_000
N_BINS = 600
EDGES = (150, 300, 450)  # first bin of each new domain
DROP = 0.4  # contact factor for a pair separated by one boundary


def motif(strand):
    best = "".join(max("ACGT", key=lambda b: CTCF_PFM[b][i]) for i in range(len(CTCF_PFM["A"])))
    return best if strand == "+" else best[::-1].translate(str.maketrans("ACGT", "TGCA"))


def domain_of(i):
    return sum(i >= e for e in EDGES)


def expected(i, j):
    crossed = abs(domain_of(i) - domain_of(j))
    return 1000.0 / max(j - i, 1) * DROP ** min(crossed, 1)


@pytest.fixture(scope="module")
def domains_cool(tmp_path_factory):
    cooler = pytest.importorskip("cooler")
    pd = pytest.importorskip("pandas")
    path = tmp_path_factory.mktemp("hic") / "domains.cool"
    bins = pd.DataFrame(
        {
            "chrom": "chr1",
            "start": np.arange(N_BINS) * RESOLUTION,
            "end": (np.arange(N_BINS) + 1) * RESOLUTION,
        }
    )
    i, j = np.triu_indices(N_BINS, 1)
    keep = j - i <= 300
    i, j = i[keep], j[keep]
    counts = np.array([expected(a, b) for a, b in zip(i, j, strict=True)])
    pixels = pd.DataFrame({"bin1_id": i, "bin2_id": j, "count": counts})
    cooler.create_cooler(str(path), bins, pixels, dtypes={"count": "float64"}, ordered=True)
    return path


@pytest.fixture(scope="module")
def domains_genome(tmp_path_factory):
    """A reverse CTCF site just left and a forward site just right of each domain edge."""
    root = tmp_path_factory.mktemp("genome")
    sequence = bytearray(b"AT" * (N_BINS * RESOLUTION // 2))
    peaks = []
    for e in EDGES:
        edge = e * RESOLUTION
        for strand, start in (("-", edge - 100 - len(motif("-"))), ("+", edge + 100)):
            sequence[start : start + len(motif(strand))] = motif(strand).encode()
            peaks.append(f"chr1\t{start - 50}\t{start + 70}\tp\t0\t.\t5.0\n")
    lines = [sequence[i : i + 80].decode() for i in range(0, len(sequence), 80)]
    (root / "genome.fa").write_text(">chr1\n" + "\n".join(lines) + "\n")
    (root / "ctcf.bed").write_text("".join(peaks))
    return root


def ctcf(root):
    return ["--ctcf", root / "ctcf.bed", "--genome", root / "genome.fa"]


def fit(domains_cool, tmp_path, *extra):
    out = tmp_path / "prior"
    argv = ["fit-prior", "--cooler", domains_cool, "--no-balanced", "--max-distance", "2000000"]
    argv += ["--species", "pig", "--assembly", "test", "--tissue", "liver", "-o", out, *extra]
    assert main([str(a) for a in argv]) == 0
    return out, json.loads((out / "manifest.json").read_text())


def test_ctcf_boundaries_and_beta_fitted_on_hic(domains_cool, domains_genome, tmp_path):
    out, prior = fit(domains_cool, tmp_path, *ctcf(domains_genome))
    assert prior["boundary_source"] == "ctcf_chip_motif"
    positions = [int(r["position0"]) for r in read_table(out / "boundaries.tsv")]
    assert positions == [e * RESOLUTION for e in EDGES]
    # One full-strength boundary multiplies contact by DROP.
    assert prior["beta"] == pytest.approx(-math.log(DROP), rel=0.1)
    assert prior["gamma"] == pytest.approx(1.0, rel=0.1)


def test_without_ctcf_the_prior_is_distance_only(domains_cool, tmp_path):
    out, prior = fit(domains_cool, tmp_path)
    assert "beta" not in prior and not (out / "boundaries.tsv").exists()


@pytest.mark.parametrize("option", ["--ctcf", "--genome"])
def test_constraint_needs_both_peaks_and_genome(domains_genome, option):
    from types import SimpleNamespace

    from pace_livestock.io.prior_fit import map_boundaries

    args = SimpleNamespace(ctcf=None, genome=None)
    setattr(args, option[2:], domains_genome / "x")
    with pytest.raises(PaceError, match="both --ctcf"):
        map_boundaries(args)


def test_no_boundaries_fits_distance_only(domains_cool, domains_genome, tmp_path):
    out, prior = fit(domains_cool, tmp_path, *ctcf(domains_genome), "--no-boundaries")
    assert "beta" not in prior and not (out / "boundaries.tsv").exists()


def test_prior_attenuates_contact_across_a_boundary():
    prior = {"a": 1.0, "gamma": 1.0, "d_ref": 5000, "d_min": 5000, "beta": 1.0}
    prior["boundary_cap"] = 1.0
    index = BoundaryIndex([{"chrom": "c", "position0": 50_000, "strength": 1.0}])
    within = contact_prior("c", 10_000, 40_000, prior, index)
    across = contact_prior("c", 40_000, 70_000, prior, index)
    assert across == pytest.approx(within * math.exp(-1.0))
    two = BoundaryIndex([{"chrom": "c", "position0": p, "strength": 1.0} for p in (50_000, 60_000)])
    # Crossed strength is capped, so a second boundary adds no further loss.
    assert contact_prior("c", 40_000, 70_000, prior, two) == pytest.approx(across)


@pytest.mark.parametrize("left,right,boundaries", [("-", "+", 1), ("+", "-", 0)])
def test_ctcf_boundary_needs_reverse_then_forward_sites(tmp_path, left, right, boundaries):
    filler = "AT" * 500
    sequence = filler + motif(left) + filler + motif(right) + filler
    (tmp_path / "genome.fa").write_text(">c\n" + sequence + "\n")
    first, second = 1000, 2000 + len(motif("+"))
    (tmp_path / "ctcf.bed").write_text(
        f"c\t{first - 50}\t{first + 70}\tp1\t0\t.\t8.0\nc\t{second - 50}\t{second + 70}\tp2\t0\t.\t4.0\n"
    )
    rows = ctcf_boundaries(tmp_path / "ctcf.bed", tmp_path / "genome.fa")
    assert len(rows) == boundaries
    if rows:
        assert first < rows[0]["position0"] < second


def test_fitted_boundary_prior_is_applied_when_scoring(domains_cool, domains_genome, tmp_path):
    from pace_livestock.config import load_config
    from pace_livestock.evidence.assets import load_asset
    from pace_livestock.evidence.resolve import load_boundaries

    out, _ = fit(domains_cool, tmp_path, *ctcf(domains_genome))
    cfg = load_config(
        None, overrides={"context": {"species": "pig", "assembly": "test", "context_id": "liver"}}
    )
    prior = load_asset(out, cfg, kind="contact_prior")
    index = load_boundaries(prior)
    within = contact_prior("chr1", 1_000_000, 1_400_000, prior, index)
    across = contact_prior("chr1", 1_300_000, 1_700_000, prior, index)
    assert across / within == pytest.approx(math.exp(-prior["beta"]), rel=0.05)


def test_reference_accepts_iupac_ambiguity_codes(tmp_path):
    from pace_livestock.io.reference import Reference

    (tmp_path / "g.fa").write_text(">c\nACGTRYKMNacgt\n")
    with Reference(tmp_path / "g.fa") as reference:
        assert reference.fetch("c", 0, 13) == "ACGTRYKMNACGT"


def boundary_cells(levels, strengths, drop, scale=1.0, n=100):
    """Cells at distances 10 kb * 2**k with the expected contact of the boundary prior."""
    cells = []
    for k, level in enumerate(levels):
        d = 10_000 * 2**k
        for m, s in enumerate(strengths[k]):
            cells.append(
                dict(
                    distance_class=k,
                    strength_class=m,
                    n_pairs=n,
                    contact_sum=scale * n * level * (d / 5000) ** -1.0 * drop**s,
                    geometric_distance=d,
                    mean_strength=s,
                )
            )
    return cells


def fit_cells(cells):
    from pace_livestock.io.prior_fit import fit_boundary_prior

    return fit_boundary_prior(cells, d_ref=5000, d_min=5000)


def test_boundary_fit_recovers_beta_and_is_scale_free():
    cells = boundary_cells([1, 1, 1, 1], [[0, 1]] * 4, DROP)
    small = fit_cells(boundary_cells([1, 1, 1, 1], [[0, 1]] * 4, DROP, scale=0.01))
    big = fit_cells(cells)
    assert big["beta_status"] == "fitted"
    assert big["beta"] == pytest.approx(-math.log(DROP), rel=1e-6)
    assert small["beta"] == pytest.approx(big["beta"], rel=1e-6)
    assert small["gamma"] == pytest.approx(big["gamma"], rel=1e-6) == pytest.approx(1.0)
    assert small["a"] == pytest.approx(0.01 * big["a"], rel=1e-6)


@pytest.mark.parametrize("scale", [1.0, 0.01])
def test_boundary_strength_confounded_with_distance_is_not_fitted(scale):
    # One strength per distance class: the strength effect cannot be told from the decay.
    prior = fit_cells(boundary_cells([1, 1, 1, 1], [[0], [0.5], [1], [1]], 1.0, scale=scale))
    assert prior["beta"] == 0 and prior["beta_status"] == "beta_unidentifiable"
    assert prior["gamma"] == pytest.approx(1.0)


def test_cells_without_contacts_count_as_blocking():
    cells = boundary_cells([1, 1, 1], [[0, 1]] * 3, 1.0)
    for r in cells:
        if r["mean_strength"] == 1:
            r["contact_sum"] = 0.0
        else:
            r["contact_sum"] *= 2
    prior = fit_cells(cells)
    # No contact crosses a boundary: the likelihood has no finite maximum.
    assert prior["beta_status"] == "beta_unbounded" and prior["beta"] == 0
    assert prior["gamma"] == pytest.approx(1.0)
    for r in cells:
        if r["mean_strength"] == 1:
            r["contact_sum"] = r["contact_sum"] + 0.5 * r["n_pairs"] / 2 ** r["distance_class"]
    partial = fit_cells(cells)
    assert partial["beta_status"] == "fitted" and 0 < partial["beta"] < 10
    for r in cells:
        if r["mean_strength"] == 1:
            r["contact_sum"] = 1e-9 * r["n_pairs"]
    steep = fit_cells(cells)
    assert steep["beta_status"] == "fitted" and steep["beta"] > 20


def test_boundary_fit_without_strength_variation():
    prior = fit_cells(boundary_cells([1, 1, 1], [[0]] * 3, 1.0))
    assert prior["beta"] == 0 and prior["beta_status"] == "no_boundary_variation"


def test_ctcf_peaks_must_share_chromosome_names_with_the_fasta(domains_genome, tmp_path, capsys):
    renamed = tmp_path / "renamed.bed"
    lines = (domains_genome / "ctcf.bed").read_text().splitlines(keepends=True)
    renamed.write_text("".join(line.replace("chr1", "1", 1) for line in lines))
    with pytest.raises(PaceError, match="chromosome names"):
        ctcf_boundaries(renamed, domains_genome / "genome.fa")
    # Some peaks on other contigs: they are left out and counted.
    renamed.write_text("".join(lines) + "chrUn\t100\t220\tp\t0\t.\t5.0\n")
    assert ctcf_boundaries(renamed, domains_genome / "genome.fa")
    assert "1 of 7 CTCF peaks" in capsys.readouterr().err
