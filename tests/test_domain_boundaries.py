"""Domain boundaries from Hi-C insulation or oriented CTCF sites, and their fitted effect."""

import json
import math

import numpy as np
import pytest

from pace_livestock.boundary_prior import CTCF_PFM, BoundaryIndex, contact_prior, ctcf_boundaries
from pace_livestock.cli import main
from pace_livestock.io.insulation import boundary_calls, insulation_scores
from pace_livestock.io.tables import read_table

RESOLUTION = 10_000
N_BINS = 600
EDGES = (150, 300, 450)  # first bin of each new domain
DROP = 0.4  # contact factor for a pair separated by one boundary


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


def test_insulation_minimum_marks_the_domain_edge():
    n, window = 200, 10
    valid = np.ones(n, dtype=bool)
    band = np.full((2 * window + 1, n), np.nan)
    for k in range(1, 2 * window + 1):
        for a in range(n - k):
            band[k, a] = (1.0 / k) * (0.3 if (a < 100) != (a + k < 100) else 1.0)
    score = insulation_scores(band, valid, window)
    calls = boundary_calls(score, window)
    assert [i for i, _ in calls] == [99]


def fit(domains_cool, tmp_path, *extra):
    out = tmp_path / "prior"
    argv = ["fit-prior", "--cooler", domains_cool, "--no-balanced", "--max-distance", "2000000"]
    argv += ["--species", "pig", "--assembly", "test", "--tissue", "liver", "-o", out, *extra]
    assert main([str(a) for a in argv]) == 0
    return out, json.loads((out / "manifest.json").read_text())


def test_hic_boundaries_and_beta_are_fitted_by_default(domains_cool, tmp_path):
    out, prior = fit(domains_cool, tmp_path, "--insulation-window", "200000")
    assert prior["boundary_source"] == "hic_insulation"
    positions = [int(r["position0"]) for r in read_table(out / "boundaries.tsv")]
    assert len(positions) == len(EDGES)
    assert all(min(abs(p - e * RESOLUTION) for p in positions) <= RESOLUTION for e in EDGES)
    # One full-strength boundary multiplies contact by DROP.
    assert prior["beta"] == pytest.approx(-math.log(DROP), rel=0.15)
    assert prior["gamma"] == pytest.approx(1.0, rel=0.1)


def test_no_boundaries_fits_distance_only(domains_cool, tmp_path):
    out, prior = fit(domains_cool, tmp_path, "--no-boundaries")
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


def motif(strand):
    best = "".join(max("ACGT", key=lambda b: CTCF_PFM[b][i]) for i in range(len(CTCF_PFM["A"])))
    return best if strand == "+" else best[::-1].translate(str.maketrans("ACGT", "TGCA"))


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


def test_fitted_boundary_prior_is_applied_when_scoring(domains_cool, tmp_path):
    from pace_livestock.config import load_config
    from pace_livestock.evidence.assets import load_asset
    from pace_livestock.evidence.resolve import load_boundaries

    out, _ = fit(domains_cool, tmp_path, "--insulation-window", "200000")
    cfg = load_config(
        None, overrides={"context": {"species": "pig", "assembly": "test", "context_id": "liver"}}
    )
    prior = load_asset(out, cfg, kind="contact_prior")
    index = load_boundaries(prior)
    within = contact_prior("chr1", 1_000_000, 1_400_000, prior, index)
    across = contact_prior("chr1", 1_300_000, 1_700_000, prior, index)
    assert across / within == pytest.approx(math.exp(-prior["beta"]), rel=0.05)


def test_motif_threshold_matches_the_random_sequence_tail():
    from pace_livestock.boundary_prior import _log_pwm, ctcf_pwm_rows, pwm_threshold

    log_pwm = _log_pwm(ctcf_pwm_rows())
    rng = np.random.default_rng(0)
    seqs = rng.integers(0, 4, size=(200_000, len(log_pwm)))
    scores = log_pwm[np.arange(len(log_pwm)), seqs].sum(axis=1)
    threshold = pwm_threshold(log_pwm, 1e-3)
    assert np.mean(scores >= threshold) == pytest.approx(1e-3, rel=0.3)


@pytest.mark.parametrize("left,right,boundaries", [("-", "+", 1), ("+", "-", 0)])
def test_builtin_motif_scan_without_ctcf_peaks(tmp_path, left, right, boundaries):
    from pace_livestock.boundary_prior import motif_genome_boundaries

    filler = "AT" * 500
    sequence = filler + motif(left) + filler + motif(right) + filler
    (tmp_path / "genome.fa").write_text(">c\n" + sequence + "\n")
    rows = motif_genome_boundaries(tmp_path / "genome.fa")
    assert len(rows) == boundaries
    assert all(r["evidence_type"] == "ctcf_motif" for r in rows)
