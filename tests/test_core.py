"""Independent hand-calculated checks of activity and relative support."""

import math

import numpy as np
import pytest

from pace_livestock.core import activity, bulk_mean, score, tss_contact
from pace_livestock.errors import PaceError
from pace_livestock.io.methylation import summarize_methylation


def fixture_edges():
    return [
        {"element_id": f"E{i + 1}", "gene_id": gene, "A_used": a, "Cbar": c}
        for gene, contacts in [("G1", [3, 2, 2]), ("G2", [1, 2, 6])]
        for i, (a, c) in enumerate(zip([4, 2, 1], contacts, strict=True))
    ]


@pytest.mark.parametrize(
    "eta,expected",
    [
        (0, [2 / 3, 2 / 9, 1 / 9, 2 / 7, 2 / 7, 3 / 7]),
        (1, [18 / 23, 4 / 23, 1 / 23, 2 / 15, 4 / 15, 9 / 15]),
    ],
)
def test_hand_calculation(eta, expected):
    rows, _ = score(fixture_edges(), eta=eta)
    # Independently derived rational values: supplied review §3.4.
    np.testing.assert_allclose([r["pace_score"] for r in rows], expected, atol=1e-10, rtol=0)


def test_measurement_aggregation_precedes_activity():
    first = activity(bulk_mean([[9, 1], [1, 9]]))
    second = activity(bulk_mean([[4, 4], [4, 4]]))
    # Review §3.1: marginal means are (5,5),(4,4), hence 5/9 (not a sum of replicate activities of 3/7).
    assert first / (first + second) == pytest.approx(5 / 9, abs=1e-10)
    assert first / (first + second) != pytest.approx(3 / 7)


def test_activity_and_tss_contact():
    # Review §3.4: sqrt(4*9)=6; .75*2+.25*6=3.
    assert activity([4, 9]) == pytest.approx(6, abs=1e-10)
    assert tss_contact([2, 6], [0.75, 0.25]) == pytest.approx(3, abs=1e-10)


def test_missing_and_true_zero():
    assert math.isnan(activity([0, math.nan]))
    assert activity([0, 9]) == 0
    assert tss_contact([3, math.nan], [1, 0]) == 3
    assert math.isnan(tss_contact([3, math.nan], [0.5, 0.5]))
    with pytest.raises(PaceError):
        activity([-1, 2])


def test_fixed_b_universe_and_eta_skip():
    edges = fixture_edges()
    edges[3]["Cbar"] = math.nan
    a, _ = score(edges, eta=0)
    b, _ = score(edges, eta=1)
    assert math.isfinite(a[0]["pace_score"])
    assert math.isnan(b[0]["pace_score"])
    assert b[0]["reason"] == "fixed_gene_set_contact_missing"


def test_zero_denominator_and_overflow():
    edges = fixture_edges()
    for r in edges:
        r["Cbar"] = 0
    rows, summaries = score(edges, eta=1)
    assert all(r["support"] == 0 and math.isnan(r["pace_score"]) for r in rows)
    assert all(r["normalization_status"] == "zero_support" for r in summaries)
    for r in edges:
        r["A_used"], r["Cbar"] = 1e300, 1e300
    rows, _ = score(edges)
    # Three equal supports normalize to 1/3 even though each raw product overflows.
    np.testing.assert_allclose([r["pace_score"] for r in rows], 1 / 3, atol=1e-10, rtol=0)
    assert all(r["support_status"] == "overflow" for r in rows)


def test_cpg_estimands():
    counts = [
        dict(
            chrom="chr1",
            dyad_start0=p,
            methylated_count=m,
            total_count=n,
            sample_id="S",
            assay="WGBS",
        )
        for p, m, n in [(1, 1, 2), (5, 8, 8)]
    ]
    regions = [dict(element_id="E", chrom="chr1", start=0, end=500)]
    result = summarize_methylation(counts, regions)[0]
    # Review §3.4: mean(1/2,8/8)=.75; (1+8)/(2+8)=.9.
    assert result["M_site"] == 0.75
    assert result["M_pooled"] == 0.9
    assert math.isnan(result["cpg_coverage_fraction"])
    with pytest.raises(PaceError):
        summarize_methylation(counts + counts, regions)
