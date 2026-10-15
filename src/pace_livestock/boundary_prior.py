"""Candidate CTCF boundaries, contact priors and Gamma-Poisson shrinkage."""

from __future__ import annotations

import math
from bisect import bisect_left, bisect_right
from collections import defaultdict
from pathlib import Path

import numpy as np

from .errors import PaceError
from .evidence.contact import distance_prior
from .io.reference import Reference
from .io.tables import integer, number, read_table
from .provenance import file_hash


class BoundaryIndex:
    """Sum point-boundary strengths strictly between two genomic anchors."""

    def __init__(self, rows):
        by_chrom = defaultdict(dict)
        for row in rows:
            chrom = row.get("chrom")
            if not chrom:
                raise PaceError("Boundary chromosome is required")
            position = integer(row.get("position0"), "boundary.position0")
            strength = number(row.get("strength"), "boundary.strength", minimum=0)
            if position in by_chrom[chrom]:
                raise PaceError("Duplicate boundary position")
            by_chrom[chrom][position] = strength
        self.data = {}
        for chrom, values in by_chrom.items():
            positions = sorted(values)
            self.data[chrom] = (positions, np.r_[0.0, np.cumsum([values[p] for p in positions])])

    def strength(self, chrom, left, right):
        low, high = sorted((left, right))
        if chrom not in self.data or low == high:
            return 0.0
        positions, cumulative = self.data[chrom]
        return float(
            cumulative[bisect_left(positions, high)] - cumulative[bisect_right(positions, low)]
        )


def load_boundaries(prior):
    if not prior.get("boundary_table"):
        if number(prior.get("beta", 0), "beta", minimum=0) > 0:
            raise PaceError("Positive boundary beta requires a boundary table")
        return BoundaryIndex([])
    root = Path(prior.get("asset_directory", ".")).resolve()
    path = (root / prior["boundary_table"]).resolve()
    if not path.is_relative_to(root):
        raise PaceError("Boundary table must be inside the prior asset directory")
    if file_hash(path) != prior.get("boundary_sha256"):
        raise PaceError("Boundary table checksum differs from the fitted prior")
    return BoundaryIndex(read_table(path, required=["chrom", "position0", "strength"]))


def contact_prior(chrom, left, right, prior, boundaries=None):
    policy = prior.get(
        "fitting_coordinate_policy", prior.get("prior_coordinate_policy", "genomic_anchors")
    )
    if policy not in ("genomic_anchors", "bin_centers"):
        raise PaceError("Unknown contact prior coordinate policy")
    if policy == "bin_centers":
        resolution = integer(prior.get("resolution"), "prior resolution", minimum=1)
        left, right = (x // resolution * resolution + resolution // 2 for x in (left, right))
    beta = number(prior.get("beta", 0), "contact prior beta", minimum=0)
    base = distance_prior(abs(left - right), prior)
    if beta == 0:
        return base
    if boundaries is None:
        boundaries = load_boundaries(prior)
    crossed = boundaries.strength(chrom, left, right)
    if prior.get("boundary_cap") is not None:
        crossed = min(crossed, number(prior["boundary_cap"], "boundary_cap", minimum=0))
    return base * math.exp(-beta * crossed)


# JASPAR MA0139.1 (CTCF, vertebrates; CC BY 4.0). On the + strand this is the forward
# motif: CTCF loop anchors pair a forward site on the left with a reverse site on the right.
CTCF_PFM = {
    "A": [87, 167, 281, 56, 8, 744, 40, 107, 851, 5, 333, 54, 12, 56, 104, 372, 82, 117, 402],
    "C": [291, 145, 49, 800, 903, 13, 528, 433, 11, 0, 3, 12, 0, 8, 733, 13, 482, 322, 181],
    "G": [76, 414, 449, 21, 0, 65, 334, 48, 32, 903, 566, 504, 890, 775, 5, 507, 307, 73, 266],
    "T": [459, 187, 134, 36, 2, 91, 11, 324, 18, 3, 9, 341, 8, 71, 67, 17, 37, 396, 59],
}


def ctcf_pwm_rows():
    return [{b: CTCF_PFM[b][i] for b in "ACGT"} for i in range(len(CTCF_PFM["A"]))]


def _log_pwm(pwm_rows):
    matrix = np.array([[number(r.get(b), f"PWM.{b}", minimum=0) for b in "ACGT"] for r in pwm_rows])
    if matrix.ndim != 2 or len(matrix) < 2 or np.any(matrix.sum(axis=1) <= 0):
        raise PaceError("PWM needs at least two nonempty A/C/G/T positions")
    matrix = (matrix + 0.01) / (matrix.sum(axis=1, keepdims=True) + 0.04)
    return np.log2(matrix / 0.25)


_LUT = np.full(256, 4, dtype=np.int8)
for _i, _base in enumerate(b"ACGT"):
    _LUT[_base] = _i


def _scores(reference, chrom, start, stop, log_pwm):
    """Both-strand log2-odds of every motif start in [start, stop); invalid starts are -inf."""
    length = len(log_pwm)
    seq = _LUT[np.frombuffer(reference.fetch(chrom, start, stop + length - 1).encode(), np.uint8)]
    valid = np.ones(stop - start, dtype=bool)
    for offset in range(length):
        valid &= seq[offset : offset + stop - start] < 4
    result = {}
    for strand, weights in (("+", log_pwm), ("-", log_pwm[::-1, ::-1])):
        score = np.zeros(stop - start)
        padded = np.column_stack([weights, np.zeros(length)])
        for offset in range(length):
            score += padded[offset, seq[offset : offset + stop - start]]
        result[strand] = np.where(valid, score, -np.inf)
    return result


def scan_motifs(fasta, pwm_rows, *, threshold, chunk_size=1_000_000):
    """Scan both strands using a supplied CTCF position-weight/count matrix.

    Scores are log2 odds against uniform background. Ambiguous bases cannot match.
    Chunk overlap includes every motif once, including across FASTA line boundaries.
    """
    threshold = number(threshold, "motif threshold")
    chunk_size = integer(chunk_size, "chunk_size", minimum=1)
    log_pwm = _log_pwm(pwm_rows)
    maximum = float(log_pwm.max(axis=1).sum())
    if maximum <= 0 or threshold < 0 or threshold > maximum:
        raise PaceError("Motif threshold must be between zero and the maximum log2 odds")
    length = len(log_pwm)
    rows = []
    with Reference(fasta) as reference:
        for chrom, size in reference.sizes.items():
            for start in range(0, max(0, size - length + 1), chunk_size):
                stop = min(start + chunk_size, size - length + 1)
                for strand, score in _scores(reference, chrom, start, stop, log_pwm).items():
                    for i in np.flatnonzero(score >= threshold):
                        rows.append(
                            dict(
                                chrom=chrom,
                                start=start + int(i),
                                end=start + int(i) + length,
                                strand=strand,
                                score=float(score[i]),
                                strength=min(1.0, float(score[i] / maximum)),
                            )
                        )
    return rows


def ctcf_sites(peaks, fasta, *, min_score_fraction=0.6):
    """Strongest CTCF motif in each ChIP peak, oriented by strand.

    Site strength is the peak signal scaled by its 95th percentile (capped at 1), or the
    relative motif score when the peaks carry no signal. Peaks without a motif are dropped.
    """
    log_pwm = _log_pwm(ctcf_pwm_rows())
    maximum = float(log_pwm.max(axis=1).sum())
    length = len(log_pwm)
    signals = [p["signal"] for p in peaks if p.get("signal") is not None]
    scale = float(np.quantile(signals, 0.95)) if signals else 0.0
    sites = []
    with Reference(fasta) as reference:
        for peak in peaks:
            chrom, start = peak["chrom"], int(peak["start"])
            if chrom not in reference.sizes:
                continue
            stop = min(int(peak["end"]), reference.sizes[chrom]) - length + 1
            if stop <= start:
                continue
            scores = _scores(reference, chrom, start, stop, log_pwm)
            strand = max(scores, key=lambda k: scores[k].max())
            i = int(np.argmax(scores[strand]))
            best = float(scores[strand][i])
            if best < min_score_fraction * maximum:
                continue
            strength = (
                min(1.0, peak["signal"] / scale)
                if scale > 0 and peak.get("signal") is not None
                else min(1.0, best / maximum)
            )
            sites.append(
                dict(
                    chrom=chrom,
                    start=start + i,
                    end=start + i + length,
                    strand=strand,
                    strength=max(strength, 0.0),
                )
            )
    return sites


def read_ctcf_peaks(path):
    """BED or narrowPeak; signal from column 7 (narrowPeak) or column 5 when numeric."""
    from .io.bed_gtf import text_open

    peaks = []
    with text_open(path) as handle:
        for n, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise PaceError(f"{path}:{n}: CTCF peaks need chrom, start and end")
            signal = None
            for column in (6, 4):
                if len(fields) > column:
                    try:
                        signal = float(fields[column])
                        break
                    except ValueError:
                        pass
            start = integer(fields[1], "CTCF peak start")
            end = integer(fields[2], "CTCF peak end")
            if end <= start:
                raise PaceError(f"{path}:{n}: CTCF peak end must exceed start")
            peaks.append(dict(chrom=fields[0], start=start, end=end, signal=signal))
    return peaks


def ctcf_boundaries(peaks_path, fasta, *, max_gap_bp=1_000_000):
    """Domain boundaries between adjacent reverse (-) and forward (+) CTCF sites."""
    sites = ctcf_sites(read_ctcf_peaks(peaks_path), fasta)
    rows = motif_boundaries(sites, max_gap_bp=max_gap_bp)
    for row in rows:
        row["evidence_type"] = "ctcf_chip_motif"
    return unique_positions(rows)


def pwm_threshold(log_pwm, pvalue, *, resolution=100):
    """Smallest log2-odds score reached by a random uniform sequence with probability <= pvalue."""
    scores = np.rint(np.asarray(log_pwm) * resolution).astype(int)
    dist, offset = np.ones(1), 0
    for row in scores:
        lo, hi = int(row.min()), int(row.max())
        new = np.zeros(len(dist) + hi - lo)
        for value in row:
            new[value - lo : value - lo + len(dist)] += 0.25 * dist
        dist, offset = new, offset + lo
    tail = np.cumsum(dist[::-1])[::-1]
    index = int(np.argmax(tail <= pvalue)) if np.any(tail <= pvalue) else len(dist) - 1
    return (index + offset) / resolution


def motif_genome_boundaries(fasta, *, pvalue=1e-5, max_gap_bp=1_000_000):
    """Boundaries from a genome-wide scan of the built-in CTCF motif, without ChIP data.

    Sites with p <= pvalue on either strand are paired as in ctcf_boundaries; a site's
    strength is its relative motif score.
    """
    rows_pwm = ctcf_pwm_rows()
    threshold = pwm_threshold(_log_pwm(rows_pwm), pvalue)
    sites = scan_motifs(fasta, rows_pwm, threshold=threshold)
    rows = motif_boundaries(sites, max_gap_bp=max_gap_bp)
    for row in rows:
        row["evidence_type"] = "ctcf_motif"
    return unique_positions(rows)


def unique_positions(rows):
    best = {}
    for row in rows:
        key = row["chrom"], row["position0"]
        if key not in best or row["strength"] > best[key]["strength"]:
            best[key] = row
    return [best[k] for k in sorted(best)]


def motif_boundaries(motifs, *, chip=None, max_gap_bp=1_000_000):
    """Divergent (- then +) adjacent nonoverlapping sites suggest candidate boundaries."""
    max_gap_bp = integer(max_gap_bp, "max_gap_bp", minimum=1)
    peaks = defaultdict(list)
    for peak in chip or []:
        start, end = integer(peak["start"], "peak start"), integer(peak["end"], "peak end")
        if end <= start:
            raise PaceError("Invalid CTCF peak interval")
        peaks[peak["chrom"]].append((start, end))
    # Prefix maxima permit overlap queries even for nested or overlapping ChIP peaks.
    peak_index = {}
    for chrom, intervals in peaks.items():
        intervals.sort()
        peak_index[chrom] = (
            [p[0] for p in intervals],
            np.maximum.accumulate([p[1] for p in intervals]),
        )
    groups = defaultdict(list)
    for row in motifs:
        start, end = integer(row["start"], "motif start"), integer(row["end"], "motif end")
        if end <= start or row["strand"] not in ("+", "-"):
            raise PaceError("Motifs need a positive interval and explicit + or - strand")
        strength = number(row.get("strength", 1), "motif strength", minimum=0, maximum=1)
        if chip is not None:
            starts, ends = peak_index.get(row["chrom"], ([], []))
            i = bisect_left(starts, end) - 1
            if i < 0 or ends[i] <= start:
                continue
        groups[row["chrom"]].append({**row, "start": start, "end": end, "strength": strength})
    output = []
    for chrom, rows in sorted(groups.items()):
        # Collapse overlapping motif calls to their strongest representative.
        rows.sort(key=lambda r: (r["start"], r["end"], r["strand"]))
        clusters, cluster, end = [], [], -1
        for row in rows:
            if cluster and row["start"] >= end:
                clusters.append(cluster)
                cluster = []
            cluster.append(row)
            end = max(end, row["end"])
        if cluster:
            clusters.append(cluster)
        sites = [min(c, key=lambda r: (-r["strength"], r["start"], r["strand"])) for c in clusters]
        for left, right in zip(sites[:-1], sites[1:], strict=True):
            if (
                left["strand"] == "-"
                and right["strand"] == "+"
                and right["start"] - left["end"] <= max_gap_bp
            ):
                output.append(
                    dict(
                        chrom=chrom,
                        position0=(left["end"] + right["start"]) // 2,
                        strength=min(left["strength"], right["strength"]),
                        evidence_type="chip_supported_motif"
                        if chip is not None
                        else "sequence_motif",
                        boundary_status="candidate_unvalidated",
                        left_site_start=left["start"],
                        right_site_start=right["start"],
                    )
                )
    return output


def posterior_contact(raw_count, count_to_contact, prior, kappa):
    """C ~ Gamma(kappa, rate=kappa/prior); Y|C ~ Poisson(C/factor)."""
    y = number(raw_count, "raw_count", minimum=0)
    factor = number(count_to_contact, "count_to_contact", minimum=0)
    mu = number(prior, "prior", minimum=0)
    shape = number(kappa, "kappa", minimum=0)
    if y != math.floor(y) or min(factor, mu, shape) <= 0:
        raise PaceError("Shrinkage requires integer raw counts and positive factor/prior/kappa")
    expected = mu / factor
    if not math.isfinite(expected) or expected <= 0:
        raise PaceError("Expected raw count is outside the numeric range")
    r = 1 / (1 + shape / expected)
    value = r * (y * factor) + (1 - r) * mu
    variance = (y + shape) * (factor * r) ** 2
    if not math.isfinite(value) or not math.isfinite(variance):
        raise PaceError("Contact posterior exceeds the numeric range")
    return dict(
        resolved_value=value,
        reliability=r,
        expected_raw_count=expected,
        posterior_variance=variance,
        kappa=shape,
    )


def fit_kappa(counts, expected):
    y, mu = np.asarray(counts, dtype=float), np.asarray(expected, dtype=float)
    if (
        len(y) < 3
        or y.shape != mu.shape
        or not np.all(np.isfinite(y))
        or not np.all(np.isfinite(mu))
        or np.any(y < 0)
        or np.any(mu <= 0)
    ):
        raise PaceError("Kappa fitting needs at least three finite callable count pairs")
    excess = float(np.sum((y - mu) ** 2 - y) / np.sum(mu**2))
    kappa = min(1e6, max(1e-6, 1 / excess)) if excess > 0 else 1e6
    return {
        "kappa": kappa,
        "moment_excess": excess,
        "n_unique_pairs": len(y),
        "method": "gamma_poisson_moments",
        "at_bound": kappa in (1e-6, 1e6),
    }


def unique_count_pairs(rows, *, resolution, boundaries):
    """Validate and deduplicate physical sample/bin pairs before estimating dispersion."""
    output, seen = [], {}
    for row in rows:
        if row.get("measurement_status") != "observed":
            continue
        if integer(row.get("resolution"), "resolution", minimum=1) != resolution:
            raise PaceError("Contact fitting cannot mix resolutions")
        left, right = integer(row.get("anchor0"), "anchor0"), integer(row.get("tss0"), "tss0")
        left, right = sorted(
            (
                left // resolution * resolution + resolution // 2,
                right // resolution * resolution + resolution // 2,
            )
        )
        if left == right:
            continue
        y = number(row.get("raw_count"), "raw_count", minimum=0)
        factor = number(row.get("count_to_contact"), "count_to_contact", minimum=0)
        if y != math.floor(y) or factor <= 0:
            raise PaceError("Fit requires integer raw counts and positive count_to_contact")
        if not math.isclose(
            number(row.get("contact_value"), "contact_value", minimum=0),
            y * factor,
            rel_tol=1e-8,
            abs_tol=1e-12,
        ):
            raise PaceError("Raw counts and count_to_contact do not reproduce contact_value")
        if not row.get("sample_id") or not row.get("bin_pair_id") or not row.get("chrom"):
            raise PaceError("Contact fitting needs sample_id, bin_pair_id and chrom")
        # Identity is the physical measurement (sample, chromosome, resolution and the
        # two bins), never its label: a renamed copy of one pixel is still one pixel.
        key = row["sample_id"], row["chrom"], resolution, left, right
        signature = (y, factor)
        if key in seen:
            if seen[key] != signature:
                raise PaceError(
                    "Measurements of one sample/bin pair disagree (conflicting raw counts or conversion factors)"
                )
            continue
        seen[key] = signature
        output.append(
            {
                **row,
                "raw_count": y,
                "count_to_contact": factor,
                "distance_bp": right - left,
                "boundary_strength": boundaries.strength(row["chrom"], left, right),
            }
        )
    return output


def fit_contact_grid(rows, *, gamma_grid, beta_grid, d_ref=5000, d_min=5000):
    """Fit amplitude analytically and decay/boundary coefficients on a fixed grid."""
    gammas = parameter_grid(gamma_grid, "gamma", positive=True)
    betas = parameter_grid(beta_grid, "beta")
    if not rows or min(d_ref, d_min) <= 0:
        raise PaceError("Contact fit needs callable rows and positive reference distances")
    y = np.array([r["raw_count"] for r in rows])
    exposure = 1 / np.array([r["count_to_contact"] for r in rows])
    logd = np.log(np.maximum([r["distance_bp"] for r in rows], d_min) / d_ref)
    b = np.array([r["boundary_strength"] for r in rows])
    fits = []
    for gamma in gammas:
        for beta in betas:
            shape = np.exp(-gamma * logd - beta * b)
            denominator = float(np.sum(exposure * shape))
            if denominator <= 0 or y.sum() <= 0 or np.any(shape <= 0):
                raise PaceError("Contact fit has zero support or numerical underflow")
            a = float(y.sum() / denominator)
            expected = exposure * a * shape
            loss = float(np.mean(expected - y * np.log(expected)))
            fits.append(dict(gamma=gamma, beta=beta, a=a, poisson_loss=loss))
    best = min(fits, key=lambda x: (x["poisson_loss"], x["beta"], x["gamma"]))
    expected = exposure * best["a"] * np.exp(-best["gamma"] * logd - best["beta"] * b)
    correlation = float(np.corrcoef(logd, b)[0, 1]) if np.std(logd) > 0 and np.std(b) > 0 else None
    return {
        **best,
        "d_ref": d_ref,
        "d_min": d_min,
        **fit_kappa(y, expected),
        "distance_boundary_correlation": correlation,
        "grid_results": fits,
    }


def parameter_grid(values, name, *, positive=False, maximum=None):
    if not isinstance(values, list) or not values or len(values) > 50:
        raise PaceError(f"{name}_grid must contain between 1 and 50 numbers")
    result = sorted({number(x, name, minimum=0, maximum=maximum) for x in values})
    if positive and result[0] <= 0:
        raise PaceError(f"{name} must be positive")
    return result
