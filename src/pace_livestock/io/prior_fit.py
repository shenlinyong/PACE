"""Fit contact decay directly from sparse cooler pixels and valid bin opportunities."""

import math
from pathlib import Path

import numpy as np

from ..errors import PaceError
from ..evidence.contact import distance_prior, fit_binned_prior
from ..provenance import file_hash, output_directory, write_json
from .tables import write_table


def cooler_distance_bins(
    uri,
    *,
    balanced=True,
    min_distance=None,
    max_distance=5_000_000,
    n_bins=30,
    test_chromosomes=(),
    chunksize=1_000_000,
):
    try:
        import cooler
    except ImportError as exc:
        raise PaceError("Fitting a cool/mcool prior requires the io extra") from exc
    c = cooler.Cooler(str(uri))
    if c.binsize is None or c.info.get("storage-mode") != "symmetric-upper":
        raise PaceError("Prior fitting requires fixed-size symmetric-upper cooler storage")
    resolution = int(c.binsize)
    min_distance = resolution if min_distance is None else min_distance
    if min_distance < resolution or max_distance <= min_distance or n_bins < 2:
        raise PaceError(
            "Fit distances require resolution <= min_distance < max_distance and at least two bins"
        )
    unknown = set(test_chromosomes) - set(c.chromnames)
    if unknown or set(test_chromosomes) == set(c.chromnames):
        raise PaceError("Held-out chromosome names must exist and leave training chromosomes")
    bins = c.bins()[:]
    if balanced and "weight" not in bins:
        raise PaceError("Balanced fitting requires cooler weight values")
    weights = bins["weight"].to_numpy(float) if balanced else np.ones(len(bins))
    valid = np.isfinite(weights) & (weights > 0)
    edges = np.geomspace(min_distance, max_distance, n_bins + 1)
    # Force exact user endpoints despite rounding by geomspace.
    edges[0], edges[-1] = min_distance, max_distance
    chrom_codes = np.empty(len(bins), dtype=int)
    records, sums = {}, {}
    for code, chrom in enumerate(c.chromnames):
        first, last = c.extent(chrom)
        chrom_codes[first:last] = code
        mask = valid[first:last].astype(float)
        n = len(mask)
        fft_size = 1 << max(1, (2 * n - 1).bit_length())
        transform = np.fft.rfft(mask, n=fft_size)
        opportunities = np.rint(np.fft.irfft(transform * transform.conj(), n=fft_size)[:n])
        opportunities = np.maximum(opportunities, 0)
        distances = np.arange(n) * resolution
        index = np.searchsorted(edges, distances, side="right") - 1
        use = (index >= 0) & (index < n_bins)
        counts = np.bincount(index[use], weights=opportunities[use], minlength=n_bins)
        logs = np.bincount(
            index[use], weights=opportunities[use] * np.log(distances[use]), minlength=n_bins
        )
        records[chrom] = [
            {
                "chrom": chrom,
                "split": "test" if chrom in test_chromosomes else "train",
                "low": float(edges[i]),
                "high": float(edges[i + 1]),
                "n_pairs": int(counts[i]),
                "log_distance_sum": float(logs[i]),
            }
            for i in range(n_bins)
        ]
        sums[chrom] = np.zeros(n_bins)
    for start in range(0, int(c.info["nnz"]), chunksize):
        pixels = c.pixels()[start : start + chunksize]
        i = pixels["bin1_id"].to_numpy(int)
        j = pixels["bin2_id"].to_numpy(int)
        values = pixels["count"].to_numpy(float)
        keep = (chrom_codes[i] == chrom_codes[j]) & valid[i] & valid[j]
        d = (j - i) * resolution
        bucket = np.searchsorted(edges, d, side="right") - 1
        keep &= (bucket >= 0) & (bucket < n_bins)
        i, j, bucket, values = i[keep], j[keep], bucket[keep], values[keep]
        values = values * weights[i] * weights[j]
        if np.any(~np.isfinite(values)) or np.any(values < 0):
            raise PaceError("Invalid contact counts in prior fit")
        for code in np.unique(chrom_codes[i]):
            select = chrom_codes[i] == code
            sums[c.chromnames[code]] += np.bincount(
                bucket[select], weights=values[select], minlength=n_bins
            )
    all_rows = []
    for chrom, rows in records.items():
        for row, total in zip(rows, sums[chrom], strict=True):
            row["contact_sum"] = float(total)
            row["mean_contact"] = float(total / row["n_pairs"]) if row["n_pairs"] else None
            row["geometric_distance"] = (
                math.exp(row["log_distance_sum"] / row["n_pairs"]) if row["n_pairs"] else None
            )
            row["fit_status"] = (
                "used"
                if row["n_pairs"] and total > 0
                else "zero_mean_excluded"
                if row["n_pairs"]
                else "empty"
            )
            all_rows.append(row)
    return all_rows, resolution


S_STEP = 0.25
S_CLASSES = 41


def _strength_class(s):
    return np.where(s <= 0, 0, np.minimum(S_CLASSES - 1, 1 + (s / S_STEP).astype(int)))


def cooler_boundary_cells(
    uri,
    boundaries,
    *,
    balanced=True,
    min_distance=None,
    max_distance=5_000_000,
    n_bins=30,
    test_chromosomes=(),
    chunksize=1_000_000,
):
    """Valid bin pairs and contact sums by distance class and crossed boundary strength.

    The strength crossed by bins i < j is the summed strength of boundaries strictly
    between their centres, as the prior applies it. Unstored pixels count as zeros.
    """
    import cooler

    c = cooler.Cooler(str(uri))
    resolution = int(c.binsize)
    min_distance = resolution if min_distance is None else min_distance
    bins = c.bins()[:]
    weights = bins["weight"].to_numpy(float) if balanced else np.ones(len(bins))
    valid = np.isfinite(weights) & (weights > 0)
    edges = np.geomspace(min_distance, max_distance, n_bins + 1)
    edges[0], edges[-1] = min_distance, max_distance
    by_chrom = {}
    for row in boundaries:
        by_chrom.setdefault(row["chrom"], []).append(
            (int(row["position0"]), float(row["strength"]))
        )
    below, upto = np.zeros(len(bins)), np.zeros(len(bins))
    chrom_codes = np.empty(len(bins), dtype=int)
    train = np.zeros(len(c.chromnames), dtype=bool)
    shape = (n_bins, S_CLASSES)
    n_pairs, contact, log_d, crossed = (np.zeros(shape) for _ in range(4))
    for code, chrom in enumerate(c.chromnames):
        first, last = c.extent(chrom)
        chrom_codes[first:last] = code
        train[code] = chrom not in test_chromosomes
        points = sorted(by_chrom.get(chrom, []))
        positions = np.array([p for p, _ in points], dtype=float)
        cumulative = np.r_[0.0, np.cumsum([s for _, s in points])]
        centres = bins["start"].to_numpy(float)[first:last] + resolution // 2
        below[first:last] = cumulative[np.searchsorted(positions, centres, side="left")]
        upto[first:last] = cumulative[np.searchsorted(positions, centres, side="right")]
        if not train[code]:
            continue
        ok = valid[first:last]
        n = last - first
        for step in range(max(1, int(np.ceil(min_distance / resolution))), n):
            d = step * resolution
            if d >= max_distance:
                break
            k = np.searchsorted(edges, d, side="right") - 1
            if k < 0 or k >= n_bins:
                continue
            both = ok[:-step] & ok[step:]
            s = below[first + step : last] - upto[first : last - step]
            m = _strength_class(s)
            n_pairs[k] += np.bincount(m[both], minlength=S_CLASSES)
            crossed[k] += np.bincount(m[both], weights=s[both], minlength=S_CLASSES)
            log_d[k] += np.bincount(m[both], minlength=S_CLASSES) * np.log(d)
    for start in range(0, int(c.info["nnz"]), chunksize):
        pixels = c.pixels()[start : start + chunksize]
        i = pixels["bin1_id"].to_numpy(int)
        j = pixels["bin2_id"].to_numpy(int)
        values = pixels["count"].to_numpy(float)
        keep = (chrom_codes[i] == chrom_codes[j]) & valid[i] & valid[j] & train[chrom_codes[i]]
        k = np.searchsorted(edges, (j - i) * resolution, side="right") - 1
        keep &= (k >= 0) & (k < n_bins) & ((j - i) * resolution < max_distance)
        i, j, k, values = (
            i[keep],
            j[keep],
            k[keep],
            values[keep] * weights[i[keep]] * weights[j[keep]],
        )
        m = _strength_class(below[j] - upto[i])
        np.add.at(contact, (k, m), values)
    cells = []
    for k in range(n_bins):
        for m in range(S_CLASSES):
            if n_pairs[k, m]:
                cells.append(
                    dict(
                        distance_class=k,
                        strength_class=m,
                        n_pairs=int(n_pairs[k, m]),
                        contact_sum=float(contact[k, m]),
                        geometric_distance=float(np.exp(log_d[k, m] / n_pairs[k, m])),
                        mean_strength=float(crossed[k, m] / n_pairs[k, m]),
                    )
                )
    return cells, resolution


BOUNDARY_CAP = 1.0
BETA_MAX_DISTANCE = 1_000_000


def fit_boundary_prior(
    cells, *, d_ref, d_min, cap=BOUNDARY_CAP, beta_max_distance=BETA_MAX_DISTANCE
):
    """Fit C(d, S) = a (d/d_ref)^-gamma exp(-beta min(S, cap)).

    beta is the log-contact drop per unit of crossed boundary strength between pairs at
    the same distance (one free level per distance class), using pairs within
    beta_max_distance, where enhancer-promoter links lie and the boundary effect is not
    yet saturated. gamma and a then come from the boundary-corrected distance levels.
    """
    used = [r for r in cells if r["contact_sum"] > 0]
    for r in used:
        r["capped_strength"] = min(r["mean_strength"], cap)
    near = [r for r in used if r["geometric_distance"] <= beta_max_distance]
    classes = sorted({r["distance_class"] for r in near})
    beta, status = 0.0, "no_boundary_variation"
    if near and np.ptp([r["capped_strength"] for r in near]) > 0:
        index = {k: i for i, k in enumerate(classes)}
        design = np.zeros((len(near), len(classes) + 1))
        for t, r in enumerate(near):
            design[t, index[r["distance_class"]]] = 1.0
            design[t, -1] = -r["capped_strength"]
        y = np.log([r["contact_sum"] / r["n_pairs"] for r in near])
        w = np.sqrt([r["n_pairs"] for r in near])
        beta = float(np.linalg.lstsq(design * w[:, None], y * w, rcond=None)[0][-1])
        status = "fitted" if beta > 0 else "no_boundary_effect"
        beta = max(beta, 0.0)
    levels = {}
    for r in used:
        level = levels.setdefault(r["distance_class"], [0.0, 0.0, 0.0, 0])
        level[0] += r["contact_sum"]
        level[1] += r["n_pairs"] * math.exp(-beta * r["capped_strength"])
        level[2] += r["n_pairs"] * math.log(r["geometric_distance"])
        level[3] += r["n_pairs"]
    x = np.array([v[2] / v[3] - math.log(d_ref) for v in levels.values()])
    y = np.log([v[0] / v[1] for v in levels.values()])
    w = np.sqrt([v[3] for v in levels.values()])
    if len(x) < 2:
        raise PaceError("Contact prior fit needs at least two nonzero distance classes")
    coef = np.linalg.lstsq(np.column_stack([np.ones(len(x)), x]) * w[:, None], y * w, rcond=None)[0]
    a, gamma = math.exp(coef[0]), -float(coef[1])
    if gamma <= 0 or not math.isfinite(a):
        raise PaceError("Fitted data do not follow a positive decreasing power-law prior")
    return {
        "a": a,
        "gamma": gamma,
        "beta": beta,
        "beta_status": status,
        "boundary_cap": cap,
        "beta_fit_max_distance": beta_max_distance,
        "d_ref": d_ref,
        "d_min": d_min,
    }


def map_boundaries(args):
    """Domain boundaries from CTCF: ChIP/CUT&Tag peaks with --ctcf, else the built-in motif.

    Both need the reference FASTA (--genome). Without it there are no boundaries and the
    prior is distance only.
    """
    from ..boundary_prior import ctcf_boundaries, motif_genome_boundaries

    if args.ctcf and not args.genome:
        raise PaceError("--ctcf needs --genome (FASTA) to orient the CTCF motifs")
    if args.ctcf:
        return ctcf_boundaries(args.ctcf, args.genome), "ctcf_chip_motif"
    if args.genome:
        return motif_genome_boundaries(args.genome), "ctcf_motif"
    return [], None


def fit_cooler_command(args):
    rows, resolution = cooler_distance_bins(
        args.cooler,
        balanced=args.balanced,
        min_distance=args.min_distance,
        max_distance=args.max_distance,
        n_bins=args.distance_bins,
        test_chromosomes=args.test_chromosomes,
    )
    # Pool by distance before regression; a chromosome is never counted as an independent pixel.
    grouped = {}
    for row in rows:
        if row["split"] != "train":
            continue
        group = grouped.setdefault(
            (row["low"], row["high"]),
            {
                "low": row["low"],
                "high": row["high"],
                "n_pairs": 0,
                "contact_sum": 0.0,
                "log_distance_sum": 0.0,
            },
        )
        for key in ("n_pairs", "contact_sum", "log_distance_sum"):
            group[key] += row[key]
    pooled = []
    for row in grouped.values():
        n = row["n_pairs"]
        row.update(
            mean_contact=row["contact_sum"] / n if n else None,
            geometric_distance=math.exp(row["log_distance_sum"] / n) if n else None,
            fit_status="used"
            if n and row["contact_sum"] > 0
            else "zero_mean_excluded"
            if n
            else "empty",
        )
        pooled.append(row)
    d_min = max(resolution, args.minimum_distance)
    fitted = fit_binned_prior(pooled, d_ref=args.reference_distance, d_min=d_min)
    boundaries, boundary_source, cells = [], None, []
    if not getattr(args, "no_boundaries", True):
        boundaries, boundary_source = map_boundaries(args)
    if boundaries:
        cells, _ = cooler_boundary_cells(
            args.cooler,
            boundaries,
            balanced=args.balanced,
            min_distance=args.min_distance,
            max_distance=args.max_distance,
            n_bins=args.distance_bins,
            test_chromosomes=args.test_chromosomes,
        )
        joint = fit_boundary_prior(cells, d_ref=args.reference_distance, d_min=d_min)
        fitted.update(joint)
    source_file = Path(args.cooler.split("::", 1)[0]).resolve()
    source = {
        "sha256": file_hash(source_file),
        "cooler_group": args.cooler.split("::", 1)[1] if "::" in args.cooler else "/",
    }
    manifest = {
        **fitted,
        "model_id": args.model_id,
        "kind": "contact_prior",
        "is_synthetic": args.synthetic,
        "species": args.species,
        "assembly": args.assembly,
        "context_id": args.tissue,
        "target_level": args.target_level,
        "resolution": resolution,
        "scale": args.scale,
        "normalization_id": args.normalization_id,
        "balancing": "balanced" if args.balanced else "unbalanced",
        "window_id": "bin_pair",
        "prior_source": "fitted_cooler",
        "training_sources": [source],
        "training_chromosomes": sorted({r["chrom"] for r in rows if r["split"] == "train"}),
        "test_chromosomes": sorted(args.test_chromosomes),
        "calibration_sources": [],
        "test_sources": [source] if args.test_chromosomes else [],
        "validation": {},
        "fit_range": [
            resolution if args.min_distance is None else args.min_distance,
            args.max_distance,
        ],
        "zero_pixel_policy": "all_callable_bin_opportunities_included",
        # Distances were measured between bins, so the prior is applied at bin centres.
        "fitting_coordinate_policy": "bin_centers",
    }
    if boundaries:
        manifest.update(boundary_table="boundaries.tsv", boundary_source=boundary_source)
    for row in rows:
        row["predicted"] = (
            distance_prior(row["geometric_distance"], manifest)
            if row["geometric_distance"]
            else None
        )
    with output_directory(args.out) as dest:
        if boundaries:
            fields = ["chrom", "position0", "strength", "evidence_type"]
            write_table(dest / "boundaries.tsv", [{f: r[f] for f in fields} for r in boundaries])
            write_table(dest / "boundary_cells.tsv", cells)
            manifest["boundary_sha256"] = file_hash(dest / "boundaries.tsv")
        write_json(dest / "manifest.json", manifest)
        write_table(dest / "distance_bins.tsv", rows)
        write_json(
            dest / "fit_report.json",
            {
                "validation": "chromosome_holdout_contact_decay_diagnostic"
                if args.test_chromosomes
                else "not_held_out",
                "functional_validation": "not_performed",
                "callable_pairs": sum(r["n_pairs"] for r in rows),
                "zero_mean_bins_excluded_from_log_fit": sum(
                    r["fit_status"] == "zero_mean_excluded" for r in rows
                ),
            },
        )
    return {
        "output": args.out,
        "gamma": manifest["gamma"],
        "beta": manifest.get("beta", 0.0),
        "boundaries": len(boundaries),
        "a": manifest["a"],
        "resolution": resolution,
    }
