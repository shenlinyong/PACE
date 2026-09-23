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
    min_distance=5000,
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
    fitted = fit_binned_prior(
        pooled, d_ref=args.reference_distance, d_min=max(resolution, args.minimum_distance)
    )
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
        "fit_range": [args.min_distance, args.max_distance],
        "zero_pixel_policy": "all_callable_bin_opportunities_included",
    }
    for row in rows:
        row["predicted"] = (
            distance_prior(row["geometric_distance"], manifest)
            if row["geometric_distance"]
            else None
        )
    with output_directory(args.out) as dest:
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
        "a": manifest["a"],
        "resolution": resolution,
    }
