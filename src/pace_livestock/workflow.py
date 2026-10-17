"""One-step prediction: peaks + GTF + bigWig (+ Hi-C or a prior) -> PACE scores."""

from __future__ import annotations

import itertools
import sys
import time
from pathlib import Path
from types import SimpleNamespace

from .errors import PaceError
from .provenance import output_directory

CONTACT_HELP = (
    "Choose how contact is measured: --hic FILE (Hi-C .cool/.mcool), --prior DIR (a prior "
    "fitted with 'pace fit-prior' on Hi-C from the same species and assembly), or "
    "--distance-only (contact from the element-TSS distance alone)"
)


def log(message: str) -> None:
    print(f"[pace {time.strftime('%H:%M:%S')}] {message}", file=sys.stderr, flush=True)


def hic_resolution(path: str, requested: int | None) -> int:
    """Resolution to use: the requested one, or the only one stored in a .cool file."""
    if requested:
        return requested
    try:
        import cooler
    except ImportError as exc:
        raise PaceError("Hi-C input requires pip install 'pace-livestock[io]'") from exc
    if "::" not in path and path.endswith(".mcool"):
        stored = sorted(int(Path(u).name) for u in cooler.fileops.list_coolers(path))
        raise PaceError(
            f"{path} stores several resolutions ({', '.join(map(str, stored))}); "
            "choose one with --hic-resolution"
        )
    return int(cooler.Cooler(path).binsize)


def check_predict(args):
    tracks = [
        (assay, path) for option, assay in ASSAY_OPTIONS for path in getattr(args, option) or []
    ]
    if not tracks:
        raise PaceError("Give at least one activity bigWig: --atac, --dnase or --h3k27ac")
    assays = {assay for assay, _ in tracks}
    if "ATAC" in assays and "DNase" in assays:
        raise PaceError("Use ATAC or DNase as the accessibility assay, not both")
    if not (args.hic or args.prior or args.abc_prior):
        raise PaceError("No contact information. " + CONTACT_HELP)
    if bool(args.ctcf) != bool(args.genome):
        raise PaceError("The CTCF boundary constraint needs both --ctcf peaks and --genome FASTA")
    if args.ctcf and not args.hic:
        raise PaceError(
            "--ctcf and --genome add CTCF domain boundaries whose effect is fitted on --hic; "
            "for --prior, give them to 'pace fit-prior'"
        )
    if args.abc_prior and (args.hic or args.prior):
        raise PaceError("--distance-only replaces Hi-C and fitted priors; do not combine them")
    if args.strict_contacts and not args.hic:
        raise PaceError("--strict-contacts applies to --hic input only")
    if args.strict_contacts and args.prior:
        raise PaceError(
            "--strict-contacts uses measured contacts only; do not combine with --prior"
        )
    if args.bin_level_contacts and not args.hic:
        raise PaceError("--bin-level-contacts applies to --hic input only")
    if args.no_shrinkage and not args.hic:
        raise PaceError("--no-shrinkage applies to --hic input only")
    for _, path in tracks:
        if not Path(path).is_file():
            raise PaceError(f"bigWig not found: {path}")
    for path in [*args.peaks, args.gtf, args.chrom_sizes, args.hic and args.hic.split("::")[0]]:
        if path and not Path(path).exists():
            raise PaceError(f"File not found: {path}")
    return tracks


ASSAY_OPTIONS = (("atac", "ATAC"), ("dnase", "DNase"), ("h3k27ac", "H3K27ac"))


def unique_sample_ids(tracks):
    from .cli import stem

    ids, seen = [], set()
    for assay, path in tracks:
        name = stem(path)
        if name in seen:
            name = f"{assay}_{name}"
        base, n = name, 1
        while name in seen:
            n += 1
            name = f"{base}_{n}"
        seen.add(name)
        ids.append(name)
    return ids


def predict(args):
    from .chunked import run_by_chromosome
    from .cli import cooler_uri, stem
    from .config import load_config
    from .io.prior_fit import FIT_MAX_DISTANCE, fit_cooler_command
    from .io.tables import read_table
    from .operations import prepare_command
    from .pipeline import compute, write_results
    from .preparation import merge_tables
    from .run_options import infer_contact_scale

    tracks = check_predict(args)
    prior = args.prior
    chrom_sizes = args.chrom_sizes or tracks[0][1]
    resolution = hic_resolution(args.hic, args.hic_resolution) if args.hic else None
    with output_directory(args.out) as dest:
        prepared = dest / "prepared"
        prepared.mkdir()
        log("1/4 building candidate elements and gene TSSs")
        catalog_dir = prepared / "catalog"
        catalog = prepare_command(
            {
                "kind": "catalog",
                "bed": args.peaks,
                "gtf": args.gtf,
                "chrom_sizes": chrom_sizes,
                "width": args.width,
                "offset": args.offset,
                "radius": args.radius,
                "include_promoters": args.include_promoters,
                "aliases": args.aliases,
                "gene_types": args.gene_types,
                "skip_unlisted_chroms": args.skip_unlisted_chroms,
            },
            catalog_dir,
        )
        log(f"    {catalog['units']} elements, {catalog['candidates']} element-gene pairs")
        activity_tables = []
        for (assay, path), sample in zip(tracks, unique_sample_ids(tracks), strict=True):
            log(f"2/4 quantifying {assay} signal: {path}")
            out = prepared / f"activity_{sample}"
            prepare_command(
                {
                    "kind": "bigwig",
                    "track": path,
                    "units": str(catalog_dir / "units.tsv"),
                    "sample_id": sample,
                    "assay": assay,
                    "unit": "normalized_signal",
                    "normalization_id": "bigwig_signal",
                    "window_id": f"grid:{args.width}:mean",
                    "missing_is_measured_zero": args.missing_as_zero,
                    "minimum_callable_fraction": args.min_callable_fraction,
                },
                out,
            )
            activity_tables.append(str(out / "observed_activity.tsv"))
            warn_sparse_track(out / "observed_activity.tsv", args.missing_as_zero)
        merge_tables(
            SimpleNamespace(
                inputs=activity_tables, table="observed_activity", out=prepared / "activity"
            )
        )
        contacts = None
        if args.hic:
            log(f"3/4 extracting Hi-C contacts at {resolution} bp: {args.hic}")
            balanced = not args.no_balance
            sample = stem(args.hic)
            prepare_command(
                {
                    "kind": "cooler",
                    "contact": cooler_uri(args.hic, resolution),
                    "catalog_dir": str(catalog_dir),
                    "resolution": resolution,
                    "balanced": balanced,
                    "missing_pixels_are_zero": True,
                    "sample_id": sample,
                    "source_id": sample,
                    "scale": "balanced_contact" if balanced else "raw_contact",
                    "normalization_id": "cooler_weight" if balanced else "none",
                },
                prepared / "hic",
            )
            contacts = str(prepared / "hic" / "observed_contacts.tsv")
            if not (args.prior or args.strict_contacts):
                log(
                    "    fitting the distance decay of this Hi-C map"
                    + (
                        " and the CTCF boundary effect"
                        if args.ctcf and not args.no_boundaries
                        else ""
                    )
                )
                fit_cooler_command(
                    SimpleNamespace(
                        cooler=cooler_uri(args.hic, resolution),
                        balanced=balanced,
                        min_distance=None,
                        max_distance=min(args.radius, FIT_MAX_DISTANCE),
                        distance_bins=30,
                        reference_distance=5000,
                        minimum_distance=5000,
                        test_chromosomes=[],
                        scale="balanced_contact" if balanced else "raw_contact",
                        normalization_id="cooler_weight" if balanced else "none",
                        model_id="same_map_distance_prior",
                        no_boundaries=args.no_boundaries,
                        ctcf=args.ctcf,
                        genome=args.genome,
                        synthetic=False,
                        species=args.species,
                        assembly=args.assembly,
                        tissue=args.tissue,
                        target_level="individual",
                        out=prepared / "contact_prior",
                    )
                )
                prior = prepared / "contact_prior"
        else:
            log("3/4 no Hi-C: contact from the distance-decay prior")
        # Hi-C contacts are shrunk pair by pair toward the distance prior (Gamma-Poisson),
        # so bin pairs supported by few reads follow the prior and well-sampled ones the data.
        shrink = bool(contacts) and not (args.no_shrinkage or args.strict_contacts)
        contact = {"mode": "shrinkage" if shrink else "observed" if contacts else "prior_only"}
        if shrink:
            contact.update(reliability="per_pair", reliability_source="gamma_poisson")
        contact["within_bin_distance"] = not args.bin_level_contacts
        if args.strict_contacts:
            # Measured bin-pair contacts only: no pseudocount, and same-bin or masked pairs stay NA.
            contact.update(
                near_diagonal_policy="unresolved",
                pseudocount="none",
                allow_prior_fallback=False,
                within_bin_distance=False,
            )
        if prior:
            contact["prior_path"] = str(Path(prior).resolve())
            # Masked (unbalanced) bins and missing pixels take the fitted expectation,
            # as ABC does; every such pair is labelled contact_prior in the outputs.
            contact["allow_prior_fallback"] = bool(contacts) and not args.strict_contacts
            # Case B: the prior is fitted on Hi-C of another tissue of the same species.
            contact["allow_cross_context_prior"] = bool(args.prior)
        if args.abc_prior:
            contact["prior_preset"] = "abc_human"
        cfg = load_config(
            None,
            overrides={
                "run_id": args.run_id,
                "target_level": "individual",
                "context": {
                    "species": args.species,
                    "assembly": args.assembly,
                    "context_id": args.tissue,
                },
                "inputs": {
                    "units": str(catalog_dir / "units.tsv"),
                    "region_membership": str(catalog_dir / "region_membership.tsv"),
                    "promoters": str(catalog_dir / "promoters.tsv"),
                    "candidates": str(catalog_dir / "candidates.tsv"),
                    "observed_activity": str(prepared / "activity" / "observed_activity.tsv"),
                    "observed_contacts": contacts,
                },
                "catalog": {
                    "profile": "canonical_grid",
                    "width_bp": args.width,
                    "offset_bp": args.offset,
                    "include_promoter_units": args.include_promoters,
                    "chrom_sizes_path": str(catalog_dir / "chrom_sizes.tsv"),
                    "candidate_radius_bp": args.radius,
                },
                "activity": {"panel": sorted({assay for assay, _ in tracks})},
                "contact": contact,
                "scoring": {"partial_policy": args.partial_policy},
                "promoters": {"weights": args.tss_weights},
            },
        )
        if not args.abc_prior:
            infer_contact_scale(cfg)
        if shrink and args.chunk_pairs:
            # Chromosome chunks share one kappa, fitted on every chromosome.
            cfg["contact"]["kappa"] = genome_kappa(cfg)
            log(f"    contact shrinkage: kappa = {cfg['contact']['kappa']:.3g}")
        log("4/4 scoring element-gene pairs")
        if args.chunk_pairs:
            result = run_by_chromosome(
                cfg,
                args.out,
                max_pairs=args.chunk_pairs,
                dest=dest,
                relative_to=dest,
                log=log,
                threads=args.threads,
            )
            qc = result["qc"]
        else:
            result = compute(cfg)
            write_results(dest, cfg, result, relative_to=dest)
            qc = result["qc"]
            del result
        status = {}
        for row in read_table(dest / "gene_summary.tsv", required=["normalization_status"]):
            status[row["normalization_status"]] = status.get(row["normalization_status"], 0) + 1
    log(f"done: {args.out}/scores.tsv.gz")
    return {
        "output": args.out,
        "scores": str(Path(args.out) / "scores.tsv.gz"),
        "n_candidates": qc["n_candidates"],
        "n_scoreable": qc["n_scoreable"],
        "genes_by_status": status,
        "contact_mode": cfg["contact"]["mode"],
        "activity_panel": cfg["activity"]["panel"],
    }


def genome_kappa(cfg):
    """Gamma-Poisson kappa from all unique bin pairs of the run, read as a stream."""
    from .evidence.assets import load_asset
    from .evidence.resolve import run_kappa
    from .io.tables import integer, iter_table, read_table

    prior = load_asset(cfg["contact"]["prior_path"], cfg, kind="contact_prior")
    if prior.get("kappa") is not None:
        return prior["kappa"]
    units = {
        r["element_id"]: {"chrom": r["chrom"], "anchor0": integer(r["anchor0"], "anchor0")}
        for r in read_table(cfg["inputs"]["units"], required=["element_id", "chrom", "anchor0"])
    }
    promoters = {
        r["promoter_id"]: {"tss0": integer(r["tss0"], "tss0")}
        for r in read_table(cfg["inputs"]["promoters"], required=["promoter_id", "tss0"])
    }
    rows = iter_table(cfg["inputs"]["observed_contacts"])
    first = next(rows, None)
    if first is None:
        raise PaceError("No Hi-C contacts to fit the shrinkage kappa")
    return run_kappa(
        itertools.chain([first], rows),
        units,
        promoters,
        prior,
        resolution=integer(first["resolution"], "resolution", minimum=1),
    )


def warn_sparse_track(table, missing_as_zero):
    """Point out bigWigs that leave out zero-signal regions (MACS bedGraph conversions)."""
    from .io.tables import read_table

    rows = read_table(table, required=["measurement_status", "callable_fraction"])
    empty = sum(r["measurement_status"] != "observed" for r in rows)
    partly = sum(float(r["callable_fraction"] or 0) < 1 for r in rows)
    if rows and not missing_as_zero and partly / len(rows) > 0.05:
        log(
            f"    note: {partly} of {len(rows)} elements are not fully covered by the bigWig "
            f"({empty} not at all). Uncovered bases are ignored (partly covered elements are "
            "averaged over covered bases, empty ones stay NA and withhold their genes). If "
            "this bigWig omits zero-signal regions (e.g. converted from a MACS bedGraph), "
            "rerun with --missing-as-zero"
        )
