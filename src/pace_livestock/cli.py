"""Command-line interface for measured activity and contact support."""

from __future__ import annotations

import argparse
import json
import sys

from . import __version__
from .errors import PaceError
from .provenance import clean, output_policy
from .run_options import add_run_options, config_from_args

COMMANDS = {
    "run": "Compute PACE from measured activity and declared contact evidence",
    "validate": "Validate input tables, scope and scientific input contracts",
    "capabilities": "Inspect available contact assets and validation scope",
    "prepare": "Convert genomic annotations and experimental tracks into standard tables",
    "fit-eta": "Learn the optional allocation exponent from independent functional labels",
    "fit-contact-prior": "Fit a distance prior from measured contacts, including true zeros",
    "compare": "Compare measured runs on common candidate denominators",
    "stability": "Measure replicate consistency on common denominators",
    "train": "Train a separate grouped classifier from measured features and functional labels",
    "predict-ml": "Apply a frozen classifier without refitting",
    "benchmark": "Evaluate independent functional labels and explicit baselines",
}


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    measured_alias = bool(argv and argv[0] == "measured")
    if measured_alias:
        argv = ["run", *argv[1:]]
    elif argv and argv[0].startswith("-") and argv[0] not in ("--help", "-h", "--version"):
        argv = ["run", *argv]
    parser = argparse.ArgumentParser(
        prog="pace",
        allow_abbrev=False,
        description="PACE: regulatory support from measured activity and contact evidence.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Start with the offline example:\n"
            "  pace demo --out demo-results\n\n"
            "Analyze measured data:\n"
            "  pace run --config examples/measured/config.yaml --out results\n"
            "  pace measured --help\n\n"
            "Activity requires ATAC, DNase or H3K27ac measurements.\n"
            "Contact may use compatible measurements or an explicitly supplied distance prior.\n"
            "Use --force to preserve and replace an existing PACE output directory."
        ),
    )
    parser.add_argument("--version", action="version", version=f"PACE {__version__}")
    sub = parser.add_subparsers(dest="command", required=True, metavar="COMMAND")
    for name, description in COMMANDS.items():
        p = sub.add_parser(name, help=description, description=description, allow_abbrev=False)
        if name in ("run", "fit-eta", "validate", "capabilities"):
            if measured_alias and name == "run":
                p.prog = "pace measured"
            add_run_options(p, output=name in ("run", "fit-eta"))
        else:
            p.add_argument("--config", required=True, help="YAML; paths resolve relative to it")
            p.add_argument("--out", required=True, help="New output directory")
    demo = sub.add_parser("demo", help="Run the offline measured-data example", allow_abbrev=False)
    demo.add_argument("--regime", choices=["measured"], default="measured")
    demo.add_argument("--out", required=True, help="New output directory")
    custom_commands(sub)
    for child in sub.choices.values():
        if any(a.dest == "out" for a in child._actions):
            child.add_argument(
                "--force",
                action="store_true",
                help="Replace an existing PACE result after preserving it as a sibling backup",
            )
    if not argv:
        parser.print_help()
        return 0
    args = parser.parse_args(argv)
    try:
        with output_policy(force=getattr(args, "force", False)):
            result = dispatch(args)
        print(json.dumps(clean(result), ensure_ascii=False, indent=2, allow_nan=False))
        return 0
    except (PaceError, OSError, KeyError, TypeError, ValueError) as exc:
        print(f"PACE error: {exc}", file=sys.stderr)
        return 2


def dispatch(args):
    command = args.command
    if command in (
        "init",
        "prepare-pairs",
        "merge-tables",
        "normalize-activity",
        "prepare-promoter-weights",
    ):
        from .preparation import (
            init_project,
            merge_tables,
            normalize_activity,
            prepare_pairs,
            promoter_weights,
        )

        return {
            "init": init_project,
            "prepare-pairs": prepare_pairs,
            "merge-tables": merge_tables,
            "normalize-activity": normalize_activity,
            "prepare-promoter-weights": promoter_weights,
        }[command](args)
    if command == "fit-prior":
        if args.config:
            if args.cooler:
                raise PaceError("Choose --config or --cooler, not both")
            from .operations import fit_contact_command

            fit_contact_command(args.config, args.out)
            return {"output": args.out}
        if not all((args.cooler, args.species, args.assembly, args.tissue)):
            raise PaceError(
                "fit-prior needs --cooler, --species, --assembly and --tissue, or --config"
            )
        from .io.prior_fit import fit_cooler_command

        return fit_cooler_command(args)
    if command == "demo":
        from .demo import demo

        return demo(args.regime, args.out)
    if command in ("validate", "capabilities", "run", "fit-eta"):
        from .evidence.assets import capabilities
        from .pipeline import compute, run

        cfg = config_from_args(args)
        if command == "fit-eta" and not cfg["allocation"]["labels_path"]:
            raise PaceError("fit-eta requires --eta-labels or allocation.labels_path")
        if command == "capabilities":
            return capabilities(cfg)
        if command == "validate":
            result = compute(cfg)
            return {"valid": True, "n_candidates": len(result["scores"]), "qc": result["qc"]}
        result = run(cfg, args.out)
        return {
            "output": args.out,
            "n_candidates": len(result["scores"]),
            "n_scoreable": result["qc"]["n_scoreable"],
            "eta": result["eta_calibration"]["eta"],
            "eta_status": result["eta_calibration"]["status"],
        }
    if command in ("prepare", "fit-contact-prior"):
        from .operations import fit_contact_command, prepare_command

        (prepare_command if command == "prepare" else fit_contact_command)(args.config, args.out)
    elif command in ("compare", "stability"):
        from .evaluation.compare import compare_command, stability_command

        (compare_command if command == "compare" else stability_command)(args.config, args.out)
    elif command == "train":
        from .learning.model import train_classifier

        train_classifier(args.config, args.out)
    elif command == "predict-ml":
        from .evaluation.benchmark import predict_ml_command

        predict_ml_command(args.config, args.out)
    elif command == "benchmark":
        from .evaluation.benchmark import benchmark_command

        benchmark_command(args.config, args.out)
    else:
        raise PaceError(f"Command {command} is not implemented")
    return {"output": args.out, "command": command, "status": "complete"}


def custom_commands(sub):
    descriptions = {
        "init": "Create a short project configuration and input templates",
        "prepare-pairs": "Build unique element-promoter query pairs from a prepared catalog",
        "merge-tables": "Combine replicate tables with duplicate-key checks",
        "normalize-activity": "Normalize raw window counts by filtered library size and window length",
        "prepare-promoter-weights": "Freeze TSS weights from measured promoter signals",
        "fit-prior": "Fit contact decay directly from cool/mcool, including unstored zero pixels",
    }
    parsers = {}
    for name, description in descriptions.items():
        p = sub.add_parser(name, help=description, description=description, allow_abbrev=False)
        p.add_argument("--out", required=True)
        parsers[name] = p
    for name in ("init", "fit-prior"):
        p = parsers[name]
        for flag in ("species", "assembly", "tissue"):
            p.add_argument("--" + flag, required=name == "init")
        p.add_argument(
            "--target-level", choices=["individual", "population_mean"], default="individual"
        )
    p = parsers["init"]
    p.add_argument("--catalog-dir")
    p.add_argument(
        "--panel", nargs="+", choices=["ATAC", "DNase", "H3K27ac"], default=["ATAC", "H3K27ac"]
    )
    p.add_argument("--contact-prior")
    p.add_argument("--prior-preset", choices=["abc_human"])
    parsers["prepare-pairs"].add_argument("--catalog-dir", required=True)
    p = parsers["merge-tables"]
    p.add_argument("--inputs", nargs="+", required=True)
    p.add_argument(
        "--table",
        choices=["observed_activity", "observed_contacts", "expression", "methylation"],
        required=True,
    )
    p = parsers["normalize-activity"]
    for flag in ("counts", "library-sizes", "units"):
        p.add_argument("--" + flag, required=True)
    p.add_argument("--normalization-id", default="CPM_density_v1")
    p.add_argument("--window-id")
    p = parsers["prepare-promoter-weights"]
    for flag in ("promoters", "signals", "normalization-id"):
        p.add_argument("--" + flag, required=True)
    p.add_argument("--assay", choices=["ATAC", "DNase", "H3K4me3", "CAGE"], required=True)
    p.add_argument("--zero-policy", choices=["error", "equal"], default="error")
    p = parsers["fit-prior"]
    p.add_argument("--config")
    p.add_argument("--cooler", help="cool file or mcool::/resolutions/N URI")
    p.add_argument("--balanced", action=argparse.BooleanOptionalAction, default=True)
    p.add_argument("--min-distance", type=int, default=5000)
    p.add_argument("--max-distance", type=int, default=5_000_000)
    p.add_argument("--distance-bins", type=int, default=30)
    p.add_argument("--reference-distance", type=int, default=5000)
    p.add_argument("--minimum-distance", type=int, default=5000)
    p.add_argument("--test-chromosomes", nargs="*", default=[])
    p.add_argument("--scale", default="cooler_native")
    p.add_argument("--normalization-id", default="cooler_native")
    p.add_argument("--model-id", default="fitted_contact_prior")
    p.add_argument(
        "--synthetic", action="store_true", help="Label a toy cooler asset for demonstration only"
    )


if __name__ == "__main__":
    raise SystemExit(main())
