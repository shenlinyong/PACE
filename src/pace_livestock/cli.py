"""Command-line interface for measured activity and contact support."""

from __future__ import annotations

import argparse
import json
import sys

from . import __version__
from .errors import PaceError
from .provenance import clean
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
        prog="PACE",
        allow_abbrev=False,
        description="PACE: regulatory support from measured activity and contact evidence.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Start with the offline example:\n"
            "  PACE demo --out demo-results\n\n"
            "Analyze measured data:\n"
            "  PACE run --config examples/measured/config.yaml --out results\n"
            "  PACE measured --help\n\n"
            "Activity requires ATAC, DNase or H3K27ac measurements.\n"
            "Contact may use compatible measurements or an explicitly supplied distance prior.\n"
            "Existing output directories are never overwritten."
        ),
    )
    parser.add_argument("--version", action="version", version=f"PACE {__version__}")
    sub = parser.add_subparsers(dest="command", required=True, metavar="COMMAND")
    for name, description in COMMANDS.items():
        p = sub.add_parser(name, help=description, description=description, allow_abbrev=False)
        if name in ("run", "fit-eta", "validate", "capabilities"):
            if measured_alias and name == "run":
                p.prog = "PACE measured"
            add_run_options(p, output=name in ("run", "fit-eta"))
        else:
            p.add_argument("--config", required=True, help="YAML; paths resolve relative to it")
            p.add_argument("--out", required=True, help="New output directory")
    demo = sub.add_parser("demo", help="Run the offline measured-data example", allow_abbrev=False)
    demo.add_argument("--regime", choices=["measured"], default="measured")
    demo.add_argument("--out", required=True, help="New output directory")
    if not argv:
        parser.print_help()
        return 0
    args = parser.parse_args(argv)
    try:
        result = dispatch(args)
        print(json.dumps(clean(result), ensure_ascii=False, indent=2, allow_nan=False))
        return 0
    except (PaceError, OSError, KeyError, TypeError, ValueError) as exc:
        print(f"PACE error: {exc}", file=sys.stderr)
        return 2


def dispatch(args):
    command = args.command
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


if __name__ == "__main__":
    raise SystemExit(main())
