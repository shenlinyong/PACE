"""Thin command-line layer: scientific errors return nonzero, outputs are transactional."""

from __future__ import annotations

import argparse
import json
import sys

from . import __version__
from .errors import PaceError
from .provenance import clean

COMMANDS = {
    "validate": "Validate schemas, provenance, assets and scientific input contracts",
    "capabilities": "Inspect implementation, weights and validation scope",
    "prepare": "Convert genomic files into canonical standard tables",
    "run": "Compute the canonical PACE score in the configured evidence regime",
    "fit-contact-prior": "Fit a zero-inclusive binned distance prior",
    "train-sequence": "Train the quantitative central-target reference CNN",
    "predict-sequence": "Predict quantitative activity from prepared individual windows",
    "fit-fusion": "Fit a target-specific log-space fusion calibrator",
    "prepare-genome": "Construct callable fixed-target individual sequence windows",
    "compare": "Compare runs after recomputing a common denominator",
    "variant-effects": "Run single-variant REF/ALT sequence scenarios",
    "stability": "Measure replicate consistency on common denominators",
    "train": "Train grouped elastic-net classification and a base-only control",
    "predict-ml": "Apply a frozen classifier without refitting",
    "benchmark": "Evaluate independent functional labels and explicit baselines",
}


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="pace-livestock",
        description="PACE — auditable regulatory support for livestock research",
    )
    parser.add_argument("--version", action="version", version=f"PACE {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)
    for name, description in COMMANDS.items():
        p = sub.add_parser(name, help=description, description=description)
        p.add_argument(
            "--config", required=True, help="YAML file; input paths resolve relative to this file"
        )
        if name not in ("validate", "capabilities"):
            p.add_argument("--out", required=True, help="New output directory (never overwritten)")
    d = sub.add_parser("demo", help="Run a fully offline synthetic example")
    d.add_argument("--regime", choices=["measured", "hybrid", "genome_only"], required=True)
    d.add_argument("--out", required=True)
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
    if command in ("validate", "capabilities", "run"):
        from .config import load_config
        from .evidence.assets import capabilities
        from .pipeline import compute, run

        if command == "capabilities":
            return capabilities(load_config(args.config))
        if command == "validate":
            result = compute(load_config(args.config))
            return {"valid": True, "n_candidates": len(result["scores"]), "qc": result["qc"]}
        result = run(args.config, args.out)
        return {
            "output": args.out,
            "n_candidates": len(result["scores"]),
            "n_scoreable": result["qc"]["n_scoreable"],
        }
    if command in (
        "prepare",
        "fit-contact-prior",
        "fit-fusion",
        "prepare-genome",
        "predict-sequence",
    ):
        from . import operations

        functions = {
            "prepare": operations.prepare_command,
            "fit-contact-prior": operations.fit_contact_command,
            "fit-fusion": operations.fit_fusion_command,
            "prepare-genome": operations.prepare_genome_command,
        }
        if command == "predict-sequence":
            operations.prepare_genome_command(args.config, args.out, predict=True)
        else:
            functions[command](args.config, args.out)
    elif command == "train-sequence":
        from .sequence.training import train_sequence

        train_sequence(args.config, args.out)
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
    elif command == "variant-effects":
        from .sequence.effects import variant_effects_command

        variant_effects_command(args.config, args.out)
    else:
        raise PaceError(f"Command {command} is not implemented")
    return {"output": args.out, "command": command, "status": "complete"}


if __name__ == "__main__":
    raise SystemExit(main())
