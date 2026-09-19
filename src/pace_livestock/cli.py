"""Thin command-line layer: scientific errors return nonzero, outputs are transactional."""

from __future__ import annotations

import argparse
import json
import sys

from . import __version__
from .errors import PaceError
from .provenance import clean
from .run_options import MODES, add_run_options, config_from_args

COMMANDS = {
    "validate": "Validate schemas, provenance, assets and scientific input contracts",
    "capabilities": "Inspect implementation, weights and validation scope",
    "prepare": "Convert genomic files into canonical standard tables",
    "run": "Compute the canonical PACE score in the configured evidence regime",
    "fit-eta": "Calibrate eta from functional labels and export scores plus a reusable artifact",
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

COMMON_COMMANDS = {"demo", "run", "validate", "capabilities"}


def print_user_help() -> None:
    """Print the task-oriented first-run guide, without exposing internals."""
    print(
        """PACE — score regulatory candidates from three kinds of evidence

Choose the kind of data you have:

  PACE measured --help     Measured activity/contact data
  PACE hybrid --help       Measured data plus sequence/genome predictions
  PACE genome --help       A reference genome plus variants for one individual

Try it safely first (choose the mode and output directory explicitly):

  PACE demo --regime measured --out demo-results

The usual workflow is:

  1. choose measured, hybrid, or genome
  2. run that mode's --help to see only its inputs
  3. provide --config FILE (recommended) and --out DIRECTORY

Examples:

  PACE measured --config examples/measured/config.yaml --out results
  PACE hybrid --config examples/hybrid/config.yaml --out results
  PACE genome --config examples/genome/config.yaml --out results

Use PACE --version for the installed version.
"""
    )


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv in (["--help"], ["-h"]):
        print_user_help()
        return 0
    requested_mode = argv[0] if argv and argv[0] in MODES else None
    if argv and argv[0] in MODES:
        argv = ["run", "--mode", argv[0], *argv[1:]]
    elif argv and argv[0].startswith("-") and argv[0] not in ("--help", "-h", "--version"):
        argv = ["run", *argv]
    parser = argparse.ArgumentParser(
        prog=f"PACE {requested_mode}" if requested_mode else "PACE",
        description=(
            "PACE — score candidate regulatory elements from measured, hybrid, or genome data.\n\n"
            "Start here:\n"
            "  PACE demo                         Run a safe offline example\n"
            "  PACE run --help                  Score your own prepared data\n"
            "  PACE validate --help             Check inputs before a run\n"
            "  PACE capabilities --help         See supported evidence and limits\n\n"
            "Outputs are written to a new directory; existing results are never overwritten."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  PACE demo --out demo-results\n"
            "  PACE measured --config examples/measured/config.yaml --out results\n"
            "\n"
            "For advanced workflows, run PACE <command> --help."
        ),
    )
    parser.add_argument("--version", action="version", version=f"PACE {__version__}")
    sub = parser.add_subparsers(dest="command", required=True, metavar="COMMAND")
    for name, description in COMMANDS.items():
        label = ("常用：" if name in COMMON_COMMANDS else "高级：") + description
        p = sub.add_parser(
            name,
            help=label,
            description=description,
            formatter_class=argparse.RawDescriptionHelpFormatter if requested_mode and name == "run" else argparse.HelpFormatter,
        )
        if name in ("run", "fit-eta", "validate", "capabilities"):
            if requested_mode and name == "run":
                examples = {
                    "measured": (
                        "Measured mode: score candidates using observed activity and contact tables.\n\n"
                        "Recommended (config file):\n"
                        "  PACE measured --config examples/measured/config.yaml --out results\n\n"
                        "Common direct inputs (one animal/tissue):\n"
                        "  PACE measured --catalog-dir data "
                        "--activity animal_001_liver_ATAC_H3K27ac.tsv "
                        "--contacts animal_001_liver_HiC_5kb.tsv --out results\n"
                        "  activity = ATAC-seq/H3K27ac; contacts = Hi-C or Prom-Hi-C"
                    ),
                    "hybrid": (
                        "Hybrid mode: combine measured evidence with sequence/genome predictions.\n\n"
                        "Recommended (config file):\n"
                        "  PACE hybrid --config examples/hybrid/config.yaml --out results\n\n"
                        "Common direct inputs (one animal/tissue):\n"
                        "  PACE hybrid --catalog-dir data "
                        "--activity animal_001_liver_ATAC_H3K27ac.tsv "
                        "--contacts animal_001_liver_HiC_5kb.tsv "
                        "--reference cattle_ARS-UCD1.2.fa --vcf animal_001.vcf.gz "
                        "--sequence-model models/sequence --out results\n"
                        "  optional annotation: animal_001_liver_RNAseq.tsv"
                    ),
                    "genome": (
                        "Genome mode: score one individual from reference, variants, and callable sites.\n\n"
                        "Recommended (config file):\n"
                        "  PACE genome --config examples/genome/config.yaml --out results\n\n"
                        "Common direct inputs (one animal):\n"
                        "  PACE genome --catalog-dir data "
                        "--reference cattle_ARS-UCD1.2.fa --vcf animal_001.vcf.gz "
                        "--callable animal_001_callable.bed --out results\n"
                        "  genome inputs = FASTA + VCF/BCF + callable sites; no ATAC/RNA file is required"
                    ),
                }[requested_mode]
                p.description = examples
                p.prog = f"PACE {requested_mode}"
            add_run_options(p, output=name in ("run", "fit-eta"), help_mode=requested_mode)
            continue
        p.add_argument(
            "--config", required=True, help="YAML file; input paths resolve relative to this file"
        )
        if name not in ("validate", "capabilities"):
            p.add_argument("--out", required=True, help="New output directory (never overwritten)")
    d = sub.add_parser(
        "demo",
        help="常用：Run a fully offline example (recommended first command)",
        description="Run a fully offline example using synthetic data. This does not need input files.",
    )
    d.add_argument(
        "--regime",
        choices=["measured", "hybrid", "genome_only"],
        required=True,
        help="Evidence mode to demonstrate",
    )
    d.add_argument("--out", required=True, help="New output directory")
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
