"""Command-line interface: every command takes ordinary options, like bedtools or samtools."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from . import __version__
from .errors import PaceError
from .provenance import clean, output_policy
from .run_options import add_run_options, config_from_args

OVERVIEW = """\
PACE: enhancer-gene prediction for livestock and poultry (Activity-by-Contact family).

One command from peaks, GTF and bigWig files:
  predict        Build candidates, quantify signal, extract Hi-C and score in one step

Step by step:
  catalog        Candidate elements and gene TSSs from peaks (BED) and a GTF
  activity       ATAC/DNase/H3K27ac signal per element from a bigWig
  contacts       Element-TSS contacts from a .cool/.mcool Hi-C file
  merge          Combine activity or contact tables (replicates, assays)
  run            Score element-gene pairs from prepared tables
  validate       Check prepared tables without writing results

No or sparse Hi-C:
  fit-prior      Fit a distance-decay contact prior from a .cool/.mcool file
  boundaries     Call candidate CTCF boundaries from motifs or a genome FASTA
  prior          Package a contact prior from known parameters
  fit-hic        Fit distance decay, boundary attenuation and dispersion from contacts
  fuse           Score with per-pair Gamma-Poisson shrinkage of sparse contacts
  fit-labels     Tune gamma, beta and eta on eQTL fine-mapping (weak labels)

Downstream analysis:
  compare        Compare two runs (animals, treatments) on common denominators
  stability      Replicate consistency across several runs
  benchmark      Evaluate against functional labels (CRISPRi etc.) and baselines
  fit-eta        Learn the optional allocation exponent from functional labels
  train          Train a separate classifier on functional labels
  predict-ml     Apply a trained classifier to a run

Other data and helpers:
  features       Overlap a BED (CTCF, H3K4me1 peaks, motifs) with elements
  methylation    Summarize CpG methylation counts per element
  expression     Convert transcript TPM to gene expression
  pairs          List element-TSS pairs that need a contact value
  normalize-activity   Counts per million per bp from raw window counts
  promoter-weights     TSS weights from promoter signal (CAGE, ATAC, H3K4me3)
  init           Write empty input-table templates for a new project
  demo           Run the bundled synthetic example

Quick start - choose A, B or C by the Hi-C you have
(peaks, a GTF and at least one ATAC/DNase/H3K27ac bigWig are always required):

  A. Hi-C from this sample (best)
     pace predict -b peaks.bed -g genes.gtf --atac atac.bw --h3k27ac k27ac.bw \\
         --hic sample.mcool --hic-resolution 10000 \\
         --species pig --assembly Sscrofa11.1 --tissue liver -o A_out

  B. No Hi-C for this sample, but any Hi-C of the same species and assembly
     (another tissue, animal or public dataset): fit gamma once, then reuse it
     pace fit-prior --cooler other.mcool::/resolutions/10000 \\
         --species pig --assembly Sscrofa11.1 --tissue liver -o pig_prior
     pace predict -b peaks.bed -g genes.gtf --atac atac.bw \\
         --prior pig_prior \\
         --species pig --assembly Sscrofa11.1 --tissue liver -o B_out

  C. No Hi-C for this species at all: human ABC power law (gamma = 1.02),
     labelled "unvalidated" in the outputs
     pace predict -b peaks.bed -g genes.gtf --atac atac.bw --abc-prior \\
         --species pig --assembly Sscrofa11.1 --tissue liver -o C_out

  No parameter values need to be set by hand in A, B or C.
  Test the installation first:  pace demo -o demo_out

Run 'pace COMMAND --help' for the options of a command.
"""

RUN_COMMANDS = {
    "run": "Score element-gene pairs from prepared activity and contact tables",
    "validate": "Check prepared tables and settings without writing results",
    "fit-eta": "Learn the optional allocation exponent from independent functional labels",
    "fuse": "Score with per-pair Gamma-Poisson shrinkage of sparse Hi-C contacts",
    "capabilities": "Inspect available contact assets and validation scope",
}
# Former names keep working.
ALIASES = {
    "measured": "run",
    "merge-tables": "merge",
    "prepare-pairs": "pairs",
    "prepare-promoter-weights": "promoter-weights",
}
HIDDEN = {"prepare", "capabilities", "fit-contact-prior"}


class Formatter(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    """Show defaults next to each option and keep hand-written examples intact."""

    def _get_help_string(self, action):
        default = action.default
        if default is None or default is False or default == argparse.SUPPRESS:
            return action.help
        if action.help is None or action.help == argparse.SUPPRESS or action.nargs == 0:
            return action.help
        if "default" in action.help.lower() or action.option_strings == []:
            return action.help
        return action.help + " (default: %(default)s)"


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    shown_alias = argv[0] if argv and argv[0] in ALIASES else None
    if shown_alias:
        argv = [ALIASES[shown_alias], *argv[1:]]
    elif argv and argv[0].startswith("-") and argv[0] not in ("--help", "-h", "--version"):
        argv = ["run", *argv]
    parser = build_parser(shown_alias)
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


def build_parser(shown_alias=None):
    parser = argparse.ArgumentParser(
        prog="pace",
        allow_abbrev=False,
        usage="pace COMMAND [options]",
        description=OVERVIEW,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"PACE {__version__}")
    sub = parser.add_subparsers(dest="command", required=True, metavar="COMMAND", help="")
    sub.help = argparse.SUPPRESS

    def command(name, description, epilog=None):
        p = sub.add_parser(
            name,
            description=description,
            epilog=epilog,
            allow_abbrev=False,
            formatter_class=Formatter,
        )
        p.prog = f"pace {name}"
        if shown_alias and ALIASES.get(shown_alias) == name:
            p.prog = f"pace {shown_alias}"
        return p

    predict_parser(command)
    catalog_parser(command)
    activity_parser(command)
    contacts_parser(command)
    for name, description in RUN_COMMANDS.items():
        p = command(name, description, RUN_EPILOGS.get(name))
        add_run_options(p, output=name in ("run", "fit-eta", "fuse"))
    prior_parsers(command)
    analysis_parsers(command)
    helper_parsers(command)
    legacy = command("prepare", "Legacy: prepare inputs from a YAML description")
    legacy.add_argument("--config", required=True, help="YAML file with a 'kind' key")
    legacy.add_argument("-o", "--out", required=True, help="New output directory")
    for child in sub.choices.values():
        if any(a.dest == "out" for a in child._actions):
            child.add_argument(
                "--force",
                action="store_true",
                help="Replace an existing PACE output directory (the old one is kept as a backup)",
            )
    return parser


def add_out(p, help="New output directory"):
    p.add_argument("-o", "--out", required=True, metavar="DIR", help=help)


def add_legacy_config(p):
    p.add_argument("--config", help=argparse.SUPPRESS)


def add_context(p, *, required=True):
    group = p.add_argument_group("biological context")
    group.add_argument("--species", required=required, help="Species, e.g. pig, cattle, chicken")
    group.add_argument("--assembly", required=required, help="Genome assembly, e.g. Sscrofa11.1")
    group.add_argument("--tissue", required=required, help="Tissue or cell type, e.g. liver")
    group.add_argument(
        "--target-level",
        choices=["individual", "population_mean"],
        default="individual",
        help="One animal, or an equal-weight mean of several animals",
    )
    return group


# ---------------------------------------------------------------- predict


PREDICT_EPILOG = """\
Pick ONE contact source according to your data:

  A. --hic FILE    Hi-C of this sample (.cool/.mcool), best choice.
                   The distance power law of the same map is fitted
                   automatically and fills masked Hi-C bins.
     pace predict -b peaks.bed -g genes.gtf --atac atac.bw --h3k27ac k27ac.bw \\
         --hic pig1.mcool --hic-resolution 10000 \\
         --species pig --assembly Sscrofa11.1 --tissue liver -o pig1_liver

  B. --prior DIR   No Hi-C for this sample, but Hi-C of the same species and
                   assembly exists (other tissue/animal/public data).
                   Fit its gamma once with 'pace fit-prior', then reuse it.
     pace fit-prior --cooler other.mcool::/resolutions/10000 \\
         --species pig --assembly Sscrofa11.1 --tissue liver -o pig_prior
     pace predict -b peaks.bed -g genes.gtf --atac atac.bw --prior pig_prior \\
         --species pig --assembly Sscrofa11.1 --tissue liver -o pig1_liver
     (prior from another tissue: fit with that tissue's name, then score with
      'pace run --allow-cross-context-prior'; see the README)

  C. --abc-prior   No Hi-C for this species at all. Uses the human ABC power
                   law (gamma = 1.02); outputs are labelled
                   unvalidated_for_target_context.
     pace predict -b peaks.bed -g genes.gtf --atac atac.bw --abc-prior \\
         --species pig --assembly Sscrofa11.1 --tissue liver -o pig1_liver

No parameter values need to be set by hand in any of the three.
Several bigWigs for one assay are treated as replicates of the same animal.
The output folder holds the scores (scores.tsv.gz) and the prepared tables
(prepared/), so each step can be inspected or re-run with 'pace run'.
"""


def predict_parser(command):
    p = command(
        "predict",
        "Predict enhancer-gene links in one step: peaks + GTF + bigWig (+ Hi-C) -> scores",
        PREDICT_EPILOG,
    )
    inputs = p.add_argument_group("required inputs")
    inputs.add_argument(
        "-b",
        "--peaks",
        nargs="+",
        required=True,
        metavar="BED",
        help="Candidate regions, e.g. merged MACS peaks (several files are united)",
    )
    inputs.add_argument(
        "-g", "--gtf", required=True, help="Gene annotation with transcript lines (.gtf[.gz])"
    )
    add_out(p, "New output directory (results and prepared tables)")
    activity = p.add_argument_group("activity (at least one assay)")
    for assay in ("atac", "dnase", "h3k27ac"):
        activity.add_argument(
            "--" + assay, nargs="+", metavar="BIGWIG", help=f"{ASSAYS[assay]} bigWig file(s)"
        )
    activity.add_argument(
        "--min-callable-fraction",
        type=float,
        metavar="F",
        default=0.0,
        help="Elements with a smaller fraction of bases covered by the bigWig become NA",
    )
    activity.add_argument(
        "--missing-as-zero",
        action="store_true",
        help="Treat bases absent from the bigWig as measured zero signal",
    )
    contact = p.add_argument_group("contact (choose one)")
    contact.add_argument("--hic", metavar="COOL", help="Hi-C .cool or .mcool file")
    contact.add_argument(
        "--prior", metavar="DIR", help="Contact prior from 'pace fit-prior' (no Hi-C needed)"
    )
    contact.add_argument(
        "--abc-prior",
        action="store_true",
        help="Use the human ABC distance power law (unvalidated in livestock)",
    )
    contact.add_argument(
        "--hic-resolution", type=int, metavar="BP", help="Hi-C bin size, e.g. 5000 or 10000"
    )
    contact.add_argument(
        "--no-balance",
        action="store_true",
        help="Use raw Hi-C counts; the cooler has no balancing weights",
    )
    contact.add_argument(
        "--strict-contacts",
        action="store_true",
        help="With --hic: leave pairs in masked Hi-C bins NA instead of using the "
        "distance prior fitted on the same map (genes with such pairs are then withheld)",
    )
    add_context(p)
    cat = p.add_argument_group("candidate elements")
    add_catalog_options(cat, chrom_sizes_required=False)
    score = p.add_argument_group("scoring")
    score.add_argument("--run-id", default="pace", help="Name recorded in the outputs")
    score.add_argument(
        "--chunk-pairs",
        type=int,
        default=1_500_000,
        help="Score whole chromosomes in chunks of at most this many candidate pairs to "
        "bound memory; results equal a single run. 0 scores everything at once",
    )
    score.add_argument(
        "-t",
        "--threads",
        type=int,
        default=1,
        help="Score this many chromosome chunks in parallel (memory grows with each)",
    )
    score.add_argument(
        "--partial-policy",
        choices=["withhold", "conditional"],
        default="withhold",
        help="Genes with unscored candidates: withhold scores, or score the measured subset",
    )
    score.add_argument(
        "--tss-weights",
        choices=["equal", "provided"],
        default="equal",
        help="Weights of the distinct TSSs of a gene",
    )


ASSAYS = {"atac": "ATAC", "dnase": "DNase", "h3k27ac": "H3K27ac"}


# ---------------------------------------------------------------- prepare steps


def add_catalog_options(group, *, chrom_sizes_required=True):
    group.add_argument(
        "-c",
        "--chrom-sizes",
        required=chrom_sizes_required,
        metavar="FILE",
        help="Chromosome lengths: chrom.sizes, genome.fa.fai or any bigWig of the assembly"
        + ("" if chrom_sizes_required else " (default: read from the first bigWig)"),
    )
    group.add_argument(
        "-w", "--width", type=int, default=500, help="Element width in bp (fixed genome grid)"
    )
    group.add_argument("--offset", type=int, default=0, help="Grid offset in bp")
    group.add_argument(
        "-r",
        "--radius",
        type=int,
        default=5_000_000,
        help="Pair each gene with elements up to this distance (bp) from its TSS",
    )
    group.add_argument(
        "--gene-types",
        nargs="+",
        metavar="TYPE",
        help="Keep only these gene/transcript biotypes, e.g. protein_coding (Ensembl) or mRNA (NCBI)",
    )
    group.add_argument(
        "--no-promoters",
        dest="include_promoters",
        action="store_false",
        help="Do not add promoters of other genes as competing candidates",
    )
    group.add_argument(
        "--skip-unlisted-chroms",
        action="store_true",
        help="Ignore peaks and genes on chromosomes missing from the sizes file",
    )
    group.add_argument(
        "--aliases", metavar="TSV", help="Chromosome renaming table (columns: alias, canonical)"
    )


def catalog_parser(command):
    p = command(
        "catalog",
        "Build candidate elements and gene TSSs from peaks (BED) and a GTF",
        "example:\n  pace catalog -b peaks.bed -g genes.gtf -c genome.fa.fai -o prepared/catalog",
    )
    p.add_argument(
        "-b", "--peaks", nargs="+", required=True, metavar="BED", help="Candidate region BED(s)"
    )
    p.add_argument("-g", "--gtf", required=True, help="Gene annotation with transcript lines")
    p.add_argument(
        "--source-id",
        nargs="+",
        help="Name for each BED file in the outputs (default: the file name)",
    )
    add_catalog_options(p)
    add_out(p)
    add_legacy_config(p)


def activity_parser(command):
    p = command(
        "activity",
        "Quantify ATAC/DNase/H3K27ac signal per candidate element from a bigWig",
        "example:\n  pace activity -i pig1_atac.bw -a ATAC -d prepared/catalog -o prepared/pig1_atac",
    )
    p.add_argument("-i", "--bigwig", required=True, help="Signal track (.bw)")
    p.add_argument(
        "-a", "--assay", required=True, choices=["ATAC", "DNase", "H3K27ac"], help="Assay"
    )
    where = p.add_mutually_exclusive_group(required=True)
    where.add_argument("-d", "--catalog-dir", help="Output folder of 'pace catalog'")
    where.add_argument("-u", "--units", help="units.tsv of a catalog")
    p.add_argument("-s", "--sample-id", help="Sample name (default: bigWig file name)")
    p.add_argument(
        "--normalization-id",
        default="bigwig_signal",
        help="Label for how the bigWig was normalized; use the same label for comparable files",
    )
    p.add_argument("--unit", default="normalized_signal", help="Label for the signal unit")
    p.add_argument("--window-id", help="Label for the summary window (default: grid:WIDTH:mean)")
    p.add_argument(
        "--missing-as-zero",
        action="store_true",
        help="Treat bases absent from the bigWig as measured zero signal",
    )
    p.add_argument(
        "--min-callable-fraction",
        type=float,
        metavar="F",
        default=0.0,
        help="Elements with a smaller covered fraction become NA (0.8 is a common choice)",
    )
    add_out(p)
    add_legacy_config(p)


def contacts_parser(command):
    p = command(
        "contacts",
        "Extract element-TSS contact frequencies from a .cool/.mcool Hi-C file",
        "example:\n  pace contacts -i pig1.mcool -r 10000 -d prepared/catalog -o prepared/pig1_hic\n\n"
        "Use contact frequencies (balanced or raw), never observed/expected or p-values.",
    )
    p.add_argument("-i", "--hic", required=True, help=".cool, .mcool or cooler URI (file::group)")
    p.add_argument("-r", "--resolution", type=int, required=True, help="Bin size in bp")
    where = p.add_mutually_exclusive_group(required=True)
    where.add_argument("-d", "--catalog-dir", help="Output folder of 'pace catalog'")
    where.add_argument("--pairs", help="pairs.tsv from 'pace pairs'")
    p.add_argument("-s", "--sample-id", help="Sample name (default: Hi-C file name)")
    p.add_argument("--source-id", help="Source name (default: the sample name)")
    p.add_argument(
        "--no-balance", action="store_true", help="Use raw counts instead of balanced contacts"
    )
    p.add_argument(
        "--missing-as-missing",
        action="store_true",
        help="Treat pixels absent from the matrix as unmeasured instead of zero contacts",
    )
    p.add_argument("--scale", help="Contact scale label (default: balanced_contact or raw_contact)")
    p.add_argument(
        "--normalization-id", help="Normalization label (default: cooler_weight or none)"
    )
    add_out(p)
    add_legacy_config(p)


# ---------------------------------------------------------------- priors


def add_asset_options(p, *, measurement=True, measurement_required=True):
    # Requirements are checked after parsing so that legacy YAML input keeps working.
    add_context(p, required=False)
    p.add_argument("--model-id", default="contact_prior", help="Name of the prior")
    p.add_argument(
        "--synthetic", action="store_true", help="Mark the asset as synthetic (demonstration)"
    )
    if measurement:
        m = p.add_argument_group("contact measurement definition")
        need = " (required)" if measurement_required else ""
        m.add_argument("--scale", help="Contact scale label" + need)
        m.add_argument("--resolution", type=int, help="Bin size in bp" + need)
        m.add_argument("--normalization-id", help="Normalization label")
        m.add_argument("--balancing", help="Balancing label, e.g. balanced or unbalanced")
        m.add_argument("--window-id", default="bin_pair", help="Contact window label")


def prior_parsers(command):
    p = command(
        "fit-prior",
        "Fit a distance-decay contact prior from a .cool/.mcool file (valid zero pixels included)",
        "example:\n  pace fit-prior --cooler pig_liver.mcool::/resolutions/10000 "
        "--species pig --assembly Sscrofa11.1 --tissue liver -o pig_prior",
    )
    p.add_argument("--cooler", help="cool file or mcool::/resolutions/N URI")
    for flag in ("species", "assembly", "tissue"):
        p.add_argument("--" + flag)
    p.add_argument(
        "--target-level", choices=["individual", "population_mean"], default="individual"
    )
    p.add_argument("--balanced", action=argparse.BooleanOptionalAction, default=True)
    p.add_argument("--min-distance", type=int, help="Default: one contact bin")
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
    add_out(p)
    add_legacy_config(p)

    p = command(
        "fit-contact-prior",
        "Fit a binned distance prior from a table of measured contacts",
    )
    p.add_argument("--contacts", dest="data", required=False, help="Contact table")
    p.add_argument("--bin-edges", nargs="+", type=float, help="Distance bin edges in bp")
    p.add_argument("--d-ref", type=float, help="Reference distance in bp")
    p.add_argument("--d-min", type=float, help="Minimum distance in bp")
    add_asset_options(p, measurement_required=False)
    add_out(p)
    add_legacy_config(p)

    p = command(
        "boundaries",
        "Call candidate CTCF boundaries from motif hits or by scanning a genome FASTA",
        "examples:\n  pace boundaries --motifs ctcf_motifs.tsv --chip ctcf_peaks.bed -o boundaries\n"
        "  pace boundaries --fasta genome.fa --pwm ctcf_pwm.tsv --threshold 12 -o boundaries",
    )
    src = p.add_mutually_exclusive_group()
    src.add_argument("--motifs", help="Motif table: chrom start end strand strength")
    src.add_argument("--fasta", help="Genome FASTA to scan with --pwm")
    p.add_argument("--pwm", help="CTCF PWM table with columns A C G T")
    p.add_argument("--threshold", type=float, help="Log2-odds motif score threshold")
    p.add_argument("--chip", help="CTCF ChIP peaks (BED with header chrom start end)")
    p.add_argument("--max-gap", dest="max_gap_bp", type=int, help="Maximum motif gap in bp")
    add_out(p)
    add_legacy_config(p)

    for name, fit in (("prior", False), ("fit-hic", True)):
        p = command(
            name,
            "Fit distance decay, CTCF boundary attenuation and dispersion from Hi-C counts"
            if fit
            else "Package a contact prior a*(d/d_ref)^-gamma*exp(-beta*boundaries) from known values",
        )
        if fit:
            p.add_argument("--contacts", dest="data", help="Contact table from 'pace contacts'")
            p.add_argument("--gamma-grid", nargs="+", type=float, help="Candidate gamma values")
            p.add_argument("--beta-grid", nargs="+", type=float, help="Candidate beta values")
            p.add_argument("--test-chromosomes", nargs="+", help="Chromosomes held out from fit")
        else:
            p.add_argument("--a", type=float, help="Amplitude a")
            p.add_argument("--gamma", type=float, help="Power-law exponent gamma")
            p.add_argument("--beta", type=float, default=0.0, help="Boundary attenuation beta")
            p.add_argument("--kappa", type=float, help="Gamma-Poisson concentration")
        p.add_argument("--d-ref", type=float, help="Reference distance in bp")
        p.add_argument("--d-min", type=float, help="Distances below this are set to it (bp)")
        p.add_argument("--boundaries", help="boundaries.tsv from 'pace boundaries'")
        p.add_argument("--pairs", help="Optional pairs.tsv to evaluate the prior on")
        add_asset_options(p, measurement_required=not fit)
        add_out(p)
        add_legacy_config(p)

    p = command(
        "fit-labels",
        "Tune gamma, beta and eta against eQTL fine-mapping PIPs with held-out chromosomes",
        "example:\n  pace fit-labels --run results/pig1 --eqtl eqtl_pip.tsv "
        "--test-chromosomes 17 18 -o weak_fit",
    )
    p.add_argument("--run", dest="run_config", help="A 'pace run' output folder (or its YAML)")
    p.add_argument(
        "--eqtl", dest="labels", help="Variant PIP table: variant_id chrom pos0 gene_id pip"
    )
    p.add_argument("--gamma-grid", nargs="+", type=float, default=[1.0])
    p.add_argument("--beta-grid", nargs="+", type=float, default=[0.0])
    p.add_argument("--eta-grid", nargs="+", type=float, default=[0.0, 0.5, 1.0])
    p.add_argument("--test-chromosomes", nargs="+", help="Chromosomes kept for final evaluation")
    p.add_argument(
        "--aggregation",
        choices=["independent_variants", "independent_signals"],
        default="independent_signals",
        help="How PIPs are combined per gene",
    )
    p.add_argument("--abc-gamma", type=float, help="ABC baseline gamma")
    add_out(p)
    add_legacy_config(p)


# ---------------------------------------------------------------- analysis


def analysis_parsers(command):
    p = command(
        "compare",
        "Compare two runs (animals, tissues, treatments) on their common candidate set",
        "example:\n  pace compare --left results/pig1 --right results/pig2 -o pig2_vs_pig1",
    )
    p.add_argument("--left", help="Reference run folder")
    p.add_argument("--right", help="Run folder compared with --left (differences are right-left)")
    p.add_argument("--min-common-units", dest="minimum_common_units", type=int, default=2)
    p.add_argument("--allow-eta-difference", action="store_true")
    p.add_argument("--allow-evidence-difference", action="store_true")
    add_out(p)
    add_legacy_config(p)

    p = command("stability", "Measure consistency of top predictions across replicate runs")
    p.add_argument(
        "--replicates", help="Table: run_path donor_id replicate_type (biological/technical)"
    )
    p.add_argument("--min-common-units", dest="minimum_common_units", type=int, default=2)
    add_out(p)
    add_legacy_config(p)

    p = command(
        "benchmark",
        "Evaluate a run against functional labels (e.g. CRISPRi) and simple baselines",
    )
    p.add_argument("--run", dest="run_config", help="A 'pace run' output folder (or its YAML)")
    p.add_argument("--labels", help="Functional label table")
    p.add_argument("--membership", dest="region_membership", help="region_membership.tsv")
    p.add_argument("--stratify", nargs="+", help="Label columns to report separately")
    p.add_argument(
        "--external",
        nargs=4,
        action="append",
        metavar=("NAME", "PATH", "VERSION", "SETTINGS"),
        help="Another method's scores (element_id gene_id score); repeatable",
    )
    p.add_argument(
        "--threshold",
        nargs=4,
        action="append",
        metavar=("METHOD", "VALUE", "SPLIT", "SOURCE"),
        help="Pre-registered threshold from a train/calibration split; repeatable",
    )
    add_out(p)
    add_legacy_config(p)

    p = command("train", "Train a separate classifier on functional labels")
    p.add_argument("--data", help="Training table")
    p.add_argument("--model-id", default="pace_classifier")
    for flag in ("species", "assembly", "tissue"):
        p.add_argument("--" + flag)
    p.add_argument("--synthetic", action="store_true", help="Synthetic training data")
    p.add_argument("--feature-contract", help="ml_feature_contract.json of the run")
    p.add_argument("--extra-features", nargs="+", help="Additional feature names")
    p.add_argument(
        "--penalty",
        nargs=2,
        type=float,
        action="append",
        metavar=("L1", "L2"),
        help="Elastic-net penalties to try; repeatable",
    )
    p.add_argument("--folds", type=int, default=3)
    p.add_argument("--seed", type=int, default=17)
    p.add_argument("--calibrate", action="store_true", help="Fit a probability calibration")
    add_out(p)
    add_legacy_config(p)

    p = command("predict-ml", "Apply a trained classifier to a run without refitting")
    p.add_argument("--run", help="A 'pace run' output folder")
    p.add_argument("--model", help="Output folder of 'pace train'")
    p.add_argument("--features", help="Optional extra features table")
    add_out(p)
    add_legacy_config(p)


# ---------------------------------------------------------------- helpers


def helper_parsers(command):
    p = command(
        "merge",
        "Combine tables of the same type (replicates, assays) with duplicate checks",
        "example:\n  pace merge -t observed_activity -i atac/observed_activity.tsv "
        "k27ac/observed_activity.tsv -o prepared/activity",
    )
    p.add_argument(
        "-t",
        "--table",
        choices=["observed_activity", "observed_contacts", "expression", "methylation"],
        required=True,
    )
    p.add_argument("-i", "--inputs", nargs="+", required=True)
    add_out(p)

    p = command("pairs", "List element-TSS pairs that need a contact value")
    p.add_argument("-d", "--catalog-dir", required=True)
    add_out(p)

    p = command("features", "Overlap a BED file (CTCF, H3K4me1, motifs) with candidate elements")
    p.add_argument("-b", "--bed", help="BED file with a header line chrom start end")
    where = p.add_mutually_exclusive_group()
    where.add_argument("-d", "--catalog-dir")
    where.add_argument("-u", "--units")
    p.add_argument("--source-id", help="Default: BED file name")
    p.add_argument("--evidence-id", help="Default: the source id")
    p.add_argument("--feature-prefix", help="Feature name prefix, e.g. CTCF")
    p.add_argument("--entity-type", choices=["element", "promoter"], default="element")
    p.add_argument("--motif-strands", action="store_true", help="Record motif strand counts")
    add_out(p)
    add_legacy_config(p)

    p = command("methylation", "Summarize CpG methylation counts per element")
    p.add_argument("--counts", help="CpG table: chrom dyad_start0 methylated_count total_count ...")
    where = p.add_mutually_exclusive_group()
    where.add_argument("-d", "--catalog-dir")
    where.add_argument("-u", "--units")
    p.add_argument("--min-coverage", dest="minimum_coverage", type=int, default=1)
    p.add_argument("--reference-cpg", help="Optional element_id n_cpg table")
    add_out(p)
    add_legacy_config(p)

    p = command("expression", "Convert transcript TPM to gene-level expression")
    p.add_argument("--tpm", dest="expression", help="Table: transcript_id sample_id tpm")
    where = p.add_mutually_exclusive_group()
    where.add_argument("-d", "--catalog-dir", help="Uses its transcript_mapping.tsv")
    where.add_argument("--transcript-mapping")
    add_out(p)
    add_legacy_config(p)

    p = command(
        "normalize-activity", "Normalize raw window counts by library size and window length"
    )
    for flag in ("counts", "library-sizes", "units"):
        p.add_argument("--" + flag, required=True)
    p.add_argument("--normalization-id", default="CPM_density_v1")
    p.add_argument("--window-id")
    add_out(p)

    p = command("promoter-weights", "Freeze TSS weights from measured promoter signals")
    for flag in ("promoters", "signals", "normalization-id"):
        p.add_argument("--" + flag, required=True)
    p.add_argument("--assay", choices=["ATAC", "DNase", "H3K4me3", "CAGE"], required=True)
    p.add_argument("--zero-policy", choices=["error", "equal"], default="error")
    add_out(p)

    p = command("init", "Write a project folder with empty input-table templates")
    add_context(p)
    p.add_argument("--catalog-dir")
    p.add_argument(
        "--panel", nargs="+", choices=["ATAC", "DNase", "H3K27ac"], default=["ATAC", "H3K27ac"]
    )
    p.add_argument("--contact-prior")
    p.add_argument("--prior-preset", choices=["abc_human"])
    add_out(p)

    p = command("demo", "Run the bundled synthetic example (tests the installation)")
    p.add_argument("--regime", choices=["measured"], default="measured", help=argparse.SUPPRESS)
    add_out(p)


RUN_EPILOGS = {
    "run": """\
examples:
  pace run -d prepared/catalog --activity prepared/activity/observed_activity.tsv \\
      --contacts prepared/pig1_hic/observed_contacts.tsv \\
      --species pig --assembly Sscrofa11.1 --tissue liver -o results/pig1_liver

  # several animals averaged: give a samples table with donor_id per sample
  pace run -d prepared/catalog --activity all_activity.tsv --samples samples.tsv \\
      --target-level population_mean --species pig --assembly Sscrofa11.1 \\
      --tissue liver --contact-mode prior_only --contact-prior pig_prior -o results/liver

Without --samples, every sample is taken as a replicate of one animal.
""",
}


# ---------------------------------------------------------------- dispatch


def options(args, *names, **renamed):
    """Collect argparse values into an operation mapping; None means 'use the default'."""
    result = {name: getattr(args, name) for name in names}
    result.update({key: getattr(args, name) for key, name in renamed.items()})
    return result


def legacy_or(args, mapping):
    if getattr(args, "config", None):
        return args.config
    return mapping


def units_from(args):
    if getattr(args, "catalog_dir", None):
        return str(Path(args.catalog_dir) / "units.tsv")
    return getattr(args, "units", None)


def asset_options(args):
    return {
        "species": args.species,
        "assembly": args.assembly,
        "context_id": args.tissue,
        "target_level": args.target_level,
        "model_id": args.model_id,
        "is_synthetic": bool(args.synthetic),
        "scale": args.scale,
        "resolution": args.resolution,
        "normalization_id": args.normalization_id,
        "balancing": args.balancing,
        "window_id": args.window_id,
    }


def stem(path) -> str:
    name = Path(str(path).split("::", 1)[0]).name
    for suffix in (".gz", ".bw", ".bigwig", ".bigWig", ".mcool", ".cool", ".bed", ".tsv"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    return name


def dispatch(args):
    command = args.command
    if command == "predict":
        from .workflow import predict

        return predict(args)
    if command in ("catalog", "activity", "contacts", "features", "methylation", "expression"):
        from .operations import prepare_command

        return {
            "output": args.out,
            **prepare_command(legacy_or(args, prepare_options(args)), args.out),
        }
    if command == "prepare":
        from .operations import prepare_command

        prepare_command(args.config, args.out)
        return {"output": args.out, "command": command, "status": "complete"}
    if command in ("boundaries", "fit-hic", "prior", "fit-labels"):
        from .contact_cli import contact_command

        return contact_command(command, legacy_or(args, contact_options(args)), args.out)
    if command in ("merge", "pairs", "normalize-activity", "promoter-weights", "init"):
        from .preparation import (
            init_project,
            merge_tables,
            normalize_activity,
            prepare_pairs,
            promoter_weights,
        )

        return {
            "init": init_project,
            "pairs": prepare_pairs,
            "merge": merge_tables,
            "normalize-activity": normalize_activity,
            "promoter-weights": promoter_weights,
        }[command](args)
    if command == "fit-prior":
        if args.config:
            if args.cooler:
                raise PaceError("Choose --config or --cooler, not both")
            from .operations import fit_contact_command

            fit_contact_command(args.config, args.out)
            return {"output": args.out}
        if not all((args.cooler, args.species, args.assembly, args.tissue)):
            raise PaceError("fit-prior needs --cooler, --species, --assembly and --tissue")
        from .io.prior_fit import fit_cooler_command

        return fit_cooler_command(args)
    if command == "fit-contact-prior":
        from .operations import fit_contact_command

        mapping = {**asset_options(args), **options(args, "data", "bin_edges", "d_ref", "d_min")}
        fit_contact_command(legacy_or(args, mapping), args.out)
        return {"output": args.out, "command": command, "status": "complete"}
    if command == "demo":
        from .demo import demo

        return demo(args.regime, args.out)
    if command in RUN_COMMANDS:
        return run_command(args)
    if command in ("compare", "stability"):
        from .evaluation.compare import compare_command, stability_command

        if command == "compare":
            mapping = options(
                args,
                "left",
                "right",
                "minimum_common_units",
                "allow_eta_difference",
                "allow_evidence_difference",
            )
            compare_command(legacy_or(args, mapping), args.out)
        else:
            mapping = options(args, "replicates", "minimum_common_units")
            stability_command(legacy_or(args, mapping), args.out)
    elif command == "train":
        from .learning.model import train_classifier

        train_classifier(legacy_or(args, train_options(args)), args.out)
    elif command == "predict-ml":
        from .evaluation.benchmark import predict_ml_command

        predict_ml_command(legacy_or(args, options(args, "run", "model", "features")), args.out)
    elif command == "benchmark":
        from .evaluation.benchmark import benchmark_command

        mapping = options(args, "run_config", "labels", "region_membership")
        mapping["stratify"] = args.stratify
        if args.external:
            mapping["external_methods"] = [
                dict(zip(("name", "path", "version", "configuration"), e, strict=True))
                for e in args.external
            ]
        if args.threshold:
            mapping["thresholds"] = {
                method: {"value": float(value), "source_split": split, "source_id": source}
                for method, value, split, source in args.threshold
            }
        benchmark_command(legacy_or(args, mapping), args.out)
    else:
        raise PaceError(f"Command {command} is not implemented")
    return {"output": args.out, "command": command, "status": "complete"}


def run_command(args):
    from .evidence.assets import capabilities
    from .pipeline import compute, run

    command = args.command
    if command == "fuse":
        if args.contact_mode not in (None, "shrinkage") or args.contact_reliability not in (
            None,
            "per_pair",
        ):
            raise PaceError(
                "fuse always uses --contact-mode shrinkage --contact-reliability per_pair"
            )
        args.contact_mode, args.contact_reliability = "shrinkage", "per_pair"
    cfg = config_from_args(args)
    if command == "fuse" and (
        cfg["contact"]["mode"] != "shrinkage" or cfg["contact"]["reliability"] != "per_pair"
    ):
        raise PaceError("fuse needs contact.mode=shrinkage and reliability=per_pair")
    if command == "fit-eta" and not cfg["allocation"]["labels_path"]:
        raise PaceError("fit-eta requires --eta-labels")
    if command == "capabilities":
        return capabilities(cfg)
    if command == "validate":
        result = compute(cfg)
        return {"valid": True, "n_candidates": len(result["scores"]), "qc": result["qc"]}
    if getattr(args, "by_chromosome", False):
        from .chunked import run_by_chromosome

        result = run_by_chromosome(cfg, args.out, max_pairs=args.chunk_pairs, threads=args.threads)
        return {
            "output": args.out,
            "n_candidates": result["qc"]["n_candidates"],
            "n_scoreable": result["qc"]["n_scoreable"],
            "chunks": len(result["qc"]["execution"]["chunks"]),
            "eta": result["eta_calibration"]["eta"],
            "eta_status": result["eta_calibration"]["status"],
        }
    result = run(cfg, args.out)
    return {
        "output": args.out,
        "n_candidates": len(result["scores"]),
        "n_scoreable": result["qc"]["n_scoreable"],
        "eta": result["eta_calibration"]["eta"],
        "eta_status": result["eta_calibration"]["status"],
    }


def prepare_options(args):
    command = args.command
    if command == "catalog":
        return {
            "kind": "catalog",
            "bed": args.peaks,
            "gtf": args.gtf,
            "chrom_sizes": args.chrom_sizes,
            "source_id": args.source_id,
            "width": args.width,
            "offset": args.offset,
            "radius": args.radius,
            "include_promoters": args.include_promoters,
            "aliases": args.aliases,
            "gene_types": args.gene_types,
            "skip_unlisted_chroms": args.skip_unlisted_chroms,
        }
    if command == "activity":
        units = units_from(args)
        window = args.window_id
        if window is None:
            window = default_window_id(units)
        return {
            "kind": "bigwig",
            "track": args.bigwig,
            "units": units,
            "sample_id": args.sample_id or stem(args.bigwig),
            "assay": args.assay,
            "unit": args.unit,
            "normalization_id": args.normalization_id,
            "window_id": window,
            "missing_is_measured_zero": args.missing_as_zero,
            "minimum_callable_fraction": args.min_callable_fraction,
        }
    if command == "contacts":
        balanced = not args.no_balance
        sample = args.sample_id or stem(args.hic)
        return {
            "kind": "cooler",
            "contact": cooler_uri(args.hic, args.resolution),
            "pairs": args.pairs,
            "catalog_dir": args.catalog_dir,
            "resolution": args.resolution,
            "balanced": balanced,
            "missing_pixels_are_zero": not args.missing_as_missing,
            "sample_id": sample,
            "source_id": args.source_id or sample,
            "scale": args.scale or ("balanced_contact" if balanced else "raw_contact"),
            "normalization_id": args.normalization_id or ("cooler_weight" if balanced else "none"),
        }
    if command == "features":
        source = args.source_id or (stem(args.bed) if args.bed else None)
        return {
            "kind": "bed_features",
            "bed": args.bed,
            "units": units_from(args),
            "source_id": source,
            "evidence_id": args.evidence_id or source,
            "feature_prefix": args.feature_prefix or source,
            "entity_type": args.entity_type,
            "motif_strands": args.motif_strands,
        }
    if command == "methylation":
        return {
            "kind": "methylation",
            "counts": args.counts,
            "units": units_from(args),
            "minimum_coverage": args.minimum_coverage,
            "reference_cpg": args.reference_cpg,
        }
    mapping = args.transcript_mapping
    if args.catalog_dir:
        mapping = str(Path(args.catalog_dir) / "transcript_mapping.tsv")
    return {"kind": "rna", "expression": args.expression, "transcript_mapping": mapping}


def default_window_id(units):
    """grid:WIDTH:mean when every element has one width, else provided_regions:mean."""
    from .io.tables import read_table

    if not units or not Path(units).is_file():
        return None
    widths = {int(r["end"]) - int(r["start"]) for r in read_table(units, required=["start", "end"])}
    return f"grid:{widths.pop()}:mean" if len(widths) == 1 else "provided_regions:mean"


def cooler_uri(path, resolution):
    """Point an .mcool file at the requested resolution; other URIs are used as given."""
    if path and "::" not in path and path.endswith(".mcool") and resolution:
        return f"{path}::/resolutions/{resolution}"
    return path


def contact_options(args):
    command = args.command
    if command == "boundaries":
        return options(args, "motifs", "fasta", "pwm", "threshold", "chip", "max_gap_bp")
    if command == "fit-labels":
        mapping = options(
            args,
            "run_config",
            "labels",
            "gamma_grid",
            "beta_grid",
            "eta_grid",
            "test_chromosomes",
            "aggregation",
            "abc_gamma",
        )
        return mapping
    mapping = {**asset_options(args), **options(args, "d_ref", "d_min", "boundaries", "pairs")}
    if command == "fit-hic":
        mapping.update(options(args, "data", "gamma_grid", "beta_grid", "test_chromosomes"))
        infer_measurement(mapping)
    else:
        mapping.update(options(args, "a", "gamma", "beta", "kappa"))
    return mapping


def infer_measurement(mapping):
    """Take unset scale/resolution/normalization labels from the contact table itself."""
    from .io.tables import read_table

    if not mapping.get("data") or not Path(mapping["data"]).is_file():
        return
    rows = read_table(mapping["data"])
    for key in ("scale", "resolution", "normalization_id", "balancing", "window_id"):
        if mapping.get(key) is None or (key == "window_id" and mapping[key] == "bin_pair"):
            values = {r.get(key) for r in rows} - {None}
            if len(values) == 1:
                value = values.pop()
                mapping[key] = int(value) if key == "resolution" else value


def train_options(args):
    contract = None
    if args.feature_contract:
        contract = json.loads(Path(args.feature_contract).read_text(encoding="utf-8"))
    context = None
    if args.species or args.assembly or args.tissue:
        context = {"species": args.species, "assembly": args.assembly, "context_id": args.tissue}
    return {
        "data": args.data,
        "model_id": args.model_id,
        "is_synthetic": bool(args.synthetic),
        "context": context,
        "feature_contract": contract,
        "extra_features": args.extra_features,
        "penalties": [list(p) for p in args.penalty] if args.penalty else None,
        "folds": args.folds,
        "seed": args.seed,
        "calibrate": args.calibrate,
    }


if __name__ == "__main__":
    raise SystemExit(main())
