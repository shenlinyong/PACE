# PACE

[![Software tests](https://github.com/shenlinyong/PACE/actions/workflows/ci.yml/badge.svg)](https://github.com/shenlinyong/PACE/actions/workflows/ci.yml)

**Activity–contact modelling for livestock enhancer–gene research.**

PACE combines quantitative regulatory activity and contact to a gene's physical
promoters into a relative support score. It supports measured, hybrid and
sequence-only evidence with explicit quality, missingness and provenance checks.
Author and maintainer: **申林用 (Linyong Shen, shenlinyong), Northwest A&F University**.

[中文完整使用说明](docs/USER_GUIDE.zh-CN.md) · [中文首页](README.zh-CN.md) ·
[Installation: Conda / pip / Docker](docs/INSTALLATION.md) ·
[Three-mode tutorial](docs/TUTORIAL.md) · [Equations](docs/```math
\boxed{
\mathrm{PACE}(E,G)=
\frac{A_\star(E)\,\overline C(E,G)\,[B(E,G)]^{\eta_{\mathrm{used}}}}
{\displaystyle\sum_{e\in\mathcal E^{\mathrm{score}}(G)}
 A_\star(e)\,\overline C(e,G)\,[B(e,G)]^{\eta_{\mathrm{used}}}}
}
```.md)

## Download and install

The commands below run on a Linux server with Git and Conda installed:

```bash
git clone https://github.com/shenlinyong/PACE.git
cd PACE
conda env create -f environment.yml
conda activate pace
PACE --version
```

Alternatively, use Python >=3.11 in an isolated environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install '.[io,ml]'
```

Or build the included Docker image locally:

```bash
docker build -t pace:local .
docker run --rm pace:local --version
```

Real CNN inference/training requires the optional `sequence` dependencies. The
[installation guide](docs/INSTALLATION.md) includes wget source downloads, CPU
PyTorch, Docker volume mounts and development installation. No published PyPI
package, registry image or validated livestock weight download is assumed.

## Run the bundled examples

From the repository root, after installation:

```bash
PACE measured --config examples/measured/config.yaml --out results/measured
PACE hybrid --config examples/hybrid/config.yaml --out results/hybrid
PACE genome --config examples/genome_only/config.yaml --out results/genome
```

These are complete **synthetic** fixtures that run offline without a GPU. They
check installation and demonstrate the workflow; they are not biological models.
Choose a new output directory each time. `PACE`, `pace` and `pace-livestock` are
aliases for the same implementation.

For your own data, read the [step-by-step tutorial](docs/TUTORIAL.md),
[three complete project configurations](docs/USER_GUIDE.zh-CN.md), and
[input preparation guide](docs/input_preparation.md). Paths in YAML resolve from
the YAML file's directory; direct CLI paths resolve from the working directory.

## Choose the evidence mode

| Mode | Activity | Contact | Requirements beyond the candidate catalog |
|---|---|---|---|
| `measured` | Qualified observed signals in a fixed assay panel | Observations, or an explicitly configured applicable prior | Actual sample/source metadata and quantitative measurements |
| `hybrid` | Qualified observations and applicable quantitative predictions; calibrated fusion where available | Observations, prior or explicit reliability shrinkage | Matching sequence model and, for fusion, a matching calibrator |
| `genome` / `genome_only` | Quantitative predictions from reference or individual sequences | Applicable distance prior | Matching trained weights, reference and prior; individual work also requires genotypes, callability and ploidy |

Every mode uses fixed candidate units, distinct physical TSSs, the same activity
panel throughout a run and the same scoring kernel. RNA, H3K4me1, H3K4me3,
H3K27me3, H3K9me3, CTCF and WGBS/RRBS have [explicit interfaces](docs/MULTIOMICS.md).
They are annotations by default and may enter a separately validated classifier;
they are not arbitrary multipliers in the primary score.

## The model

```math
\boxed{
\mathrm{PACE}(E,G)=
\frac{A_\star(E)\,\overline C(E,G)\,[B(E,G)]^{\eta_{\mathrm{used}}}}
{\displaystyle\sum_{e\in\mathcal E^{\mathrm{score}}(G)}
 A_\star(e)\,\overline C(e,G)\,[B(e,G)]^{\eta_{\mathrm{used}}}}
}
```

The numerator is the current unit's activity × multi-TSS contact × optional
allocation factor. The denominator sums the same support across this gene's
actually scoreable planned candidates. Activity is the equal geometric mean of
the selected quantitative ATAC/DNase/H3K27ac signals. TSS contact uses fixed promoter
weights and compatible measurement scales. Optional allocation describes how an
unit's contact is distributed across its fixed candidate genes.

Automatic allocation uses a numerical exponent of zero without sufficient
applicable functional labels and successful grouped validation. An accepted
estimate lies in [0,1]. Details are in [calibration](docs/eta_calibration.md).

The [equation manual](docs/```math
\boxed{
\mathrm{PACE}(E,G)=
\frac{A_\star(E)\,\overline C(E,G)\,[B(E,G)]^{\eta_{\mathrm{used}}}}
{\displaystyle\sum_{e\in\mathcal E^{\mathrm{score}}(G)}
 A_\star(e)\,\overline C(e,G)\,[B(e,G)]^{\eta_{\mathrm{used}}}}
}
```.md) provides **the fully expanded total formula,
fully expanded formulas for all three modes, every symbol, and a numerical example**.
A score is a relative support share, not a causal probability or expression effect.
Missing evidence stays NA; observed zeros stay zero. No arbitrary denominator
pseudocount or gene-expression multiplier is added.

## Inspect the output before interpreting a score

| File | What to inspect |
|---|---|
| `scores.tsv.gz` | Every candidate, `pace_score`, activity, contact, support/log support, missing reasons and normalization status |
| `gene_summary.tsv` | Planned/scoreable candidates, coverage and actual denominator |
| `resolved_activity.tsv`, `resolved_contacts.tsv` | Source selection, prediction/fusion/prior, scale and status |
| `multiomics_features.tsv.gz` | RNA, auxiliary marks, methylation and other named annotations |
| `evidence.tsv`, `sources.tsv` | Traceable source and derived evidence identities |
| `eta_calibration.json` | Actual parameter, fitting/validation decision or fallback reason |
| `qc_report.json`, `report.md` | Coverage, limitations and unavailable evidence |
| `run_manifest.json`, `resolved_config.yaml` | Input/model hashes, environment, configuration and scientific contracts |
| `ml_feature_contract.json` | Measurement and feature semantics required by a compatible classifier |

Both complete and partial score sets can sum to one. `partial` means conditional
on the remaining scoreable candidates, not repaired missing data. Use
[`PACE compare`](docs/comparison.md) to recompute a common background when comparing
runs. Changes in relative support alone do not establish activity or expression
changes. There is no universal PACE significance threshold.

## Why the design is useful for livestock

The framework declares species, assembly and tissue; handles incomplete assay
panels without treating missing data as zeros; averages distinct promoters;
separates technical replicates from animals; and checks contact resolution and
normalization before combining data. Genomic reconstruction requires explicit
chromosome ploidy, callability and actual carried alleles. Frozen model scope and
grouped functional validation limit inappropriate reuse and repeated-locus leakage.
These are practical responses to heterogeneous livestock data, not evidence of
universal cross-species performance.

**The repository includes working software and synthetic examples, but no
independently validated real livestock weights or general biological accuracy
claim.** Measured data analysis can be used without a sequence model. A new
individual's WGS is informative for this workflow only when an applicable model
and candidate atlas already exist. A species with no functional training data
requires separately validated transfer. See [limitations](docs/limitations.md).

## Documentation and reproducibility

- [All parameters and defaults](docs/parameters.md), [CLI](docs/cli.md), [table dictionary](docs/data_dictionary.md).
- [Input preparation](docs/input_preparation.md), [multiomics interfaces](docs/MULTIOMICS.md), [training](docs/training.md).
- [Comparisons and benchmarks](docs/comparison.md), [worked calculations](docs/WORKED_EXAMPLES.md), [troubleshooting](docs/TROUBLESHOOTING.md).
- [Validation](docs/validation.md), [contributing](CONTRIBUTING.md), [citation](CITATION.cff), [MIT license](LICENSE).

Record the exact commit with `git rev-parse HEAD`. Save the configuration, manifest,
QC and all model identities with the analysis. File reproducible issues through
[GitHub Issues](https://github.com/shenlinyong/PACE/issues). No software-paper DOI is
claimed by this repository.
