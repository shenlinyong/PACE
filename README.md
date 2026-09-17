# PACE

[![Software tests](https://github.com/shenlinyong/PACE/actions/workflows/ci.yml/badge.svg)](https://github.com/shenlinyong/PACE/actions/workflows/ci.yml)

**Activity–contact support for livestock enhancer–gene research.**

PACE combines quantitative activity, promoter contact and optional target allocation
into a gene-normalized relative support score. The same model powers **measured**,
**hybrid** and **genome-only** analysis. Developed by **申林用 (Linyong Shen,
shenlinyong), Northwest A&F University (西北农林科技大学)**.

[中文说明](README.zh-CN.md) · [Formula](docs/FORMULA.md) · [Command line](docs/cli.md) ·
[Input preparation](docs/input_preparation.md) · [Training](docs/training.md) ·
[Validation](docs/validation.md)

## Current formula

```math
\boxed{
\mathrm{PACE}(E,G)=
\frac{A_\star(E)\,\overline C(E,G)\,[B(E,G)]^{\eta_{\mathrm{used}}}}
{\displaystyle\sum_{e\in\mathcal E^{\mathrm{score}}(G)}
 A_\star(e)\,\overline C(e,G)\,[B(e,G)]^{\eta_{\mathrm{used}}}}
}
```

The denominator contains only the support of the actually scoreable candidates of
that gene. The fixed activity panel uses an equal geometric mean; contact is
reliability-adjusted and averaged over distinct TSSs using declared promoter weights.
`B` is the contact share over an element's fixed candidate-gene set.

Default `--eta auto` uses **eta=0** when applicable functional calibration is absent
or has insufficient information. `--eta-labels` estimates a continuous value in
**[0,1]** from eligible training/calibration perturbations; test labels never fit eta.
`--eta-model` reuses a frozen calibration. A predeclared `--eta 0.35` is also valid.
At eta=0, B is skipped. True zeros stay zero; required missing evidence stays NA.
A zero gene-support total gives NA. [Definitions and boundary cases](docs/FORMULA.md).

## Install and run

Python 3.11 or newer is required:

```bash
git clone https://github.com/shenlinyong/PACE.git
cd PACE
python -m venv .venv
source .venv/bin/activate
python -m pip install .
PACE --version
PACE run --help
```

An isolated prefix installation is available through `./install.sh`; see
[installation](docs/INSTALLATION.md). Install `'.[io]'` for bigWig/cool/BCF adapters
or `'.[sequence]'` for real CNN training/inference. The command aliases `PACE`,
`pace` and `pace-livestock` all use the same implementation.

These commands run the bundled **synthetic** fixtures:

```bash
PACE measured --config examples/measured/config.yaml --out results/measured
PACE hybrid --config examples/hybrid/config.yaml --out results/hybrid
PACE genome --config examples/genome_only/config.yaml --out results/genome
```

YAML is optional: `PACE --mode measured|hybrid|genome` accepts direct file options,
including `--catalog-dir`, `--activity`, `--contacts`, `--reference`,
`--sequence-model`, `--contact-prior` and genotype options.
The [command guide](docs/cli.md) contains complete real-input command templates.
CLI paths resolve from the working directory; YAML paths resolve from the YAML file.
Existing output directories are never overwritten.

## Choose a mode

| Mode | Activity evidence | Contact evidence | Required assets |
|---|---|---|---|
| `measured` | Qualified measurements for a fixed assay panel | Qualified observations or an explicit compatible prior | No sequence model required |
| `hybrid` | Measurements and applicable quantitative sequence predictions | Observations and/or prior | Sequence weights for prediction; matching calibrator for fusion |
| `genome` / `genome_only` | Quantitative reference or individual sequence predictions | Compatible contact prior | Species, assembly, tissue and target-matched weights and prior |

All modes use the same canonical units, bulk marginal aggregation and numerical
kernel. RNA, auxiliary marks and methylation are annotations by default; an
independent classifier may use named features when suitable functional labels exist.
They are not extra multipliers in the primary formula.

## Read the results

| Output | Contents |
|---|---|
| `scores.tsv.gz` | All candidate edges, activity/contact/allocation, support, score, missingness and normalization status |
| `gene_summary.tsv` | Actual gene denominators and candidate coverage |
| `resolved_activity.tsv`, `resolved_contacts.tsv` | Measurements, predictions, fusion and prior provenance |
| `eta_calibration.json` | Actual eta, fit/fallback status, eligible labels and reusable scope |
| `qc_report.json`, `run_manifest.json`, `resolved_config.yaml` | QC, input/model/software hashes, environment and configuration |

PACE is a **relative support share**, not a causal probability or expression effect.
Partial normalization may sum to one while omitting unscoreable candidates. Use
`PACE compare` to recompute a common denominator for compatible sample comparisons.
There is no universal link threshold or built-in biologically optimal livestock prior.

The repository supplies functioning training/inference code and synthetic tests.
**Validated real livestock weights and independent biological performance evidence
are not bundled.** See [limitations](docs/limitations.md).

## Documentation and citation

[Model details](docs/model.md) · [Parameters](docs/parameters.md) ·
[Data dictionary](docs/data_dictionary.md) · [Eta calibration](docs/eta_calibration.md) ·
[Comparisons](docs/comparison.md) · [Worked calculations](docs/WORKED_EXAMPLES.md) ·
[Troubleshooting](docs/TROUBLESHOOTING.md) · [Contributing](CONTRIBUTING.md)

Record the exact software commit in research outputs. Authorship is recorded in
[CITATION.cff](CITATION.cff); no software-paper DOI is claimed. Superseded interfaces
are available through Git history, with their [migration context](docs/migration.md).
