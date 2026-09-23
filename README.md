# PACE

PACE estimates relative enhancer–gene regulatory support from **measured activity**
and declared promoter-contact evidence in livestock research.

[中文](README.zh-CN.md) · [Full Chinese manual](docs/USER_GUIDE.zh-CN.md) ·
[Formula](docs/FORMULA.md) · [Inputs](docs/data_dictionary.md) · [Parameters](docs/parameters.md)

## Download and install

```bash
git clone https://github.com/shenlinyong/PACE.git
cd PACE
conda env create -f environment.yml
conda activate pace
pace --version
```

Alternatively, use Python >=3.11:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install '.[io,ml]'
```

For Docker:

```bash
docker build -t pace:local .
docker run --rm --user "$(id -u):$(id -g)" -v "$PWD":/work pace:local \
  demo --out /work/results/docker_demo
```

See [installation](docs/INSTALLATION.md) for a wget source archive, paths and upgrades.
No GPU or deep-learning framework is required.

## First run

```bash
pace demo --out results/demo
pace validate --config examples/measured/config.yaml
pace run --config examples/measured/config.yaml --out results/measured
pace measured --help
```

The bundled data are synthetic software examples. Use experimental data with
`execution_profile: research` for scientific analysis. Output directories must be new
unless `--force` is used to replace a recognized PACE output with a retained backup.

## Required evidence

| Component | Input and role |
|---|---|
| Activity | Measured ATAC, DNase or H3K27ac; one fixed supported single/two-assay panel |
| Contact | Compatible Hi-C or other contact assays; an explicitly supplied distance prior is optional |
| Candidate catalog | Fixed units, distinct TSSs and candidate element–gene pairs |
| Samples | Donor, biological/technical replicate, species, assembly and tissue metadata |
| Additional omics | RNA, histone marks, CTCF and methylation annotations; optional separate classifier |

Use [input preparation](docs/input_preparation.md) for BED/GTF, bigWig and cool/mcool.
PACE starts from processed experimental data; it does not align raw reads or call peaks.
Missing activity remains unavailable. `pace run` and `pace measured` call the same
measured-activity implementation. [Migration](docs/migration.md) describes removed APIs.

## Total formula

```math
\boxed{
\mathrm{PACE}(E,G)=
\frac{A_\star(E)\,\overline C(E,G)\,[B(E,G)]^{\eta_{\mathrm{used}}}}
{\displaystyle\sum_{e\in\mathcal E(G)}
 A_\star(e)\,\overline C(e,G)\,[B(e,G)]^{\eta_{\mathrm{used}}}}
}
```

A_star is the fixed-panel geometric mean of measured signals. Cbar integrates
distinct promoters with fixed weights. B is optional target allocation; eta defaults
to zero without eligible independent functional validation. The primary denominator
includes all planned candidates for the gene. If support is incomplete, the primary
score is withheld by default; a separate conditional score uses only scoreable
candidates, and sensitivity bounds describe the unresolved support.
See [fully expanded equations](docs/FORMULA.md) and [worked examples](docs/WORKED_EXAMPLES.md).

## Analyze your experiment

```bash
pace run --config experiment.yaml --out results/experiment
```

A complete configuration and all required table fields are explained in the
[tutorial](docs/TUTORIAL.md). YAML paths are relative to the YAML file, whereas
direct CLI file paths are relative to the current working directory. For multiple
animals use separate individual runs or explicitly declare `target_level: population_mean`.
Replicates require correct sample metadata and comparable normalization.

## Read the results

| Output | Meaning |
|---|---|
| scores.tsv.gz | Relative support, activity, contact, reasons and coverage for every candidate |
| gene_summary.tsv | Actual denominator and complete/partial/zero-support status |
| resolved_activity.tsv / resolved_contacts.tsv | Measurements and the contact evidence actually used |
| multiomics_features.tsv.gz | Additional measured annotations |
| eta_calibration.json | Actual exponent, validation result and fallback reason |
| qc_report.json / run_manifest.json | QC, context, input hashes and reproducibility information |

PACE scores are relative support shares. They are not causal probabilities or
expression effects. Additional classifier outputs are separate. Software tests do
not demonstrate livestock predictive accuracy. See [limitations](docs/limitations.md),
[verification](docs/validation.md), [benchmarking](docs/comparison.md) and
[manuscript scope](docs/MANUSCRIPT_SCOPE.zh-CN.md).

Maintainer: Linyong Shen (申林用). [MIT license](LICENSE).

## Sparse experimental data and preparation

The default primary score requires complete planned support. Conditional scores and support-assumption intervals remain separate. Matched contact priors enable sparse-count regularization, and diagonal contacts use declared priors or recorded neighbor correction. Region summaries aggregate unique grid cells.

Use `pace init`, `prepare-pairs`, `merge-tables`, `normalize-activity`, `prepare-promoter-weights`, and `fit-prior` to prepare a project without hand-written joins. The explicit human contact-shape baseline is unvalidated for livestock and remains measured-activity-only.

[Practical workflow](docs/PRACTICAL_WORKFLOW.md) · [Robustness study protocol](docs/ROBUSTNESS_VALIDATION.md) · [Development](docs/DEVELOPMENT.md).
