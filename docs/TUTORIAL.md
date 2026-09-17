# From inputs to reproducible PACE scores

This tutorial uses the installed PACE command and the [current formula](FORMULA.md).
The [command guide](cli.md) includes full file-parameter templates for real projects.
The runnable commands below use synthetic fixtures and do not validate livestock biology.

## 1. Prepare a fixed catalog and qualified evidence

Create canonical nonoverlapping units, distinct physical promoters and frozen candidate
edges from compatible annotations. Use processed quantitative signals and contact data
with matching species, assembly, tissue, units and windows. The
[input-preparation guide](input_preparation.md) provides BED/GTF, bigWig, cooler,
CpG and RNA adapter configurations. It also documents genotype and callability inputs.

Real measured data do not require a sequence model. Hybrid/genome predictions require
appropriate quantitative weights; a FASTA file alone is insufficient. Fix the assay
panel and target level before scoring.

## 2. Validate and run

```bash
PACE validate --config examples/measured/config.yaml
PACE measured --config examples/measured/config.yaml --out results/tutorial_measured
PACE hybrid --config examples/hybrid/config.yaml --out results/tutorial_hybrid
PACE genome --config examples/genome_only/config.yaml --out results/tutorial_genome
```

Inspect `qc_report.json`, every edge's reason and each gene's actual denominator.
A successful process does not imply complete evidence or biological validation.

## 3. Calibrate eta only with appropriate functional evidence

```bash
PACE fit-eta --config examples/measured/config.yaml \
  --eta-labels examples/training/eta_labels.tsv --eta-min-genes 2 \
  --out results/tutorial_eta
PACE measured --config examples/measured/config.yaml \
  --eta-model results/tutorial_eta/eta_calibration.json \
  --out results/tutorial_frozen
```

This toy example deliberately contains two genes and explicitly sets that minimum.
The normal guard is three informative genes; neither threshold is a biological
sample-size justification. Real labels must satisfy the
[functional schema and split rules](eta_calibration.md). Without applicable fitting
evidence the automatic exponent remains zero with a recorded fallback reason.
Use one frozen compatible calibration for sample comparisons.

## 4. Evaluate and compare

`PACE compare --config comparison.yaml --out results/comparison` recomputes support
shares on common scoreable backgrounds. `PACE benchmark` evaluates independent
functional labels, coverage and declared controls. See [configuration examples](comparison.md).

A relative share may change because another element's support changed. It is not
an expression fold change or causal probability. There is no built-in universal link
cutoff. Train a separate classifier only for a declared task with appropriate data;
see [training](training.md). Keep all input/model/software hashes with research results.
