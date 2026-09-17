# PACE canonical software interface

PACE converts measured or predicted activity and contact evidence into relative enhancer–gene
support. The installable `pace-livestock` package implements the September 2026 contracts.
The older `scripts/pace.py` interface remains available for reproducing region-based analyses;
its missing-data rules and defaults differ. Do not mix their scores.

## Install and run

Python 3.11 or newer is required. From this repository:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install .
pace-livestock demo --regime measured --out results/demo_measured
pace-livestock demo --regime hybrid --out results/demo_hybrid
pace-livestock demo --regime genome_only --out results/demo_genome
```

The base package requires NumPy and PyYAML, with no GPU or runtime downloads. All three demos
use synthetic inputs and deterministic fixture weights. Training/inference with the reference
CNN requires `pip install '.[sequence]'`. Genomic binary adapters use `.[io]`.
`.[ml]` installs the optional scientific ML ecosystem; the transparent elastic-net solver
itself uses NumPy. Development and regression testing use `.[dev,io,sequence,ml]`.
Installation may access package indexes; running the bundled examples does not.

To use supplied tables:

```bash
pace-livestock validate --config examples/measured/config.yaml
pace-livestock capabilities --config examples/genome_only/config.yaml
pace-livestock run --config examples/measured/config.yaml --out results/measured
```

Paths resolve relative to the YAML file. Unknown and duplicate keys are errors. Commands refuse
existing output directories, return a nonzero exit status on invalid input, and publish outputs
only after successful computation. A partial/empty score set is a legitimate scientific result;
inspect the QC report rather than interpreting successful execution as complete evidence.

## Choose an evidence regime

| Regime | Minimum activity input | Minimum contact input | Real assets |
|---|---|---|---|
| `measured` | Fixed ATAC, DNase, H3K27ac or accessibility+H3K27ac panel | Observed contact or explicit prior | No sequence model required |
| `hybrid` | Qualified observations and/or matching quantitative predictions | Observations/prior | Applicable sequence model when predictions are used; calibrator when sources are fused |
| `genome_only` | Reference/individual windows or manifest-matched precomputed quantitative predictions | Applicable distance prior | Species/assembly/context/target-matched quantitative weights |

Each run also declares `bulk_proxy`, `individual` or `population_mean` target level, species,
assembly and tissue/context. For an individual hybrid run, observed donor IDs must match the
genotype individual ID. Prediction rows have model/evidence IDs and no invented sample IDs.

`demonstration` accepts only synthetic assets. `research` rejects synthetic assets but allows
real weights whose particular task has not been independently validated. `validated` additionally
requires task-specific reports, matching context and report checksums. It is intentionally not
a switch that turns an unvalidated model into a validated one. No real livestock weights or
biological validation reports are distributed with this package.

## Commands

Every operation takes `--config FILE --out NEW_DIRECTORY`, except `validate` and `capabilities`,
which print a report. All registered commands execute real code paths.

| Command | Purpose | Configuration contract |
|---|---|---|
| `prepare` | BED/GTF catalog, bigWig, cool/mcool, CpG or transcript RNA adapter | [Input preparation](input_preparation.md) |
| `run` / `validate` | Shared three-mode scoring pipeline and full input checks | [Parameters](parameters.md) |
| `capabilities` | Asset availability, blockers and task status | Same as run |
| `prepare-genome` | Fixed-target haplotype FASTA and mapping report | Same as run, with sequence manifest |
| `predict-sequence` | Individual windows and quantitative predictions | Same as run |
| `train-sequence` | Actual quantitative CNN training | [Model training](training.md) |
| `fit-contact-prior` | Zero-inclusive distance-bin fit | [Model training](training.md) |
| `fit-fusion` | Independent-target convex log-space calibration | [Model training](training.md) |
| `compare` | Full and conditional support-share changes | [Comparisons](comparison.md) |
| `stability` | Common-denominator replicate agreement | [Comparisons](comparison.md) |
| `variant-effects` | Separate single-variant REF/ALT signal scenarios | [Comparisons](comparison.md) |
| `train` / `predict-ml` | Grouped elastic-net and frozen inference | [Model training](training.md) |
| `benchmark` | Negative distance, single-TSS ABC-style, PACE eta 0/1 and external scores | [Comparisons](comparison.md) |

## Interpretation

The main score is `A × Cbar × B^eta` normalized within the actual scoreable candidates for each
gene. Default `eta=0`; eta 1 is a separate comparison. The assay panel is fixed for the run.
True zeros remain zero, technical absence remains NA, and required positive-pi TSS contacts
are never dropped to renormalize promoter weights. No arbitrary epsilon or residual support
is added. Log-space support preserves valid shares when raw products overflow.

A score is a relative support share, not a causal probability, expression fold change or
percentage of expression. Both complete and partial denominators may sum to one. A high score
with one remaining candidate therefore says little about biological certainty.

For comparisons, inspect delta activity, delta support, both gene totals and delta PACE together.
The software recomputes shares on a common measurable subset and keeps that conditional change
separate from a complete-background change. Sequence predictions are averaged by assay across
copies **before** activity is constructed; summing copy-specific supports is a different,
unimplemented research estimand.

See the [model contract](model.md), [data dictionary](data_dictionary.md),
[validation record](validation.md), [limitations](limitations.md), and [Chinese quick start](../README.zh-CN.md).
