# Measured-data tutorial

## 1. Install and verify

Follow [installation](INSTALLATION.md), then run:

```bash
pace demo --out results/demo
pace validate --config examples/measured/config.yaml
pace run --config examples/measured/config.yaml --out results/measured
```

These commands use synthetic data and do not need external downloads beyond installation.

## 2. Prepare an experiment

1. Use one species, assembly and tissue/experimental context. Define animal identities and replicates in samples.tsv.
2. Build fixed units, physical TSSs and candidate links from processed BED/GTF annotations.
3. Quantify normalized experimental ATAC/DNase/H3K27ac signals over the same units. Freeze the activity panel for the entire run.
4. Extract compatible contacts for each element–TSS pair. If using a distance prior, supply an applicable fitted asset and label that choice explicitly.
5. Add optional RNA, CpG counts and named epigenomic features with real provenance.

See [input preparation](input_preparation.md) for executable adapters and [data dictionary](data_dictionary.md) for exact fields.

## 3. Configure a real run

Save this as `experiment.yaml` in the project root and replace identifiers and paths with your actual data:

```yaml
regime: measured
execution_profile: research
run_id: animal1_liver
target_level: individual
context:
  species: chicken
  assembly: your_assembly_identifier
  context_id: liver
inputs:
  units: prepared/units.tsv
  promoters: prepared/promoters.tsv
  candidates: prepared/candidates.tsv
  samples: prepared/samples.tsv
  sources: prepared/sources.tsv
  observed_activity: prepared/observed_activity.tsv
  observed_contacts: prepared/observed_contacts.tsv
activity:
  panel: [ATAC, H3K27ac]
contact:
  mode: observed
  scale: depth_normalized_contact
allocation:
  eta: auto
```

The default catalog is a non-overlapping 500 bp grid including eligible promoter units. Use the matching prepared catalog settings. If the entire experiment has only H3K27ac, declare `panel: [H3K27ac]`. A missing required assay within a fixed panel remains NA.

Multiple animals cannot be combined under `individual`. Use separate runs for animal comparisons or predeclare `population_mean` for equal-donor summaries. Perform normalization and batch QC upstream.

## 4. Run and inspect

```bash
pace validate --config experiment.yaml
pace run --config experiment.yaml --out results/animal1_liver
```

Inspect scores.tsv.gz, gene_summary.tsv, qc_report.json and resolved evidence. Check candidate coverage and every partial/NA reason before ranking links. The primary scores sum to 1 only for a complete planned background with positive support. Incomplete genes have separate conditional scores and sensitivity bounds. A score of 1 is not independent functional confirmation.

## 5. Optional contact prior

Fit an appropriate prior from measured contacts using `fit-contact-prior`; see [training](training.md). In the run configuration set `contact.mode: prior_only` and `contact.prior_path` to the fitted asset. Experimental activity is still mandatory. A shrinkage analysis additionally requires a justified reliability value and source; the software does not guess them.

## 6. Validate scientific claims

Independent functional labels are needed for optional allocation fitting or supervised classification. Final evaluation data must remain outside parameter selection. Compare methods on the same candidates and report coverage. See [benchmarking](comparison.md), [formula](FORMULA.md) and [limitations](limitations.md).
