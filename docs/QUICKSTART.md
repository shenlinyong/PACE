# PACE quick start

Follow the installation and quantified-table example in the [repository README](../README.md#quick-start). Use the table interface when assay signals, scales, contact expectations and QC fields have already been prepared.

```bash
python scripts/prepare_tss.py --gtf annotation.gtf.gz --output results/tss.tsv
python scripts/pace.py --pairs pairs.tsv --activity-config activity.json --output results/scores.tsv
```

Prepare one candidate-enhancer/TSS row per sample and assembly. Coordinates are zero-based and half-open. Supply missing measurements as NA, rather than zero. The [input contract](INPUTS.md) lists required and optional fields.

For a fully specified synthetic run with bedtools available:

```bash
python scripts/smoke_test.py --output-dir results/smoke
```
