# PACE workflow tutorial

Start with the [quick start](QUICKSTART.md) and [input contract](INPUTS.md).

## Quantified candidates

Use `scripts/pace.py` for prepared candidate/TSS tables. Save the activity configuration, annotation catalogue and complete unfiltered output alongside each run. Do not use benchmark labels to construct candidate sets or input-quality values.

## Read or signal files

`scripts/calculate_pace_score.py` accepts a biosample configuration, candidate BED and gene/TSS table. The workflow interface splits neighborhood quantification and prediction between `workflow/scripts/pace_neighborhoods.py` and `workflow/scripts/pace_predict.py`.

Run `bash example/run_example_direct.sh` for the bundled five-step read example, or `python scripts/smoke_test.py --output-dir results/smoke` for independently checked synthetic inputs. Both require bedtools. The supported activity mode is `missing_geometric`.

## Optional contact measurements

Supply compatible expected contacts, reliability and source labels through `--contact_metadata`. A readable Hi-C file alone is not evidence that normalization or tissue provenance matches. Missing metadata leads to a labelled prior-only prediction.

## Snakemake

Install the optional workflow dependencies and configure references and samples before running:

```bash
snakemake --snakefile workflow/Snakefile --configfile config/config.yaml --cores 4
```

The small-file smoke workflow does not validate a complete raw-read Snakemake run. Archived supervised tools and historical expression-weighted scoring are outside the current primary workflow.
