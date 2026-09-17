# Synthetic read example

From the repository root, run `bash example/run_example_direct.sh`. Install the core/test dependencies and bedtools first. Inputs are synthetic software fixtures, not biological validation data. The script produces current PACE predictions in `example/results/Example_Sample/`.

The separate `example_quantified/` directory supports the lightweight table CLI. `scripts/smoke_test.py` supplies the independently checked eight-command example. Optional Snakemake execution is described in the main workflow tutorial.
