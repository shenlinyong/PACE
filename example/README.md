# Genomic-read example

PACE (Prediction of Activity-based regulatory Connections for Enhancers) includes synthetic accessibility/H3K27ac reads, a small genome/annotation and candidate peaks. No external download is needed.

From the repository root after creating the [Conda environment](../docs/INSTALLATION.md):

```bash
conda activate pace
bash example/run_example_direct.sh
```

The five direct stages are candidate construction, neighborhood quantification, prediction, filtering and QC. The unfiltered output at `example/results/Example_Sample/Predictions/EnhancerPredictionsAllPutative.tsv.gz` contains 12,000 candidate enhancer–gene links. The same sample directory contains `Neighborhoods/`, filtered `Predictions/` and `Metrics/` outputs.

The [tutorial](../docs/TUTORIAL.md) gives the individual commands and a transcript-TSS route. The optional [Snakemake environment](../docs/INSTALLATION.md#optional-peak-calling-and-snakemake-environment) adds MACS2 and scheduling; it may produce a different candidate set because it calls peaks rather than using the supplied narrowPeak file.

These files test software behavior, not biological accuracy. Generated results are not tracked in Git. Do not run different examples concurrently into the same output directory.
