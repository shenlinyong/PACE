# Quick start

Install the current package using [these instructions](INSTALLATION.md). From the
repository root, run the three explicitly synthetic fixtures:

```bash
PACE measured --config examples/measured/config.yaml --out results/quick_measured
PACE hybrid --config examples/hybrid/config.yaml --out results/quick_hybrid
PACE genome --config examples/genome_only/config.yaml --out results/quick_genome
```

Or generate self-contained offline inputs and results from any directory:

```bash
PACE demo --regime measured --out demo_measured
PACE demo --regime hybrid --out demo_hybrid
PACE demo --regime genome_only --out demo_genome
```

A normal run writes `scores.tsv.gz`, `gene_summary.tsv`, resolved evidence,
`eta_calibration.json`, `qc_report.json`, `run_manifest.json` and the resolved YAML.
The default automatic eta is numerically zero in these examples because no
functional fitting labels were supplied. Existing output directories are refused.

For actual data, supply matching canonical tables and real applicable model assets.
YAML is optional: [direct input commands](cli.md) use `PACE --mode ...` plus file
and context options. `PACE run --help` lists them. The [tutorial](TUTORIAL.md)
walks through preparation, calibration, fixed reuse and interpretation.
