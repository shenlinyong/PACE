# Quantified-table example

PACE (Prediction of Activity-based regulatory Connections for Enhancers) scores the supplied synthetic candidate/TSS table using `activity.json` to specify signal columns, scales and quality.

```bash
python scripts/pace.py \
  --pairs example_quantified/candidates.tsv \
  --activity-config example_quantified/activity.json \
  --output results/quickstart/predictions.tsv
```

Run from the repository root in the [Conda environment](../docs/INSTALLATION.md). Nine input rows represent three enhancers, two genes and three distinct TSSs; aggregation produces six gene-level links. Two scores remain missing because activity is unavailable. The fixture also includes measured-zero signal, absent contact and unknown annotation quality.

See the [quick start](../docs/QUICKSTART.md) for inspection commands and the [input specification](../docs/INPUTS.md) before supplying biological data. These fixtures do not establish biological performance.
