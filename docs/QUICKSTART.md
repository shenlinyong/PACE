# Quick start

PACE (Prediction of Activity-based regulatory Connections for Enhancers) can score the included numerical fixture without downloading a genome or biological dataset. Start in the repository root after [Conda installation](INSTALLATION.md).

## 1. Activate the environment

```bash
conda activate pace
python scripts/pace.py --help
```

## 2. Inspect the inputs

The [candidate table](../example_quantified/candidates.tsv) contains three enhancers, two genes and three distinct TSSs. The [activity JSON](../example_quantified/activity.json) names the assays, scales and quality columns. Each enhancer–TSS combination is a separate row; missing fields are intentional.

```bash
head -4 example_quantified/candidates.tsv
cat example_quantified/activity.json
```

## 3. Score all candidates

```bash
python scripts/pace.py \
  --pairs example_quantified/candidates.tsv \
  --activity-config example_quantified/activity.json \
  --output results/quickstart/predictions.tsv
```

Expected: `Wrote 6 gene-level edges; all scores are uncalibrated.`

```bash
python - <<'PYCODE'
import pandas as pd
p = pd.read_csv('results/quickstart/predictions.tsv', sep='\t')
print(p[['start', 'TargetGeneEnsemblID', 'n_tss', 'PACE.Score',
         'evidence_status']].to_string(index=False))
assert len(p) == 6
assert p['PACE.Score'].isna().sum() == 2
PYCODE
```

Nine TSS-level rows become six gene-level links. Two scores remain `NA` because enhancer activity was not observed. Scores for the remaining rows are conditional on the available candidate background. A large score does not eliminate missing evidence.

Expected values, rounded to six decimal places (`NaN` is pandas' display of a missing number):

| Enhancer start | Gene | Distinct TSSs | PACE score | Evidence status |
| --- | --- | --- | --- | --- |
| 1000 | g1 | 2 | 0.828036 | `provisional` |
| 1000 | g2 | 1 | 0.819960 | `provisional` |
| 0 | g2 | 1 | 0.180040 | `provisional` |
| 0 | g1 | 2 | 0.171964 | `provisional` |
| 2000 | g1 | 2 | `NA` | `insufficient` |
| 2000 | g2 | 1 | `NA` | `insufficient` |

## 4. Create a separate filtered output

```bash
python workflow/scripts/pace_filter.py \
  --predictions results/quickstart/predictions.tsv \
  --output results/quickstart/selected.tsv \
  --threshold 0.02
```

The original `predictions.tsv` remains the unfiltered record. The filter produces `selected.tsv` and `selected_Full.tsv`; these are both selected rows, with different column sets. The cutoff illustrates syntax and is not a species-specific statistical threshold.

## 5. Verify the installation

```bash
python -m pytest tests -q
python scripts/smoke_test.py --output-dir results/smoke
```

The small-file workflow prints a final `PASS` report and compares outputs with a separate numerical calculation. To process supplied genomic-read files next, run `bash example/run_example_direct.sh` or follow the [step-by-step tutorial](TUTORIAL.md).

## Choose your next route

To see the effect of a missing assay, missing contact or measured zero contact using runnable inputs, continue with [worked examples](WORKED_EXAMPLES.md).

| Your data | Next step |
| --- | --- |
| Already quantified enhancer signals, contacts and quality | Prepare the [input table](INPUTS.md#quantified-candidate-table) and use `scripts/pace.py` |
| Called accessibility peaks and BAM/tagAlign/bigWig signals | Follow the direct-file [tutorial](TUTORIAL.md#step-by-step-genomic-file-example) |
| Already aligned reads requiring peak calling | Use the optional [Snakemake workflow](TUTORIAL.md#optional-snakemake-workflow) |
| An error or unexpected empty/NA output | Check [troubleshooting](TROUBLESHOOTING.md) and the output reason codes |
