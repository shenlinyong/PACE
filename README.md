# PACE

Prediction of Activity-based regulatory Connections for Enhancers.

**Author and maintainer: shenlinyong — 申林用 (Linyong Shen), Northwest A&F University.**

PACE ranks candidate enhancer–gene links from activity and contact evidence. It retains missing measurements, combines distinct transcription start sites, and reports input-evidence quality separately from the relative score. Scores are not calibrated probabilities.

## Model

$$
\mathit{PACE}(E,G)=\frac{A(E)C(E,G)B(E,G)^\eta}{\sum_{e\in\mathcal E^{\mathrm{obs}}(G)}A(e)C(e,G)B(e,G)^\eta+U(G)}.
$$

Here, $A(E)$ is enhancer activity, $C(E,G)$ is gene-level contact, $B(E,G)$ is enhancer-centric target allocation, and $U(G)$ is independently estimated residual support. The default uses $\eta=1$ and computationally sets unknown $U(G)$ to zero, with that uncertainty retained in the output. RNA can annotate gene eligibility; it does not multiply the current primary score.

For comparison, the earlier expression-weighted model is defined as

$$
\mathit{PACE}(E,G)=W_{\mathrm{expr}}(G)\frac{A(E)C(E,G)}{\sum_{e\in\mathcal E(G)}A(e)C(e,G)}.
$$

The expression weight is applied **after** within-gene normalization and does not enter its denominator. This correction also applies to the retained legacy documentation. Historical model scores and current model scores use different definitions and should not be pooled.

## Quick start

Python 3.10 or newer is required. Run the following from the repository root:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-core.txt
python scripts/pace.py \
  --pairs example_quantified/candidates.tsv \
  --activity-config example_quantified/activity.json \
  --output results/example_predictions.tsv
```

The bundled data are synthetic software fixtures. They produce six candidate links, retaining unscored and provisional cases, and do not establish biological performance.

## Tests and small-file workflow

```bash
python -m pip install -r requirements-test.txt
python -m pytest tests -q
python scripts/smoke_test.py --output-dir results/smoke
```

The smoke workflow additionally requires `bedtools` on `PATH`. It runs TSS preparation, read quantification, both prediction interfaces, measured BEDPE contacts, filtering, QC and the quantified-table example. It compares scores with a separate calculation. Test fixtures distinguish observed zero, missing measurements, unknown quality and omitted candidates.

For the five-step bundled read example, use `bash example/run_example_direct.sh`. Full optional bioinformatics and Snakemake dependencies are listed in `requirements.txt`; the complete FASTQ-to-prediction workflow has not been validated by the small-file checks.

## Documentation

- [Equations and model boundaries](docs/FORMULA.md)
- [Methods](docs/METHODS.md)
- [Canonical notation and code fields](docs/NOTATION.md)
- [Inputs, outputs and migration](docs/INPUTS.md)
- [Quick start](docs/QUICKSTART.md) and [workflow tutorial](docs/TUTORIAL.md)
- [Validation scope](VALIDATION.md) and [changes](CHANGELOG.md)
- [Historical methods](legacy/docs/METHODS.md), retained with the corrected expression-weighted formula

The single scoring kernel is `workflow/scripts/pace_core.py`. The table interface, signal-file calculator and workflow predictor delegate to it. Optional legacy ML tools are archived separately and do not contribute to primary predictions. The initial quality cutoff and descriptive score cutoff are operational settings, not calibrated FDR thresholds.

## Attribution and licence

PACE builds on the Activity-by-Contact model and incorporates the attributed target-allocation and multiple-TSS ideas described in the methods. Third-party attribution is retained. See [AUTHORS.md](AUTHORS.md) and the [MIT licence](LICENSE).
