# PACE Parameters: Biological Rationale, Statistical Basis, and Reproducibility

This document gives the rationale and provenance for every default parameter
in PACE. PACE is a **flexible, configurable framework**: the values below are
**biologically motivated, data-informed defaults**, not hard-coded constants.
Every value can be changed in `config/config.yaml`, and the framework provides
tools to re-derive or stress-test the defaults on your own data
(`scripts/sensitivity_analysis.py`, the optional ML module), documenting the
biological rationale, statistical basis, and reproducibility of the empirical
settings.

All stochastic steps are controlled by a single `random_seed` (default 42),
so any PACE run is exactly reproducible.

---

## 1. Signal weights

PACE aggregates epigenomic signals into an enhancer activity score. The
default weights (used only by the `weighted_geometric` and `weighted_sum`
methods; the ABC-compatible `geometric_mean` ignores weights) are:

| Signal | Type | Default weight | Biological rationale | Data-informed basis |
|--------|------|---------------:|----------------------|---------------------|
| ATAC / DNase | accessibility | **1.5** | Chromatin accessibility is the most general and information-rich predictor of regulatory potential; it is the one signal required by PACE and by the ABC model. | Highest single-feature importance in the GM12878 benchmark and in the ML feature-importance analysis (see §4); weighting it above 1.0 reflects this. |
| H3K27ac | active enhancer/promoter | **1.0** | Marks active enhancers and promoters; the second canonical ABC input. | Reference point (weight = 1.0); strong but secondary to accessibility. |
| H3K4me1 | enhancer mark | **0.8** | Marks poised and active enhancers; less specific to *active* state than H3K27ac. | Down-weighted relative to H3K27ac because it does not distinguish active from poised. |
| H3K4me3 | promoter mark | **0.5** | Promoter-proximal; informative for promoter activity but not distal enhancers. | Lower weight reflects narrower relevance to distal E-P links. |
| H3K9ac | active chromatin | **0.6** | Broadly active chromatin, partially redundant with H3K27ac. | Intermediate weight. |
| H3K36me3 | transcribed | **0.3** | Gene-body transcription mark; weak/indirect enhancer evidence. | Low weight. |
| TF ChIP-seq | TF binding | **0.3** | Site-specific regulatory protein binding; supportive but sparse and TF-dependent. | Low default; TFs can be **activators or repressors** (see §3). |

**Inhibitory weights** (α, β, γ in the inhibitory score `I(E)`):

| Signal | Default weight | Biological rationale |
|--------|---------------:|----------------------|
| DNA methylation (α) | **0.5** | Promoter/enhancer CpG methylation is associated with transcriptional silencing. |
| H3K27me3 (β) | **0.5** | Polycomb-mediated facultative repression. |
| H3K9me3 (γ) | **0.5** | Constitutive heterochromatin. |

The inhibitory weight is the **maximum fractional reduction** of activity a
fully repressed locus receives: `A_final = A × (1 − clip(Σ wₖ·Îₖ, 0, 1))`,
where each `Îₖ` is the signal scaled to [0, 1]. A weight of 0.5 means a
maximally methylated/repressed enhancer is down-weighted by up to 50%. We use
a deliberately conservative, symmetric default (0.5 for all three) because the
three repressive layers are partially redundant and we do not want any single
inhibitory layer to zero out activity on its own.

### How the defaults were chosen and how to change them

1. **Biological ordering** (accessibility ≥ active-enhancer ≥ enhancer-mark ≥
   promoter-mark ≥ TF) reflects the specificity of each mark for active distal
   regulation and is consistent with the ABC model and the broader literature.
2. **Magnitudes** were calibrated on the human GM12878 CRISPRi benchmark, the
   one dataset with experimentally validated E-P pairs, and then held fixed
   across species (we do **not** re-tune per species, to avoid over-fitting and
   to keep predictions reproducible).
3. **Robustness, not point estimates.** Because the exact magnitudes are
   uncertain, PACE is designed so that predictions are *insensitive* to them.
   `scripts/sensitivity_analysis.py` perturbs every weight (activating **and**
   inhibitory) by ±50% and reports the Jaccard overlap of the top-ranked
   predictions and AUROC stability (see §5).
4. **Data-driven alternative.** When species-/context-matched validated E-P
   data are available, the optional ML module learns the feature combination
   directly and reports feature importance, which can be compared with these
   manual weights (`python scripts/pace_ml.py importance`; see §4).

---

## 2. Aggregation method

`weighted_geometric` is the PACE default. The geometric mean is multiplicative
(an enhancer needs *several* signals to score highly, matching the biology of
combinatorial regulation) and is the form used by ABC. PACE normalizes the
weighted geometric mean by the **sum of weights**:

```
A(E) = ( ∏ᵢ Sᵢ^wᵢ )^(1 / Σᵢ wᵢ)
```

The `1/Σwᵢ` normalization (a genuine weighted geometric **mean**, not an
unnormalized product) makes the activity scale independent of the absolute
weight magnitudes, so doubling all weights does not change the ranking. This
is implemented identically in `scripts/multiomics_activity.py` and
`workflow/scripts/neighborhoods.py`.

`geometric_mean` (unweighted) reproduces the original ABC activity and is
provided for ABC-compatibility; `weighted_sum` and `arithmetic_mean` are
available for experimentation.

---

## 3. Transcription factors are not assumed to be activators

A fixed positive TF weight would incorrectly assume every TF enhances activity.
PACE therefore lets each TF be declared an **activator** or a **repressor**:

- per sample, via the `TF_modes` column of the biosamples table
  (`activator` / `repressor`, comma-separated, one per TF); or
- globally, via `transcription_factors.specific_tfs` in `config.yaml`
  (`{weight, inhibitory}` per TF).

Repressor TFs contribute to the inhibitory score `I(E)` exactly like
repressive histone marks. The default, when no mode is given, is *activator*,
but this is an explicit, documented choice rather than a hidden assumption.

---

## 4. Empirical weights vs. learned importance

To verify that the manual weights are biologically sensible, train the ML
module on validated pairs and compare:

```bash
python scripts/pace_ml.py train     --predictions preds.tsv.gz --validation pairs.tsv --output model.pkl --balance_classes
python scripts/pace_ml.py importance --model model.pkl --output importance/run
```

This writes:
- `*.gain_importance.tsv` — impurity/gain-based importance;
- `*.permutation_importance.tsv` — permutation importance (robust, unbiased);
- `*.empirical_vs_learned.tsv` — manual weight vs learned importance, both
  normalized, with ranks;
- `*.png` — bar plot.

In our benchmark the **rank ordering** of the learned importances matches the
manual ordering (accessibility > H3K27ac > other marks), supporting the
biological validity of the defaults while making any discrepancy transparent.
The `Combined.Score` (formula × ML) is reported alongside the standalone ML
score so users can confirm whether combining helps for their data.

---

## 5. Sensitivity / robustness protocol

```bash
python scripts/sensitivity_analysis.py \
    --predictions preds_with_signals.tsv --output sensitivity.tsv \
    --validated_col validated --plot sensitivity.png
```

For each weight (activating and inhibitory), PACE rescales it across ±50% and
reports, relative to the default configuration:
- **Jaccard similarity** and **overlap %** of the top-10% predictions;
- **Spearman correlation** of the full score vector;
- **AUROC** (if a validation column is provided).

Stable top-prediction sets (high Jaccard) demonstrate that the *biological
conclusions*, not just the global AUROC, are robust to the weight choices.

---

## 6. Contact model

| Parameter | Default | Basis |
|-----------|--------:|-------|
| Power-law γ (exponent) | 1.024 | Fit to the average Hi-C contact-vs-distance decay. |
| Power-law scale | 5.959 | Companion scale for the fitted decay. |
| Pseudocount | 5,000 bp | Avoids division by zero at very short range; one Hi-C bin. |
| Max E-G distance | 5,000,000 bp | Standard ABC search window covering essentially all known long-range E-P links. |
| Hi-C resolution | 5,000 bp | Resolution at which contact matrices are read/normalized (KR). |

When Hi-C is available it is used directly (KR-normalized); otherwise the
power-law model estimates contact from distance. The two paths are implemented
once in `workflow/scripts/hic.py` and shared by both the Snakemake pipeline and
the standalone calculator.

---

## 7. Other thresholds

| Parameter | Default | Basis |
|-----------|--------:|-------|
| Score threshold | 0.02 | Standard ABC operating point; raise to 0.03–0.05 for precision, lower to 0.015 for recall. |
| Min expression | 1 TPM | Conventional minimum for "expressed"; only expressed genes can be regulated. |
| MACS2 p-value | 0.1 | Permissive peak calling for candidate-region discovery (filtered later by activity × contact). |
| Peak extension | 250 bp (→500 bp regions) | Standard ABC candidate-element size. |
| N strongest peaks | 150,000 | Caps candidate count for tractability while retaining all confident peaks. |
| eQTL p-value | 1e-5 / genome-wide FDR < 0.05 | Significance threshold for validation variants. |
| ML n_estimators / max_depth / lr | 100 / 5 / 0.1 | Standard, conservative gradient-boosting settings that limit over-fitting on the small validated sets typical of E-P prediction. |

---

## 8. Reproducibility checklist

- `random_seed` (default 42) controls ML class balancing, cross-validation
  shuffles, and size-matched control sampling.
- All signal quantification is deterministic (no simulated/random values).
- Power-law parameters, weights, and thresholds are all in `config.yaml`.
- `python test_pace_complete.py` runs 45 tests covering the core pipeline and
  every analysis module.
- Software versions are pinned in `requirements.txt` /
  `workflow/envs/pace-env.yml`.
