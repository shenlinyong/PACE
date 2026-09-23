# PACE: measured activity and contact support

[中文公式与逐项解读](FORMULA.zh-CN.md) · [Tutorial](TUTORIAL.md) · [Parameters](parameters.md)

PACE ranks candidate regulatory units for a gene using experimental activity and
declared contact evidence. Every run requires measured ATAC, DNase or H3K27ac.
Missing activity is retained as unavailable. The primary score is a relative share
of support in a declared candidate set, not a causal probability or a fraction of
gene expression. The current implementation is restricted to measured activity.

## 1. Total formula

```math
\boxed{
\mathrm{PACE}(E,G)=
\frac{A_\star(E)\,\overline C(E,G)\,[B(E,G)]^{\eta_{\mathrm{used}}}}
{\displaystyle\sum_{e\in\mathcal E^{\mathrm{score}}(G)}
 A_\star(e)\,\overline C(e,G)\,[B(e,G)]^{\eta_{\mathrm{used}}}}
}
```

The numerator is the support of E for G. The denominator sums the same support
over the actually scoreable members of the frozen candidate set for G.

| Symbol | Definition | Output/configuration |
|---|---|---|
| E, G | Regulatory scoring unit and target gene | element_id, gene_id |
| e | Each candidate unit for the same gene | candidates.tsv |
| M, m | Fixed assay panel and one assay | activity.panel |
| x_obs,m | Qualified normalized experimental signal, aggregated across declared replicates | resolved_activity.tsv |
| A_star | Geometric mean of the required measured assays | A_used |
| T(G), t | Distinct physical TSSs and one TSS | promoters.tsv |
| pi(t given G) | Frozen promoter weight; within-gene weights sum to 1 | pi |
| Cbar | Promoter-weighted contact support | Cbar |
| G(E), H | Fixed candidate-gene set for E and one member | candidates.tsv |
| B | Contact fraction assigned to G within G(E) | B |
| eta_used | Actual allocation exponent in [0,1] | eta_calibration.json |
| E_score(G) | Candidates with available support in this run | scoreable, normalization_universe_id |

## 2. Measured activity

```math
A_\star(E)=\left[\prod_{m\in\mathcal M}x_{\mathrm{obs},m}(E)\right]^{1/|\mathcal M|}.
```

Supported fixed panels are ATAC, DNase, H3K27ac, ATAC+H3K27ac and DNase+H3K27ac.
For the two-assay panel, ATAC=4 and H3K27ac=9 give activity 6. If either required
assay is unavailable, activity is NA; a measured zero gives zero. The panel cannot
change per element to accommodate missing data.

Normalize and quality-control each assay upstream. The software first averages
technical repeats within biological repeats, then biological repeats within each
donor, then donors with equal weight. Individual runs require one donor; use
population_mean for a declared multi-donor summary. Aggregation takes place for
each assay before calculating its geometric mean. It does not remove batch effects.

## 3. Contact and distinct promoters

The default contact mode uses measurements. Optional contact policies use an
explicit distance prior, or a declared observed/prior shrinkage weight r:

```math
\widetilde C(E,t)=rC_{\mathrm{obs}}(E,t)+(1-r)C_{\mathrm{prior}}(E,t),
\qquad
\overline C(E,G)=\sum_{t\in\mathcal T(G)}\pi(t\mid G)\widetilde C(E,t),
\qquad \sum_t\pi(t\mid G)=1.
```

Only positive-weight sources are required: r=1 selects observations, r=0 selects
the prior, and an interior weight requires both. The current reliability setting
is declared for the run, with its calibration/source recorded; it is not an
automatically inferred posterior reliability for every edge.

```math
C_{\mathrm{prior}}(E,t)=a\left[\frac{\max(d(E,t),d_{\min})}{d_{\mathrm{ref}}}\right]^{-\gamma}.
```

Prior parameters must come from a supplied, applicable fit. The software does not
supply universal animal constants. Activity remains measured under every contact
policy. Label distance-prior results explicitly; they do not establish observed
chromatin interactions. Input contacts and priors must share resolution, scale,
normalization, balancing and measurement-window definitions.

Deduplicate transcripts sharing a physical TSS. Promoter weights are supplied or
explicitly set equal. Gene TPM does not identify promoter usage. Missing contact
at a positive-weight TSS makes Cbar unavailable; remaining weights are not renormalized.
Same-bin/near-diagonal contacts use the declared prior-or-NA policy for raw and
imported evidence. A diagonal value alone is not treated as a reliable regulatory loop.

## 4. Optional target allocation

```math
B(E,G)=\frac{\overline C(E,G)}{\displaystyle\sum_{H\in\mathcal G(E)}\overline C(E,H)}.
```

B compares candidate genes of one element; the PACE denominator compares candidate
elements of one gene. The default allocation setting is auto, falling back to zero
without eligible functional data. A nonzero learned exponent requires independent
grouped validation. At zero, B is skipped completely. At a positive exponent, the
fixed candidate-gene contact set must be available; missing genes cannot be dropped
to inflate B. See [eta calibration](eta_calibration.md) for label and split rules.

## 5. Fully expanded formula with measured contacts

```math
\mathrm{PACE}(E,G)=
\frac{
\left[\prod_{m\in\mathcal M}x_{\mathrm{obs},m}(E)\right]^{1/|\mathcal M|}
\left[\sum_{t\in\mathcal T(G)}\pi(t\mid G)C_{\mathrm{obs}}(E,t)\right]
\left[
\frac{\sum_{t\in\mathcal T(G)}\pi(t\mid G)C_{\mathrm{obs}}(E,t)}
{\sum_{H\in\mathcal G(E)}\sum_{u\in\mathcal T(H)}\pi(u\mid H)C_{\mathrm{obs}}(E,u)}
\right]^{\eta_{\mathrm{used}}}
}{
\displaystyle\sum_{e\in\mathcal E^{\mathrm{score}}(G)}
\left[\prod_{m\in\mathcal M}x_{\mathrm{obs},m}(e)\right]^{1/|\mathcal M|}
\left[\sum_{t\in\mathcal T(G)}\pi(t\mid G)C_{\mathrm{obs}}(e,t)\right]
\left[
\frac{\sum_{t\in\mathcal T(G)}\pi(t\mid G)C_{\mathrm{obs}}(e,t)}
{\sum_{H\in\mathcal G(e)}\sum_{u\in\mathcal T(H)}\pi(u\mid H)C_{\mathrm{obs}}(e,u)}
\right]^{\eta_{\mathrm{used}}}
}.
```

This expression expands activity, promoter integration and allocation. For a
declared contact prior/shrinkage analysis, replace each C_obs(E,t) by
r C_obs(E,t)+(1-r) a[max(d(E,t),d_min)/d_ref]^(-gamma), and make the same substitution
for every e and u in the numerator and denominator. Keep x_obs unchanged.

With the default zero exponent, the calculation reduces to:

```math
\mathrm{PACE}_{\eta=0}(E,G)=
\frac{
\left[\prod_{m\in\mathcal M}x_{\mathrm{obs},m}(E)\right]^{1/|\mathcal M|}
\sum_{t\in\mathcal T(G)}\pi(t\mid G)\widetilde C(E,t)}
{\displaystyle\sum_{e\in\mathcal E^{\mathrm{score}}(G)}
\left[\prod_{m\in\mathcal M}x_{\mathrm{obs},m}(e)\right]^{1/|\mathcal M|}
\sum_{t\in\mathcal T(G)}\pi(t\mid G)\widetilde C(e,t)}.
```

## 6. Additional omics and interpretation

RNA-seq, H3K4me1/3, H3K27me3, H3K9me3, CTCF and methylation have supported
[interfaces](MULTIOMICS.md). They are annotations by default, with no universal
inhibitory constants or expression multiplier. A separate classifier may use measured
features with functional labels and independent validation. Its score and calibrated
probability are distinct outputs and are never averaged or multiplied into PACE.

All planned candidates are retained. Partial normalization is conditional on its
reported scoreable set; a zero total support produces NA. A single available candidate
with positive support receives 1, which alone is not evidence of a validated link.
Compute comparisons from common denominators; changes in a relative share do not
establish changes in absolute activity or expression. Complete refers to the planned
catalog, not every biological enhancer. See [worked examples](WORKED_EXAMPLES.md).
