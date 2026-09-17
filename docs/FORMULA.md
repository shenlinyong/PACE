# PACE: model equations and interpretation

[中文公式说明：逐项解读与三种模式完整展开](FORMULA.zh-CN.md)

This document defines the formula evaluated by the software, including the expanded
formula for each evidence mode. For configuration and worked commands, see the
[Chinese user manual](USER_GUIDE.zh-CN.md), [English tutorial](TUTORIAL.md), and
[complete parameter reference](parameters.md).

## 1. What the score means

For one gene, PACE compares the support supplied by each candidate regulatory unit.
It first estimates the unit's activity, then its contact with the gene's promoters,
and optionally its allocation to that gene. It divides this support by the sum for
that gene's scoreable candidate units. The score answers **which candidates have
more relative support in this specified background**. It is not the probability
that a regulatory interaction is real or a percentage of gene expression.

A fixed scoring unit may carry an enhancer role, a promoter role, or both. With
`include_promoter_units: true`, eligible promoter units are part of the same
background. Calling the denominator “all biological enhancers” would be inaccurate.

## 2. Compact total formula

```math
\boxed{
\mathrm{PACE}(E,G)=
\frac{A_\star(E)\,\overline C(E,G)\,[B(E,G)]^{\eta_{\mathrm{used}}}}
{\displaystyle\sum_{e\in\mathcal E^{\mathrm{score}}(G)}
 A_\star(e)\,\overline C(e,G)\,[B(e,G)]^{\eta_{\mathrm{used}}}}
}
```

```math
S(E,G)=A_\star(E)\,\overline C(E,G)\,[B(E,G)]^{\eta_{\mathrm{used}}}.
```

| Symbol | Meaning | Software field or source |
|---|---|---|
| E, e | Current unit and a unit traversed in the denominator | `element_id`; fixed coordinates in `units.tsv` |
| G, H | Target gene and another candidate target gene | `gene_id` |
| M, m | Fixed activity assay panel and one assay | `activity.panel` |
| x_obs,m(E) | Qualified observed quantitative signal for assay m | Aggregated `observed_activity.tsv` |
| x_hat,m(E) | Quantitative sequence-model prediction | `predictions.tsv` or local model inference |
| x_star,m(E) | Signal actually selected or fused for the assay | `resolved_activity.tsv:resolved_value` |
| A_star(E) | Activity after resolving every required assay | `scores.tsv.gz:A_used` |
| T(G), t | Distinct physical TSSs of G and one TSS | `promoters.tsv` |
| pi(t given G) | Fixed TSS weight; nonnegative and sums to one per gene | `promoters.tsv:pi` |
| C_obs, C_prior | Observed and prior contact on compatible scales | Contact tables and prior manifest |
| r(E,t) | Declared reliability of the observed contact, between 0 and 1 | Resolved contact `reliability` |
| Cbar(E,G) | TSS-weighted contact to G | `scores.tsv.gz:Cbar` |
| G(E) | Frozen candidate-gene set for E | `candidates.tsv` |
| B(E,G) | Share of E's candidate-gene contact assigned to G | `scores.tsv.gz:B` when used |
| eta_used | Allocation exponent actually used | `eta_calibration.json:eta` |
| E_score(G) | Planned candidates with computable support in this run | `scoreable`; actual normalization set |
| S(E,G) | Unnormalized support | `support` and numerically stable `log_support` |

In earlier discussions, **A-hat** denotes a prediction and **A-star** the activity
finally used. They coincide only when prediction is the selected evidence. This
software resolves each assay first and then builds activity; it does not average
an observed final activity and a predicted final activity indiscriminately.

## 3. Activity resolution

```math
A_\star(E)=\left[\prod_{m\in\mathcal M}x_{\star,m}(E)\right]^{1/|\mathcal M|}.
```

The allowed panels are ATAC, DNase, H3K27ac, ATAC+H3K27ac and DNase+H3K27ac.
For two assays this is the familiar square root of their product. For one assay it
is exactly that assay's signal. The panel is fixed for the entire run. A missing
assay does not cause a per-element change in the definition of activity.

Qualified normalized technical replicates are averaged within each biological
replicate, biological replicate means within each donor, and donor means with equal
weight for a population target. An individual target cannot pool different donors.
Raw sequencing counts must be normalized or pooled upstream under a declared protocol.

For a calibrated hybrid run, define the following assay-level fusion:

```math
x_{\star,m}(E)=s_m\left\{
\exp\left[
\lambda_m\log\left(1+\frac{x_{\mathrm{obs},m}(E)}{s_m}\right)
+(1-\lambda_m)\log\left(1+\frac{\widehat x_m(E)}{s_m}\right)
\right]-1\right\},\qquad s_m>0,\quad 0\le\lambda_m\le1.
```

Here `s_m` is the frozen scale and `lambda_m` the observed-signal weight from the
applicable assay/quality-stratum calibrator. This shift is part of the **fusion
transform**, not an added constant in activity or the PACE denominator. Calibrators
must use independently measured targets and the same signal unit, normalization,
target window and biological scope. If a valid interior fusion lacks one source,
the current resolver uses the available qualified observation, then the applicable
prediction. Without an identifiable matching calibrator, it follows this same
single-source priority; it does not guess a 50:50 blend. Explicit endpoint weights
require only the selected source. If neither source resolves, the assay remains NA.

## 4. Contact, multiple promoters and allocation

```math
\widetilde C(E,t)=r(E,t)C_{\mathrm{obs}}(E,t)
+[1-r(E,t)]C_{\mathrm{prior}}(E,t),
```

```math
\overline C(E,G)=\sum_{t\in\mathcal T(G)}\pi(t\mid G)\widetilde C(E,t),
\qquad \sum_{t\in\mathcal T(G)}\pi(t\mid G)=1,
```

```math
B(E,G)=\frac{\sum_{t\in\mathcal T(G)}\pi(t\mid G)\widetilde C(E,t)}
{\displaystyle\sum_{H\in\mathcal G(E)}\sum_{u\in\mathcal T(H)}\pi(u\mid H)\widetilde C(E,u)}.
```

Observed mode uses r=1 where observed contact is valid. Prior-only mode uses r=0.
Shrinkage uses an explicit reliability and its source. An allowed prior fallback
is recorded separately. The program never computes `0 × NA`: only a source with
positive weight is required. Near-diagonal pairs follow the explicit prior-or-NA
policy; diagonal observations are not automatically treated as reliable loops.

Physical TSSs are deduplicated, so several transcript IDs at one TSS do not multiply
its contact support. Provided TSS weights must sum to one before missing-data
filtering. Missing a positive-weight TSS makes the gene contact unresolved; weights
are not redistributed to the promoters that happen to have measurements. Gene TPM
alone cannot identify which TSS was used.

With automatic allocation, the numerical exponent is zero unless suitable
functional data pass grouped validation and the recorded stability checks. An
accepted estimate stays in [0,1]. At zero the entire B factor is skipped. See the
[short calibration specification](eta_calibration.md) for the separate estimation
procedure. This optional term does not establish a physical conservation law.

## 5. Fully expanded general formula

The following substitutes the activity product, promoter averaging and allocation
ratio into the total score. Every summation uses a fixed, explicitly recorded set.

```math
\mathrm{PACE}(E,G)=
\frac{
\left[\prod_{m\in\mathcal M}x_{\star,m}(E)\right]^{1/|\mathcal M|}
\left\{\sum_{t\in\mathcal T(G)}\pi(t\mid G)
[r(E,t)C_{\mathrm{obs}}(E,t)+(1-r(E,t))C_{\mathrm{prior}}(E,t)]\right\}
\left\{
\frac{\sum_{t\in\mathcal T(G)}\pi(t\mid G)
[r(E,t)C_{\mathrm{obs}}(E,t)+(1-r(E,t))C_{\mathrm{prior}}(E,t)]}
{\sum_{H\in\mathcal G(E)}\sum_{u\in\mathcal T(H)}\pi(u\mid H)
[r(E,u)C_{\mathrm{obs}}(E,u)+(1-r(E,u))C_{\mathrm{prior}}(E,u)]}
\right\}^{\eta_{\mathrm{used}}}
}{
\displaystyle\sum_{e\in\mathcal E^{\mathrm{score}}(G)}
\left[\prod_{m\in\mathcal M}x_{\star,m}(e)\right]^{1/|\mathcal M|}
\left\{\sum_{t\in\mathcal T(G)}\pi(t\mid G)
[r(e,t)C_{\mathrm{obs}}(e,t)+(1-r(e,t))C_{\mathrm{prior}}(e,t)]\right\}
\left\{
\frac{\sum_{t\in\mathcal T(G)}\pi(t\mid G)
[r(e,t)C_{\mathrm{obs}}(e,t)+(1-r(e,t))C_{\mathrm{prior}}(e,t)]}
{\sum_{H\in\mathcal G(e)}\sum_{u\in\mathcal T(H)}\pi(u\mid H)
[r(e,u)C_{\mathrm{obs}}(e,u)+(1-r(e,u))C_{\mathrm{prior}}(e,u)]}
\right\}^{\eta_{\mathrm{used}}}
}.
```

The branch defining each x_star is measured, calibrated hybrid, or sequence-only
as set out below. Software evaluates this expression in log space instead of
literally multiplying tiny floating-point numbers.

## 6. Mode 1: qualified observed activity and contact

For the fully observed case, x_star=x_obs and r=1:

```math
\mathrm{PACE}_{\mathrm{measured}}(E,G)=
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

This is the recommended starting point when the actual activity and contact data
are available. The software also allows measured activity with a declared contact
prior or shrinkage; use the **general formula's contact expression** in that case.
Measured does not mean that every candidate necessarily passes QC.

## 7. Mode 2: calibrated hybrid evidence

For edges whose activity is successfully fused in each assay, the full expression
is the following. Single-source fallbacks replace the corresponding product entry
with the selected x_obs or x_hat; the contact expression is unchanged.

```math
\mathrm{PACE}_{\mathrm{hybrid}}(E,G)=
\frac{
\left[\prod_{m\in\mathcal M}s_m\left\{
\exp\left[\lambda_m\log\left(1+\frac{x_{\mathrm{obs},m}(E)}{s_m}\right)
+(1-\lambda_m)\log\left(1+\frac{\widehat x_m(E)}{s_m}\right)\right]-1\right\}\right]^{1/|\mathcal M|}
\left\{\sum_{t\in\mathcal T(G)}\pi(t\mid G)[r(E,t)C_{\mathrm{obs}}(E,t)+(1-r(E,t))C_{\mathrm{prior}}(E,t)]\right\}
\left\{\frac{
\sum_{t\in\mathcal T(G)}\pi(t\mid G)[r(E,t)C_{\mathrm{obs}}(E,t)+(1-r(E,t))C_{\mathrm{prior}}(E,t)]}
{\sum_{H\in\mathcal G(E)}\sum_{u\in\mathcal T(H)}\pi(u\mid H)[r(E,u)C_{\mathrm{obs}}(E,u)+(1-r(E,u))C_{\mathrm{prior}}(E,u)]}
\right\}^{\eta_{\mathrm{used}}}
}{
\displaystyle\sum_{e\in\mathcal E^{\mathrm{score}}(G)}
\left[\prod_{m\in\mathcal M}s_m\left\{
\exp\left[\lambda_m\log\left(1+\frac{x_{\mathrm{obs},m}(e)}{s_m}\right)
+(1-\lambda_m)\log\left(1+\frac{\widehat x_m(e)}{s_m}\right)\right]-1\right\}\right]^{1/|\mathcal M|}
\left\{\sum_{t\in\mathcal T(G)}\pi(t\mid G)[r(e,t)C_{\mathrm{obs}}(e,t)+(1-r(e,t))C_{\mathrm{prior}}(e,t)]\right\}
\left\{\frac{
\sum_{t\in\mathcal T(G)}\pi(t\mid G)[r(e,t)C_{\mathrm{obs}}(e,t)+(1-r(e,t))C_{\mathrm{prior}}(e,t)]}
{\sum_{H\in\mathcal G(e)}\sum_{u\in\mathcal T(H)}\pi(u\mid H)[r(e,u)C_{\mathrm{obs}}(e,u)+(1-r(e,u))C_{\mathrm{prior}}(e,u)]}
\right\}^{\eta_{\mathrm{used}}}
}.
```

A weight is a calibrated evidence-combination parameter, not an automatically
inferred universal quality score. In particular, low-quality observations must be
marked unusable before fusion; adding a sequence prediction does not validate them.

## 8. Mode 3: reference or individual sequence evidence

For the implemented distance-prior contact model:

```math
\widehat x_m(E)=\frac1{K_E}\sum_{k=1}^{K_E}f_{\theta,m}(\mathrm{seq}_{E,k}),
\qquad C_{\mathrm{prior}}(E,t)=a
\left[\frac{\max(d_{Et},d_{\min})}{d_{\mathrm{ref}}}\right]^{-\gamma}.
```

The inputs are the validated fixed-target sequence windows for the actual copies
K_E (one or two in the implemented genotype path). f predicts a **quantitative
assay signal**, not a peak probability. The marginal signal is averaged across
copies before the geometric mean. d is the distance between the unit anchor and
TSS. The positive a, gamma, d_min and d_ref are fitted/declared in a matching prior
asset; there are no hard-coded human contact parameters.

```math
\mathrm{PACE}_{\mathrm{genome}}(E,G)=
\frac{
\left[\prod_{m\in\mathcal M}\left\{\frac1{K_E}\sum_{k=1}^{K_E}f_{\theta,m}(\mathrm{seq}_{E,k})\right\}\right]^{1/|\mathcal M|}
\left[\sum_{t\in\mathcal T(G)}\pi(t\mid G)a\left(\frac{\max(d_{Et},d_{\min})}{d_{\mathrm{ref}}}\right)^{-\gamma}\right]
\left[\frac{
\sum_{t\in\mathcal T(G)}\pi(t\mid G)a\left(\frac{\max(d_{Et},d_{\min})}{d_{\mathrm{ref}}}\right)^{-\gamma}}
{\sum_{H\in\mathcal G(E)}\sum_{u\in\mathcal T(H)}\pi(u\mid H)a\left(\frac{\max(d_{Eu},d_{\min})}{d_{\mathrm{ref}}}\right)^{-\gamma}}
\right]^{\eta_{\mathrm{used}}}
}{
\displaystyle\sum_{e\in\mathcal E^{\mathrm{score}}(G)}
\left[\prod_{m\in\mathcal M}\left\{\frac1{K_e}\sum_{k=1}^{K_e}f_{\theta,m}(\mathrm{seq}_{e,k})\right\}\right]^{1/|\mathcal M|}
\left[\sum_{t\in\mathcal T(G)}\pi(t\mid G)a\left(\frac{\max(d_{et},d_{\min})}{d_{\mathrm{ref}}}\right)^{-\gamma}\right]
\left[\frac{
\sum_{t\in\mathcal T(G)}\pi(t\mid G)a\left(\frac{\max(d_{et},d_{\min})}{d_{\mathrm{ref}}}\right)^{-\gamma}}
{\sum_{H\in\mathcal G(e)}\sum_{u\in\mathcal T(H)}\pi(u\mid H)a\left(\frac{\max(d_{eu},d_{\min})}{d_{\mathrm{ref}}}\right)^{-\gamma}}
\right]^{\eta_{\mathrm{used}}}
}.
```

This is **tissue-conditioned predicted regulatory potential**, conditional on an
applicable trained model and supplied candidate atlas. It does not measure the
animal's current chromatin state. The current package implements a distance prior,
not a sequence-based 3D contact predictor. A new individual's WGS can use an
appropriate independently trained model; no functional data anywhere in the
species means that such a model cannot be assumed. Cross-species extrapolation
requires an explicitly validated adaptation, not a renamed manifest.

## 9. Additional omics and gene expression

RNA-seq, auxiliary histone marks, methylation and CTCF enter named annotations and,
when trained and validated, a separate supervised model. They do not receive
arbitrary fixed multipliers in the primary score. This separation keeps observed
activity and structural or repressive evidence from being counted as the same
biological quantity. See [omics interfaces](MULTIOMICS.md) for inputs and commands.

A gene-only expression multiplier inside both numerator and denominator cancels
exactly for the same G. Multiplying it after normalization changes the quantity's
interpretation and cross-gene scale; this implementation deliberately keeps the
primary PACE share unchanged. Gene expression is available for separate annotation
and a predeclared classifier ablation. Expression magnitude is not a functional
positive/negative label for an enhancer–gene edge.

## 10. Boundary cases and an example

| Condition | Output |
|---|---|
| Required evidence missing or invalid | NA support and an explicit reason |
| Required evidence known, one necessary factor is exactly zero | Zero support, retained in the scoreable set |
| Positive allocation requested, required candidate-gene contact missing | Unresolved support; no silent candidate removal |
| Gene support total zero or no scoreable candidates | NA PACE |
| Only part of the planned candidate set is scoreable | `normalization_status=partial` |
| Very small or large support product | Stable log-space scoring; inspect `log_support` |

For one gene and eta=0, suppose E1 has ATAC=4, H3K27ac=9 and Cbar=2; E2 has
ATAC=1, H3K27ac=4 and Cbar=1. Their activities are 6 and 2, supports 12 and 2,
and PACE scores 12/14≈0.857 and 2/14≈0.143. The first unit has more support
within this candidate set; “85.7% chance of being causal” is not implied.

Removing E2 because its measurement is unavailable gives E1 a conditional score
of 1 with partial coverage. Its experimental evidence has not become stronger.
Use [common-background comparison](comparison.md), inspect support and coverage,
and report the exact software commit and evidence context.

Implementation: [scoring](../src/pace_livestock/core/scoring.py),
[activity resolution](../src/pace_livestock/evidence/resolve.py),
[contact prior](../src/pace_livestock/evidence/contact.py).
