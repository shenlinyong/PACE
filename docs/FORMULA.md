# PACE formulas and interpretation

[中文公式](FORMULA.zh-CN.md) · [Practical preparation](PRACTICAL_WORKFLOW.md) · [Parameters](parameters.md)

PACE requires measured ATAC, DNase or H3K27ac activity. Its primary output is a relative support share within a fixed candidate universe, not a causal probability, a gene-expression fraction, or a calibrated score comparable across all genes.

## Total formula

```math
\boxed{
\mathrm{PACE}(E,G)=
\frac{A_\star(E)\,\overline C(E,G)\,[B(E,G)]^{\eta_{\mathrm{used}}}}
{\displaystyle\sum_{e\in\mathcal E(G)}
 A_\star(e)\,\overline C(e,G)\,[B(e,G)]^{\eta_{\mathrm{used}}}}
}
```

The denominator includes the **full planned candidate set** E(G). By default, the software publishes `pace_score` only when all planned supports are resolved and the total is positive. Unresolved candidates are retained. For an incomplete gene, `pace_score` is NA; the separate `pace_score_conditional` describes the available subset. This prevents incomplete-background scores from masquerading as full-background scores. Complete still means the planned catalog, not every biological enhancer.

| Symbol | Meaning | Input/output |
|---|---|---|
| E, G | Current scoring unit and target gene | element_id, gene_id |
| e | Each planned candidate for G | candidates.tsv |
| M, m | Fixed assay panel and one assay | activity.panel |
| x_obs,m | Qualified normalized experimental signal after replicate aggregation | resolved_activity.tsv |
| A_star | Equal geometric mean of the fixed assay panel | A_used |
| T(G), t | Distinct physical TSSs and one TSS | promoters.tsv |
| pi(t given G) | Frozen promoter weights summing to one per gene | pi |
| Ctilde | Contact after the declared correction/prior policy | resolved_contacts.tsv |
| Cbar | Sum of pi times Ctilde over promoters | Cbar |
| G(E), H | Fixed candidate-gene set for an element and one gene | candidates.tsv |
| B | Cbar divided by its sum over G(E) | B |
| eta_used | Actual optional allocation exponent | eta_calibration.json |
| E_score(G) | Subset whose support can be resolved | scoreable |

## Fully expanded formula

```math
\mathrm{PACE}(E,G)=
\frac{
\left[\prod_{m\in\mathcal M}x_{\mathrm{obs},m}(E)\right]^{1/|\mathcal M|}
\left[\sum_{t\in\mathcal T(G)}\pi(t\mid G)\widetilde C(E,t)\right]
\left[
\frac{\sum_{t\in\mathcal T(G)}\pi(t\mid G)\widetilde C(E,t)}
{\sum_{H\in\mathcal G(E)}\sum_{u\in\mathcal T(H)}\pi(u\mid H)\widetilde C(E,u)}
\right]^{\eta_{\mathrm{used}}}
}{
\displaystyle\sum_{e\in\mathcal E(G)}
\left[\prod_{m\in\mathcal M}x_{\mathrm{obs},m}(e)\right]^{1/|\mathcal M|}
\left[\sum_{t\in\mathcal T(G)}\pi(t\mid G)\widetilde C(e,t)\right]
\left[
\frac{\sum_{t\in\mathcal T(G)}\pi(t\mid G)\widetilde C(e,t)}
{\sum_{H\in\mathcal G(e)}\sum_{u\in\mathcal T(H)}\pi(u\mid H)\widetilde C(e,u)}
\right]^{\eta_{\mathrm{used}}}
}.
```

This expression uses the same resolved-contact function in every numerator and denominator term. At eta=0, the entire allocation factor is omitted, including its data requirements.

## Measured activity

```math
A_\star(E)=\left[\prod_{m\in\mathcal M}x_{\mathrm{obs},m}(E)\right]^{1/|\mathcal M|}.
```

Panels: ATAC, DNase, H3K27ac, ATAC+H3K27ac, DNase+H3K27ac. A missing required assay gives NA; a measured zero remains zero. A panel cannot change per element. First average normalized technical replicates within biological replicates, then biological replicates within donors, then donors equally. `individual` requires one donor; `population_mean` explicitly pools donors.

For raw window counts, `pace normalize-activity` implements count × 10^6 / filtered-library-size / window-length. Already normalized bigWigs must not be normalized a second time. Library-size correction is not batch correction; animal, sex, age, tissue and protocol remain part of study design.

## Contact: prior, near diagonal and sparse counts

```math
P(d)=a\left[\frac{\max(d,d_{\min})}{d_{\mathrm{ref}}}\right]^{-\gamma},\qquad
p(d)=\kappa\min\{P(d),P(d_0)\}.
```

For off-diagonal measured contacts, the regularized form is:

```math
\widetilde C(E,t)=C_{\mathrm{obs}}(E,t)+p(d(E,t)).
```

`pseudocount: auto` adds this term only in observed mode when a compatible supplied prior exists. Otherwise it leaves the measurement unchanged and records that no prior was available. `powerlaw` requires a compatible prior; `none` disables it. Defaults: kappa=1, d0=5000 bp. These are algorithm settings, not demonstrated livestock optima. Missing activity is never filled.

| Situation | Resolved contact |
|---|---|
| Same bin or declared near-distance range, with compatible prior | P(d), without adding another pseudocount |
| Same bin, no prior, default prior_or_neighbor policy | Recorded maximum of valid neighboring contacts from the cooler adapter, if available |
| Near diagonal with no eligible replacement | NA; partial primary scores are withheld |
| Off-diagonal finite observation, including a sampled zero | Observation plus declared pseudocount, or the unchanged observation |
| Missing contact with allow_prior_fallback=true and matching prior | P(d), explicitly labeled as a prior |
| Missing contact without allowed fallback | NA |
| Explicit prior_only | P(d) |
| Explicit shrinkage | r C_obs + (1-r) P(d); do not add a second pseudocount |

Near-diagonal handling precedes ordinary contact-mode selection. `unresolved` deliberately retains near-diagonal NA. `prior_or_unresolved` is the earlier strict option; the default `prior_or_neighbor` also accepts the adapter's recorded neighbor correction. Same-bin geometry depends on resolution and bin boundaries, not just a distance less than 5 kb. A neighborhood maximum is a local contact surrogate, not an experimentally validated promoter loop.

Raw `observed_value`, `prior_value`, `pseudocount_value`, evidence type, prior identity and resolution reason remain in the output. `regularized` distinguishes additive prior regularization from observations and from convex shrinkage. Its reliability field is the retained observation coefficient (1), not a posterior confidence. A contact count of zero is a sampling observation, not proof of absent biological interaction. Invalid balanced bins remain unavailable unless an explicitly declared prior replaces them.

A prior mixed with measured contacts must match scale, normalization, resolution, balancing and measurement windows. Fit one directly using `pace fit-prior --cooler ...`. All callable bin opportunities, including unstored zeros, enter distance-bin means. Zero-mean bins cannot enter the logarithmic regression and are reported. Optional held-out chromosomes diagnose contact-decay fitting; they are not functional validation.

An explicit `prior_preset: abc_human` with `mode: prior_only` uses the published human ABC gamma=1.024238616787792, with a=1 and d_ref=d_min=5000 on a **relative** scale. Amplitude cancels in this prior-only score. This is an unvalidated transferred baseline, not a fitted livestock prior; it is rejected in validated/demonstration profiles and cannot be mixed with measured contact tables. No measured-activity requirement is relaxed.

## Multiple TSSs and optional allocation

```math
\overline C(E,G)=\sum_{t\in\mathcal T(G)}\pi(t\mid G)\widetilde C(E,t),\qquad
B(E,G)=\frac{\overline C(E,G)}{\sum_{H\in\mathcal G(E)}\overline C(E,H)}.
```

Deduplicate identical physical TSSs. Distinct TSSs in the same measured bin can reuse an available query for the same element and sample, without adding experimental replicates. Biological TSS identities and summed weight remain intact; `n_contact_bins` reports spatial resolution. The exact-distance prior may differ between TSSs within one bin. Missing a genuinely distinct positive-weight contact still gives NA; weights are not renormalized after missingness.

`pace prepare-promoter-weights` can freeze proportional ATAC, DNase, H3K4me3 or CAGE promoter-signal weights. These are assay-supported proxies, not necessarily transcription initiation fractions. All planned TSSs need measured signals. All-zero genes error unless an explicit equal-weight policy is requested. Gene TPM cannot identify promoter usage.

The allocation term satisfies Cbar × B^eta = Cbar^(1+eta) / (sum_H Cbar(E,H))^eta. Thus annotation and candidate-gene density matter. It is an extension, not a proven biological competition law. Default auto falls back to zero without independent grouped functional evidence. Benchmarks include `contact_power_2` as a control for merely squaring contact, alongside eta=0 and eta=1; final comparisons need independently selected settings and frozen gene annotations.

## Missing denominator: conditional score and sensitivity bounds

Let S_i be the nonnegative support of one candidate. For available support:

```math
\mathrm{PACE}_{\mathrm{conditional}}(E,G)=
\frac{S(E,G)}{\sum_{e\in\mathcal E^{\mathrm{score}}(G)}S(e,G)}.
```

`scoring.partial_policy: conditional` explicitly restores this quantity in the primary column for legacy exploration; `score_scope` still identifies it as conditional. The default is `withhold`. Cross-run comparisons recompute common denominators and separate conditional differences.

For support assumptions L_i <= S_i <= U_i:

```math
\mathrm{PACE}_{i,\mathrm{lo}}=\frac{L_i}{L_i+\sum_{j\ne i}U_j},\qquad
\mathrm{PACE}_{i,\mathrm{hi}}=\frac{U_i}{U_i+\sum_{j\ne i}L_j}.
```

Resolved supports have exact bounds. Unresolved supports default to [0,infinity); optional `inputs.support_bounds` supplies justified bounds on final A×C×B^eta support with a `bound_source`. These bounds do not impute activity or make the primary score complete. Empty `support_upper` means unbounded. Results are `pace_score_lo/hi`; they are sensitivity/identification ranges conditional on positive total support, **not statistical confidence intervals**. Correlated unknown supports can make the rectangular bounds conservative. If no positive total is possible, return NA; a sole potentially positive candidate has [1,1], which establishes no biological validity.

Example: observed supports 2 and 1, with a third unknown support. The first conditional share is 2/3; its full-background range is [0,2/3]. If independent assumptions constrain the missing support to [2,4], the range becomes [2/7,2/5]. Without justified upper limits, the software does not invent a narrow range.

## Region summaries and other omics

A peak can overlap several unique 500 bp grid cells. `region_scores.tsv` sums each member cell once for each source/region/gene; it does not create a new scoring denominator. Overlapping regions must not be summed together as independent elements. This does not reproduce ABC's summit-centered candidate construction. Regional interval sums are conservative bounds, capped at one.

RNA, H3K4me1/3, H3K27me3, H3K9me3, CTCF and methylation retain annotation and independently trained classifier interfaces. No generic expression multiplier or inhibitory constant is added to the core formula. H3K4me3/ATAC promoter weights are an explicit separate preparation choice. Software checks do not establish increased biological accuracy.

## Reference implementations

The contact pseudocount follows the minimum of distance expectation and the expectation at a fixed cap distance in the [ABC predictor](https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction/blob/main/workflow/scripts/predictor.py). Defaults are documented in the [ABC configuration](https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction/blob/main/config/config.yaml). PACE does not claim byte-for-byte equivalence to ABC: candidate construction, TSS handling, fallback policies and missing-denominator reporting differ. See [ABC comparison](ABC_COMPARISON.md) and the [validation design](ROBUSTNESS_VALIDATION.md).
