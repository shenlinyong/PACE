# PACE model

[中文公式](FORMULA.zh-CN.md) · [Data preparation](PRACTICAL_WORKFLOW.md) · [Settings](parameters.md)

PACE uses measured activity and promoter-contact evidence to rank candidate regulatory elements for each gene. Scores are relative support shares, not causal probabilities or fractions of gene expression.

## Main formula

```math
\boxed{
\mathrm{PACE}(E,G)=
\frac{A_\star(E)\,\overline C(E,G)}
{\displaystyle\sum_{e\in\mathcal E(G)} A_\star(e)\,\overline C(e,G)}
}
```

The denominator contains the **complete planned candidate set** E(G), including promoter units by default. A positive, fully resolved denominator is required for `pace_score`. If a candidate is unresolved, the primary score is NA and the available-subset score is reported separately as `pace_score_conditional`. Complete refers to the input catalog, not every biological enhancer.

| Symbol | Meaning | Input or output |
|---|---|---|
| E, G | Current element and target gene | element_id, gene_id |
| e, E(G) | One candidate and the complete candidate set for G | candidates.tsv |
| M, m | Fixed assay panel and one assay | activity.panel |
| x_m(E) | Measured signal after replicate aggregation | resolved_activity.tsv |
| A_star(E) | Geometric mean activity | A_used |
| T(G), t | Distinct physical TSSs and one TSS | promoters.tsv |
| pi(t given G) | TSS weight; sums to one per gene | promoter_weights.tsv |
| Ctilde(E,t) | Processed contact to one TSS | resolved_contacts.tsv |
| Cbar(E,G) | Weighted contact over the selected TSSs | Cbar |

## Measured activity

```math
A_\star(E)=\left[\prod_{m\in\mathcal M}x_m(E)\right]^{1/|\mathcal M|}.
```

Supported panels are ATAC, DNase, H3K27ac, ATAC+H3K27ac and DNase+H3K27ac. The panel is fixed across elements. A missing required assay gives NA; a measured zero gives zero activity by default. Technical replicates are averaged within biological replicates, biological replicates within donors, and donors equally. Signals must be normalized comparably before averaging. `individual` requires one donor; `population_mean` pools donors explicitly.

Multiplying every element's signal for assay m by the same positive constant multiplies all activities by a common factor. It cancels in a gene's normalized score. This does not remove differences in signal-to-noise, measurement windows, missingness or replicate-specific scaling.

For shallow measurements, users may set an assay-specific nonnegative pseudocount:

```math
A_{\star,\epsilon}(E)=
\left[\prod_{m\in\mathcal M}(x_m(E)+\epsilon_m)\right]^{1/|\mathcal M|}.
```

`activity.pseudocounts` defaults to `{}` (all epsilon values zero). Pseudocounts are in the units of the normalized assay signal and are applied after replicate aggregation. Raw resolved measurements are preserved; `activity_pseudocount` records the addition. NA remains NA. Rescaling an assay preserves scores only if its pseudocount is rescaled too. A fixed positive offset can increase background support; choose it independently of test outcomes and report sensitivity to zero.

`pace normalize-activity` converts raw counts to count × 10^6 / filtered-library-size / window-length. Already normalized bigWigs must not be normalized twice. Library-size normalization is not batch correction.

## Contact evidence

```math
P(d)=a\left[\frac{\max(d,d_{\min})}{d_{\mathrm{ref}}}\right]^{-\gamma},\qquad
p(d)=\kappa\min\{P(d),P(d_0)\}.
```

| Situation | Contact used |
|---|---|
| Same bin or configured near-distance range, compatible prior available | P(d) |
| Same bin, no prior, default prior_or_neighbor policy | Maximum of valid neighboring contacts recorded by the cooler adapter |
| Near diagonal with no eligible replacement | NA |
| Finite off-diagonal observation, including zero | Observed contact plus the configured pseudocount |
| Missing observation with allow_prior_fallback=true and compatible prior | P(d), labeled as a prior |
| Missing observation without allowed fallback | NA |
| Explicit prior_only mode | P(d) |
| Explicit shrinkage mode | r C_obs + (1-r) P(d) |

Near-diagonal handling runs first. Same-bin membership depends on bin boundaries and resolution; it is not equivalent to every pair within a fixed ±5 kb distance. `near_diagonal_policy: unresolved` leaves these contacts unavailable. `prior_or_unresolved` permits only a matching prior. The default `prior_or_neighbor` also accepts the recorded neighbor correction. A local maximum is a contact surrogate, not a measured regulatory loop.

`contact.pseudocount: auto` adds p(d) only to finite, off-diagonal observations in observed mode when a compatible prior is available. `powerlaw` requires that prior; `none` disables the addition. The defaults are kappa=1 and d0=5000 bp; these have not been optimized for livestock. Pseudocounts are not added again to near-diagonal replacements, prior-only contacts or shrinkage. Invalid balanced bins remain unavailable unless an explicitly allowed prior replaces them.

Outputs retain `observed_value`, `prior_value`, `pseudocount_value`, evidence type, prior identity and processing reason. A regularized contact is labeled `regularized`; its retained observation coefficient is not a confidence probability.

### Fitting and transferring a prior

`pace fit-prior --cooler ...` estimates decay from distance-bin means, including all callable zero pixels. The default lower fitting distance is the matrix resolution. Zero-mean bins cannot enter the log fit and are reported. Optional held-out chromosomes diagnose contact-decay fit; they do not validate enhancer–gene predictions. A prior mixed with observed contacts must match their scale, normalization, resolution, balancing and measurement window.

A different tissue normally fails the context check. `contact.allow_cross_context_prior: true` permits a prior from the same species and assembly, with the same target level. Its original context is retained, the target context is recorded, and the transfer is marked unvalidated. It cannot be used in the `validated` profile. Renaming a scale or tissue does not calibrate a prior; assess transferred decay on target data when possible.

The optional `prior_preset: abc_human` requires `mode: prior_only` and the research profile. It uses the [ABC reference configuration](https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction/blob/main/config/config.yaml) gamma=1.024238616787792, with a=1 and d_ref=d_min=5000 on a relative scale. This amplitude cancels in prior-only scores. It is a transferred human baseline, not a fitted livestock parameter, and cannot be mixed with measured contact tables. Measured activity is still required.

## Multiple TSSs

```math
\overline C(E,G)=\sum_{t\in\mathcal T(G)}\pi(t\mid G)\widetilde C(E,t).
```

Identical physical TSSs are deduplicated. Distinct TSSs in the same measured bin may reuse the same element/sample contact query without adding replicates. The distance prior may differ at their exact coordinates. By default all positive-weight TSSs are required; missing one gives NA. Zero-weight TSSs do not require contact.

`pace prepare-promoter-weights` calculates proportional weights from measured ATAC, DNase, H3K4me3 or CAGE promoter signals. These weights are proxies for promoter usage. Gene TPM alone cannot determine TSS usage.

Optional filtering uses **one retained set per gene for every candidate**:

- `promoters.minimum_weight` removes TSSs below the original weight threshold (default 0).
- `promoters.missing_policy: drop_missing` removes a TSS if any planned candidate lacks a resolved contact to it. The default `strict` keeps it.
- The retained original weight must reach `minimum_retained_weight` (default 0.9); otherwise the gene is left unscoreable.

```math
\pi_{\mathrm{used}}(t\mid G)=
\frac{\pi(t\mid G)}{\sum_{u\in\mathcal T_{\mathrm{keep}}(G)}\pi(u\mid G)},
\qquad t\in\mathcal T_{\mathrm{keep}}(G).
```

The original candidates and promoter units remain in the denominator catalog. `promoter_weights.tsv` records original/effective weights, dropped TSSs, missing-contact counts and retained weight. `tss_contact_scope=selected_tss_set` identifies scores for a filtered promoter definition. They do not recover regulation through excluded promoters. Bounds refer to this selected definition too. Different retained sets cannot be treated as full biological score changes in `pace compare`. With sparse maps, dropping every TSS touched by a missing contact may leave too little weight; use a justified contact prior or retain the default NA result.

## Experimental target allocation

The optional extension adds a cross-gene contact share:

```math
B(E,G)=\frac{\overline C(E,G)}{\sum_{H\in\mathcal G(E)}\overline C(E,H)},\qquad
\mathrm{PACE}_{\eta}(E,G)=
\frac{A_\star(E)\overline C(E,G)B(E,G)^\eta}
{\sum_{e\in\mathcal E(G)}A_\star(e)\overline C(e,G)B(e,G)^\eta}.
```

At eta=0 the allocation term and its data requirements are omitted, giving the main formula. `allocation.eta: auto` uses zero without eligible independent functional labels. A fixed nonzero value is an explicit experimental choice. At nonzero eta, missing contact to any candidate target prevents resolving B for that element; candidate genes are not silently removed.

Cbar × B^eta equals Cbar^(1+eta) / (sum_H Cbar)^eta. This strengthens contact contrasts and depends on gene annotation density. It is not an established biological competition law. Compare it against eta=0 and the separate `contact_power_2` control using held-out functional data. See [calibration](eta_calibration.md).

## Incomplete denominators

Writing S=A_star×Cbar (or its explicitly selected experimental extension),

```math
\mathrm{PACE}_{\mathrm{conditional}}(E,G)=
\frac{S(E,G)}{\sum_{e\in\mathcal E^{\mathrm{score}}(G)}S(e,G)}.
```

`scoring.partial_policy: conditional` also writes this quantity into the primary column for compatibility; `score_scope` still marks it conditional. The default is `withhold`. An incomplete gene with observed supports 2 and 1 has a conditional first score of 2/3. An unseen support of 97 would make its full score 0.02.

For nonnegative support ranges L_i <= S_i <= U_i:

```math
\mathrm{PACE}_{i,\mathrm{lo}}=\frac{L_i}{L_i+\sum_{j\ne i}U_j},\qquad
\mathrm{PACE}_{i,\mathrm{hi}}=\frac{U_i}{U_i+\sum_{j\ne i}L_j}.
```

Resolved supports are fixed. Unresolved supports default to [0,infinity); `inputs.support_bounds` accepts justified limits and a `bound_source`. Bounds do not impute a primary score. For supports (2,1,unknown), the first full-score range is [0,2/3]; bounding the unknown support to [2,4] narrows it to [2/7,2/5]. These are sensitivity ranges conditional on positive total support, **not confidence intervals**. A zero total yields NA. Correlated unknowns can make these ranges conservative.

## Candidates and other omics

The grid covers supplied regions and required promoters, not the whole genome. Duplicate or overlapping peaks do not duplicate cells. Broad or noisy input regions still enlarge the denominator; peak quality remains an upstream responsibility. `region_scores.tsv` sums unique member cells without changing the scoring denominator; overlapping region scores must not be added as independent units. This catalog is not identical to ABC's peak-centered catalog, so published ABC thresholds should not be transferred directly.

RNA, additional histone marks, CTCF and methylation are annotations or inputs to a separate classifier. Multiplying both numerator and denominator by the same gene-expression weight cancels algebraically. Multiplying only after normalization changes the meaning of the score.

References: [Fulco et al., 2019](https://doi.org/10.1038/s41588-019-0538-0), [Nasser et al., 2021](https://doi.org/10.1038/s41586-021-03446-x), [ABC methods](https://abc-enhancer-gene-prediction.readthedocs.io/en/stable/usage/methods.html). Human benchmark results do not establish livestock accuracy.
