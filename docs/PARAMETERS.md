> **Interface scope:** This retained guide documents the legacy region-based scripts/workflow.
> For the installable canonical-grid package and its September 2026 defaults, see [software.md](software.md).

# Parameters and their rationale

PACE (Prediction of Activity-based regulatory Connections for Enhancers) uses the same scoring rules across species. Its additional parameters describe available assays, contact reliability, alternative promoters and incomplete candidate support. They make the model particularly useful for incomplete livestock datasets; their defaults are declared starting choices rather than species-trained optima.

**Where a parameter is set matters.** The prepared-table CLI, signal-file calculator and Snakemake wrapper share a scoring kernel but expose different controls. A YAML key is not a substitute for a CLI option or a table column. The tables below describe the actual implementation.

## Start with the data condition

| Data condition | Setting to inspect | Why this setting exists | What to do for the first analysis |
| --- | --- | --- | --- |
| A planned assay is unavailable | $m$, assay JSON, $Q_A$ | Separate unmeasured input from a measured zero | Keep the assay declared and use `NA`; do not invent measurements |
| Assays use different numerical units | $a_i$, $w_i$ | Separate scale from assay importance | Choose fixed positive scales independently; start with equal assay priors |
| Contact is sparse or from a surrogate | $H$, $D$, $\lambda$, source | Let qualified observations change a distance prior in proportion to reliability | Use the prior when QC is unknown; provide matched expectations and independent reliability when available |
| Several transcripts share or change promoters | $\pi(G,t)$ | Prevent transcript count from multiplying support | Deduplicate TSSs, then assign uniform weights over the full available catalogue unless independent use data exist |
| One enhancer has several candidate genes | $B(E,G)$, $\eta$ | Account for enhancer-side target sharing | Use $\eta=1$; compare with 0 as an allocation ablation |
| Evidence or candidate coverage is incomplete | $Q$, $\mathcal E^{\mathrm{obs}}$, $U$ | Keep score, input sufficiency and unscored support distinct | Preserve unfiltered output; leave unsupported QC and residual support unknown |

These are the model extensions explained in the [six-part ABC comparison](ABC_COMPARISON.md). The cis window, candidate width, peak cap and output cutoff are preparation or selection settings, not livestock-specific innovations. [Worked examples](WORKED_EXAMPLES.md) show the effect of changing input availability and allocation.

## 1. Activity integration

| Quantity / setting | Default or range | How to set it | Reason and effect |
| --- | --- | --- | --- |
| Aggregation | `missing_geometric` | Signal YAML `activity_method`; neighborhood `--activity_method` | Shifted geometric integration keeps measured zero distinct from missing input; other aggregation modes are rejected by the primary adapters |
| Planned assays | Explicit JSON signal list; accessibility plus supplied H3K27ac in direct-file runs | `--activity-config activity.json`; `signals` entries | Defines which evidence is expected; leave a planned but unmeasured assay in JSON with `NA` values to count its absence in quality assessment |
| Assay prior $w_i$ | 1 for accessibility and H3K27ac; nonnegative, positive total | JSON `signals.<name>.weight`; direct neighborhood `--accessibility_weight`, `--H3K27ac_weight`; standalone YAML assay weights | Equal default importance avoids claiming fitted livestock weights; changing priors changes activity and requires a declared sensitivity analysis |
| Assay scale $a_i$ | 1 in example JSON; finite and positive | JSON `signals.<name>.scale` | Makes assay units explicit before `log1p`; choose scales on reference/training data and hold them fixed during testing and downsampling |
| Observation mask $m(E,i)$ | 1 for finite observations, 0 for missing | Derived from signal entries | A measured 0 remains observed; `NA` means no measurement. All effective observations missing gives undefined activity |
| Local assay quality $q(E,i)$ | Unknown unless supplied; [0,1] | JSON `quality_column` naming a table column | Allows an independently assessed assay quality to change its contribution. Unknown uses unit weight only for provisional activity; output quality remains unknown |
| Log shift | `log(1+x)` and inverse `exp(x)-1` | Fixed in `aggregate_activity` | Allows a positive layer to coexist with a zero layer while keeping an all-zero activity at zero; may alter low-signal false-positive behavior |
| Inhibitory strength $\kappa$ | 0; nonnegative | Python API `aggregate_activity(..., inhibition_strength=...)` | Optional attenuation is disabled in primary predictions because fraction-scale inhibitory evidence and validation are needed |
| Inhibitory evidence $I(E)$ | Weighted mean of supplied [0,1] values | Python API `inhibitory`, `inhibitory_weights` | Allows explicitly scaled extensions. Methylation/repressive-mark switches in the raw-file adapters are not a supported route to this extension |

Assay scales are not estimated by `scripts/pace.py`. A median positive signal over a prespecified reference candidate set is one possible scaling rule, not a universal default. Record the chosen rule and numerical values. Raw-file adapters produce reads per kb or bigWig summaries; they do not automatically make different assays, depths or tissues comparable. Use prepared tables when fixed scales and local QC are required.

## 2. Contact, promoters and score normalization

| Quantity / setting | Default or range | How to set it | Reason and effect |
| --- | --- | --- | --- |
| Candidate window | Strict distance <5,000,000 bp, same chromosome | Predictor `--max_distance`; standalone YAML `params_predict.window` | Keeps a broad cis candidate background; increasing it changes computation and score denominators. The table CLI scores the pairs supplied and does not apply a window |
| Prior exponent $\gamma$ | 1.024238616787792 | Predictor `--hic_gamma`; standalone/workflow YAML `params_predict.hic_gamma` | Controls contact decay with distance. Retained starting value; refit only on an independent, compatible reference profile |
| Prior scale | 5.9594510043736655 | Predictor `--hic_scale`; standalone/workflow YAML `params_predict.hic_scale` | Sets raw contact units. A common multiplicative scale cancels from the normalized score when residual support is zero; it still matters for raw components and nonzero residual units |
| Prior distance offset | 5,000 bp | Fixed default in the current file-adapter prediction path; custom `contact_prior` in table input | Prevents a singularity at distance zero. The retained YAML `hic_pseudocount_distance` key is not forwarded by that path |
| Observed contact $H$ | Missing | `contact_observed` column or a supported contact file | Supplies measured contact; zero can be a valid observation and must not be substituted for an absent edge |
| Expected contact $D$ | Missing; positive when supplied | `contact_expected` column or `--contact_metadata` | Places observed contact on a compatible distance-dependent scale. Estimate from eligible genomic bin pairs with matching sample, resolution and normalization |
| Reliability $\lambda$ | 0 used for contact mixing when prerequisites are missing; supplied values in [0,1] | `contact_reliability` column / metadata | Shrinks uncertain observations toward the prior. It is QC, not the observed contact count or a parameter fitted to positive labels |
| Contact source | `distance_prior` if unspecified | `contact_source`: `matched`, `surrogate`, `distance_prior`, `unknown` | Permits measured-contact mixing only for declared matched/surrogate sources and records provenance |
| Matrix resolution | 5,000 bp in file-reader CLI | `--hic_resolution`; sample-sheet `HiC_resolution` | Must match the data. The Snakemake sample validator currently requires 5 kb; use direct commands for other supported resolutions |
| TSS-use weight $\pi(G,t)$ | Equal across distinct annotated TSSs | `tss_weight` in genes or pair table | Avoids multiplying support by transcript count; independent promoter-use evidence can replace uniform weights. Assign over the full gene catalogue before windowing |
| Target allocation $B(E,G)$ | Contact share over the enhancer's candidate genes | Computed from candidate contacts | Represents competing candidate targets; incomplete annotation can bias the allocation |
| Allocation exponent $\eta$ | 1; allowed range [0,1] | Table CLI `--competition-power`; Python `ScoreConfig` | 1 uses full allocation, 0 removes this component for an ablation. Raw-file adapters currently use 1 |
| Residual support $U(G)$ | Unknown, treated as 0 computationally | `unassigned_mass` column, constant per gene | Allows an independent estimate of omitted support in raw-support units. No automatic estimator is supplied; an arbitrary missing-peak fraction is not a compatible value |

See [equations](FORMULA.md) for the order of operations. Keep promoter-proximal candidates in the scoring background when comparing distal links. A display filter should be applied after scoring.

## 3. Evidence, RNA and output selection

| Setting | Default / control | Why it exists and how to interpret it |
| --- | --- | --- |
| Activity quality | Computed from declared assay qualities, or supplied `activity_quality` | Distinguishes availability and quality from the size of the activity signal |
| TSS quality | `tss_quality`; unknown by default | Marks independently assessed promoter annotation; more transcripts are not evidence of better quality |
| Candidate coverage quality | `catalogue_quality`; unknown by default | Requires an independent catalogue-coverage assessment; cannot be inferred from the fraction of retained predictions |
| Sufficient-evidence threshold | 0.5; table CLI `--evidence-threshold` or `ScoreConfig` | All four quality components must meet it and the gene must have no unscored supplied candidates. It is an operational QC rule, not a statistical confidence cutoff |
| Expression eligibility | Predictor `--min_expression 1.0` TPM | Annotates `detected` / `below_threshold` / `unknown`; no expression multiplier enters the score. The standalone adapter currently uses its predictor default |
| Score cutoff | 0.02; `--threshold` for file predictors/filter; workflow `params_filter_predictions.threshold` | Produces a convenient filtered table. It is not a livestock-calibrated FDR threshold and must be chosen on independent validation for decision use |
| Expression-only output | Off; filter `--only_expressed` or workflow `only_expressed_genes` | An explicit downstream selection. Unknown expression must not silently be presented as measured zero |
| Resampling summaries | API requires ≥2 independently rescored runs; 5%, 50%, 95% quantiles | Describes stability across actual reruns; not a causal confidence interval. Missing edges remain in the scoring-frequency denominator |

Prior-only predictions have contact-observation quality 0. With the default evidence threshold, they are normally `provisional` even when their structural scores are finite. Lowering the threshold changes the designation, not the underlying measurement support.

## 4. Candidate preparation and workflow settings

These are practical starting settings rather than claimed model innovations.

| Setting | Default / location | Why and when to change |
| --- | --- | --- |
| Summit extension | ±250 bp; `--peakExtendFromSummit`, YAML `params_candidate.peakExtendFromSummit` | Produces approximately 500-bp starting candidates around accessible summits; overlapping intervals are not automatically merged |
| Minimum region width | 500 bp; direct candidate CLI `--minPeakWidth` | Avoids very narrow signal-counting intervals; boundary handling affects final widths |
| Peak cap | 150,000; `--nStrongestPeaks`, YAML `params_candidate.nStrongestPeaks` | Limits the candidate background and resource use; this implementation ranks narrowPeak `signalValue`. Evaluate coverage if changing the cap |
| MACS2 significance | `p=0.1`; YAML `params_macs.pval` | Permissive candidate discovery before ranking/filtering; not the significance of an enhancer–gene association |
| MACS2 effective genome size | `2.5e9` placeholder in generic YAML | Replace with an independently justified effective genome size for the exact assembly/mapping strategy; do not reuse it for every species |
| MACS2 shift / extension | `--shift -75 --extsize 150 --nomodel`, fixed workflow settings | Assumes the declared accessibility-read processing; for different fragment protocols call peaks upstream and use the direct route |
| MACS2 duplicate policy | `--keep-dup all`, fixed workflow setting | Requires a documented upstream duplicate/QC policy; the workflow is not a raw-read QC pipeline |
| Blacklist | User-supplied BED | Remove artifact-prone regions for the same assembly; no human blacklist is automatically appropriate for livestock |
| Promoter candidate regions | User-supplied BED | The adapter only marks overlaps; add promoter candidates upstream when required and keep the same policy for compared methods |
| Parallel jobs | Snakemake `--cores` | Controls scheduling. The retained `params_macs.threads` key is not passed to MACS2 as a threading option |

## 5. Which configuration entries are active?

Use the prepared-table JSON for assay scales, local assay QC and model ablations. Use standalone YAML for sample file paths, enabled activating assays, their weights, `params_predict.window`, gamma and scale. Use the predictor CLI for custom candidate windows, expression annotation and contact metadata.

The Snakemake wrapper forwards gamma/scale and sample metadata, but does **not** forward arbitrary `ScoreConfig` fields, `params_predict.window`, per-assay YAML weights or the distance-offset key. The direct neighborhood command has its own weight options. `params_neighborhoods.use_qnorm`, `params_predict.flags`, `output_options` and compatibility ML/eQTL sections do not activate new primary scoring features. Disabling `include_self_promoter` in YAML is not a substitute for explicitly controlling the candidate background. The primary model does not require or run the archived ML module.

To change a model component that a wrapper does not expose, use the documented table interface or core API and preserve that choice with the output. Do not add an unused YAML key and assume the predictions changed.

## 6. Parameter justification and reproducibility

Record assembly, candidate/TSS catalogue, assay units/scales, contact normalization, QC procedure, window and filtering threshold for every run. Choose defaults before examining test outcomes. Test sensitivity to the uncertain settings against independently held-out data and report candidate/positive coverage alongside score-based metrics. No quality threshold, prior exponent or species-specific optimum is learned by this software automatically.

Source of truth: [core kernel](../workflow/scripts/pace_core.py), [table CLI](../scripts/pace.py), [predictor](../workflow/scripts/predictor.py), [candidate construction](../workflow/scripts/peaks.py), and [workflow rules](../workflow/rules/). Scientific comparisons and attribution are in [ABC comparison](ABC_COMPARISON.md).
