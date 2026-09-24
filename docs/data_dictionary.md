# Input and output data dictionary

Tables are UTF-8 TSV (optionally gzip) with headers. BED/internal intervals are 0-based half-open; tss0 is a 0-based position. Use NA for unavailable values and 0 for true measurements of zero. All coordinates and experimental metadata must use the specified assembly and context.

## Inputs

units, promoters and candidates must be nonempty. Activity must be available from qualified experimental rows or valid imports of experimentally aggregated activity. Register every real sample and source. Additional tables are required only when referenced; do not fabricate sample identities.

### units

```text
element_id	chrom	start	end	anchor0	element_roles	canonical_catalog_id
```

### region_membership

```text
region_id	element_id	source_id	membership_rule
```

### promoters

```text
gene_id	promoter_id	chrom	tss0	strand	pi	pi_source
```

### candidates

```text
element_id	gene_id	candidate_universe_id
```

### samples

```text
sample_id	donor_id	assay	biological_replicate	technical_replicate	species	assembly	context_id	source_id
```

### observed_activity

```text
element_id	sample_id	assay	signal	measurement_status	callable_fraction	unit	normalization_id	window_id
```

### observed_contacts

```text
element_id	promoter_id	sample_id	contact_value	measurement_status	bin_pair_id	scale	resolution	source_id
```

### resolved_activity

```text
element_id	assay	observed_value	resolved_value	evidence_id	evidence_type	observation_sample_id	parent_evidence_ids	resolution_status	reason	unit	normalization_id	window_id
```

### resolved_contacts

```text
element_id	promoter_id	resolved_value	evidence_id	evidence_type	observation_sample_id	prior_id	reliability	resolved_mode	bin_pair_id	resolution_status	reason	scale
```

### features

```text
entity_type	entity_id	feature_name	value	evidence_id	status
```

### methylation

```text
chrom	dyad_start0	methylated_count	total_count	sample_id	assay
```

### expression

```text
gene_id	sample_id	tpm	status
```

### labels

```text
label_id	assayed_region_id	gene_id	context_id	perturbation_type	effect_direction	effect_size	label_status	assay_id	group_id	source_id
```

### evidence

```text
evidence_id	evidence_type	source_id	parent_evidence_ids	model_id	unit	processing_method	checksum
```

### sources

```text
source_id	path_or_accession	source_type	assembly	processing_method	normalization_id	checksum
```

### support_bounds

```text
element_id	gene_id	support_lower	support_upper	bound_source
```

Optional bounds on final unnormalized support of unresolved edges. NA upper bound means unbounded. These are explicit sensitivity assumptions, not inferred experiments or confidence limits.

## Status and measurement definitions

- measurement_status: observed, unmeasured, low_coverage, unmappable, invalid, not_applicable.
- resolution_status: resolved, unresolved, invalid.
- Activity evidence: observed or aggregate only. Imported measurements must identify the correct assay/sample and include unit, normalization_id and window_id.
- Contact evidence: observed, aggregate, contact_prior, regularized or fused. Here fused refers solely to specified contact shrinkage between observed contacts and a distance prior.
- Imported contacts require a resolution supplied by the table, matching raw contacts/prior or run configuration; scale alone is insufficient. Optional normalization_id, balancing and window_id must agree when specified.
- Normalization status: complete, partial, zero_support, empty. Complete concerns the planned candidate set.
- samples records identify donor, biological and technical repeats separately. Metadata checks do not replace batch correction.
- Generic features should declare assay, unit, normalization_id and window_id; keep the same definitions during training and prediction.

## Main outputs

| File | Meaning |
|---|---|
| scores.tsv.gz | Every candidate, primary and conditional scores, pace_score_lo/hi, score_scope, A_used, Cbar, support, coverage and reasons |
| region_scores.tsv | Sum of each source/region/gene's unique cells, without another denominator |
| gene_summary.tsv | Candidate counts, coverage and actual denominator |
| resolved_activity.tsv | Qualified experimental activity per unit/assay |
| resolved_contacts.tsv | Contact value, source, prior identity, resolution and policy |
| evidence.tsv / sources.tsv | Traceable input and derived evidence identities |
| multiomics_features.tsv.gz | RNA, histone, CTCF and methylation annotations |
| eta_calibration.json | Actual eta, applicability, fitting and fallback reasons |
| ml_feature_contract.json | Measurement definitions used for the optional classifier |
| qc_report.json / run_manifest.json | QC, context, parameters, input and software hashes |
| resolved_config.yaml / report.md | Reproducible configuration and concise interpretation |

The classifier fields pace_ml_score and pace_ml_probability remain separate from pace_score. NA is preserved when inference is unavailable or out of scope. See [multiomics](MULTIOMICS.md), [parameters](parameters.md) and [formula](FORMULA.md).

## Optional sparse-data settings and output fields

`resolved_activity.tsv` preserves the measured `resolved_value`; `activity_pseudocount`
records an optional offset used in `A_used`. Loading this table again does not apply
the offset twice.

`promoter_weights.tsv` retains all input TSS rows, `pi_original`, effective `pi`,
`tss_selection_reason`, and `n_missing_candidate_contacts`. Score and gene-summary
rows include `tss_retained_weight`, `tss_dropped_ids`, `tss_policy_status`, and
`tss_contact_scope`. `selected_tss_set` refers to a user-selected promoter definition,
not full recovery of the original annotation. `n_tss_used` counts positive effective
weights; `n_tss` counts all original TSSs.

Contacts using a prior also record `prior_source_context` and
`prior_transfer_status`. The run manifest retains the source prior, target context,
and original validation reports separately from validation applicable to the target.
