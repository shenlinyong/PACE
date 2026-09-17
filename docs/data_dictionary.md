# Standard table dictionary

All tables are UTF-8 TSV with a header. `NA` or an empty cell means missing; strings such as
`nan`, `inf`, malformed numbers, duplicate primary keys and unknown foreign keys are rejected.
BED coordinates are zero-based half-open; TSS is a zero-based base coordinate. All arrays use
float64 in the formula kernel. Additional columns are retained; unknown YAML keys are errors.

The following required headers are generated from `pace_livestock.schemas.SCHEMAS`.
Requirements are conditional on the table being used; genome-only does not need fabricated
samples or empty observed-assay files.

## units.tsv

```text
element_id	chrom	start	end	anchor0	element_roles	canonical_catalog_id
```

## region_membership.tsv

```text
region_id	element_id	source_id	membership_rule
```

## promoters.tsv

```text
gene_id	promoter_id	chrom	tss0	strand	pi	pi_source
```

## candidates.tsv

```text
element_id	gene_id	candidate_universe_id
```

## samples.tsv

```text
sample_id	donor_id	assay	biological_replicate	technical_replicate	species	assembly	context_id	source_id
```

## observed_activity.tsv

```text
element_id	sample_id	assay	signal	measurement_status	callable_fraction	unit	normalization_id	window_id
```

## observed_contacts.tsv

```text
element_id	promoter_id	sample_id	contact_value	measurement_status	bin_pair_id	scale	resolution	source_id
```

## resolved_activity.tsv

```text
element_id	assay	observed_value	predicted_value	resolved_value	evidence_id	evidence_type	observation_sample_id	parent_evidence_ids	model_id	calibrator_id	fusion_weight	resolution_status	reason	unit	window_id
```

## resolved_contacts.tsv

```text
element_id	promoter_id	resolved_value	evidence_id	evidence_type	observation_sample_id	prior_id	reliability	resolved_mode	bin_pair_id	resolution_status	reason	scale
```

## predictions.tsv

```text
element_id	assay	predicted_value	model_id	unit	normalization_id	window_id	status
```

## features.tsv

```text
entity_type	entity_id	feature_name	value	evidence_id	status
```

## methylation.tsv

```text
chrom	dyad_start0	methylated_count	total_count	sample_id	assay
```

## expression.tsv

```text
gene_id	sample_id	tpm	status
```

## labels.tsv

```text
label_id	assayed_region_id	gene_id	context_id	perturbation_type	effect_direction	effect_size	label_status	assay_id	group_id	source_id
```

## evidence.tsv

```text
evidence_id	evidence_type	source_id	parent_evidence_ids	model_id	unit	processing_method	checksum
```

## sources.tsv

```text
source_id	path_or_accession	source_type	assembly	processing_method	normalization_id	checksum
```

## Values, identities and provenance

- `measurement_status`: observed, unmeasured, low_coverage, unmappable, invalid, not_applicable.
- `resolution_status`: resolved, unresolved, invalid. Physical structure is reported separately.
- `evidence_type`: observed, sequence_prediction, contact_prior, fused, aggregate.
- `normalization_status`: complete, partial, zero_support, empty. Complete concerns the planned
  candidate set only, not discovery of all biological enhancers.
- `target_level`: individual or population_mean; main `estimand` is bulk_proxy.
- Canonical `window_id`: `grid:<width>:mean`, e.g. `grid:500:mean`. Imported resolved quantitative
  activity must additionally include `normalization_id`; it is checked against source/model units.
- `labels.label_status`: enhancing_positive, powered_negative, or an explicit excluded status
  such as low_power or not_tested. Enhancing positives require effect_direction=down.
- `features.entity_type`: element, promoter, gene, edge. Edge `entity_id` uses `element_id|gene_id`;
  include element_id/gene_id columns when exporting structured E–G features. Use distinct feature
  names for different biological levels. Status determines whether the value is usable.

Observed rows reference actual samples and their assay/context. Predictions reference model
assets with no experimental sample. Derived evidence carries parent IDs, method and row hash;
run outputs include an expanded evidence.tsv and sources.tsv so parent IDs remain resolvable.
Input-file and model-manifest SHA256 hashes are in run_manifest.json. Actual catalogs/candidates/
promoters are hashed from contents, independent of caller-supplied human-readable IDs.
`software_sha256` hashes the installed Python implementation independently of installation path,
so two source revisions sharing a package version remain distinguishable.

## Output-only fields

`eta_calibration.json` records the numeric exponent, status, applicability scope and
calibration provenance. `run_manifest.json` repeats it under `allocation`, and the
comparison contract stores the numeric eta. Functional eta calibration uses a dedicated
element-level table; see [its schema](eta_calibration.md#required-functional-label-table).

`scores.tsv.gz` retains all candidates with A_used, Cbar, B, support, log_support, denominator,
log_denominator, pace_score, scoreable, support_status, normalization_status, actual normalization
ID, distance_bp, n_tss, sources, structural_status and reasons. Raw overflow is NA with a finite
log_support and status=overflow; true zero is support=0 and log_support=NA in TSV (mathematically
negative infinity). The reader recovers that zero without conflating it with missing evidence.

`pace_ml_score` and `pace_ml_probability` are separate optional fields, NA unless computed with
an applicable model. An uncalibrated classifier only fills its score. Gene totals, coverage and
actual denominator identities are reported in gene_summary.tsv. Comparison fields are documented
in [comparison.md](comparison.md); all deltas are right minus left.
