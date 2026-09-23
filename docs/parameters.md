# Run parameters

The only activity regime is measured. Choose a fixed assay panel, contact policy and target level. Default contact uses observations; a supplied prior or shrinkage policy must be explicitly selected. Other experimental omics are annotations or separate ML features.

The following defaults are generated from `pace_livestock.config.DEFAULTS`. Supply biological context and actual file paths. YAML paths are relative to the configuration file; null means not provided.

```yaml
schema_version: pace-1
run_id: pace
regime: measured
execution_profile: research
estimand: bulk_proxy
target_level: individual
context:
  species: null
  assembly: null
  context_id: null
inputs:
  units: null
  region_membership: null
  promoters: null
  candidates: null
  samples: null
  observed_activity: null
  observed_contacts: null
  resolved_activity: null
  resolved_contacts: null
  features: null
  methylation: null
  expression: null
  labels: null
  evidence: null
  sources: null
catalog:
  profile: canonical_grid
  width_bp: 500
  offset_bp: 0
  include_promoter_units: true
  chrom_sizes_path: null
  candidate_radius_bp: 5000000
activity:
  panel:
  - ATAC
  - H3K27ac
  combine: geometric_equal
  missing_policy: unresolved
  minimum_callable_fraction: 0.0
  replicate_aggregation: equal_donor_mean
contact:
  mode: observed
  scale: depth_normalized_contact
  prior_path: null
  near_diagonal_policy: prior_or_unresolved
  near_diagonal_bp: 0
  allow_prior_fallback: false
  reliability: null
  reliability_source: null
  resolution: null
  normalization_id: null
  balancing: null
  window_id: null
promoters:
  weights: provided
allocation:
  eta: auto
  missing_policy: fixed_gene_set
  labels_path: null
  calibrator_path: null
  minimum_genes: 3
  minimum_groups: 3
  validation_folds: 5
  minimum_positive_fraction: 0.8
multiomics:
  mode: annotate
  model_path: null
methylation:
  minimum_coverage: 1
  promoter_upstream_bp: 2000
  promoter_downstream_bp: 500
  reference_cpg_path: null
comparison:
  full_delta_requires_complete: true
  allow_conditional_intersection: true
  minimum_common_units: 2
output:
  format: tsv_gz
  retain_all_candidates: true
seed: 17
```

## Decisions that affect interpretation

| Setting | Interpretation |
|---|---|
| target_level | individual requires one donor; population_mean aggregates donors equally |
| activity.panel | Fixed measured assay set; no automatic per-element deletion of missing layers |
| minimum_callable_fraction | Minimum usable track coverage, not a multiplier of activity |
| catalog.profile | canonical_grid uses the declared width and offset; provided_regions uses supplied non-overlapping regions |
| include_promoter_units | Include eligible promoter units in the planned denominator; do not drop them silently |
| contact.mode | observed, prior_only or shrinkage; activity remains measured in every case |
| contact.reliability | Explicit observed-contact weight for shrinkage; needs reliability_source |
| near_diagonal_policy | Use a supplied compatible prior or leave the contact unresolved |
| promoters.weights | Provided fixed weights or explicit equal weights over distinct physical TSSs |
| allocation.eta | auto defaults to zero without eligible validation; fixed [0,1] values are explicit sensitivity settings |
| allocation labels/model | Use functional labels or a frozen calibrator, not both |
| multiomics.mode | annotate by default; ml requires an applicable independently trained model |
| comparison | Recompute common denominators; partial comparisons remain conditional |

No minimum sample-count setting by itself establishes adequate statistical power. Preserve the actual learned parameters, scientific context and source commit with the analysis. See [CLI](cli.md), [labels](eta_calibration.md) and [inputs](data_dictionary.md).
