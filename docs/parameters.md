# Canonical run parameters

The defaults below are generated from `pace_livestock.config.DEFAULTS`. Input paths resolve
relative to the YAML file; `null` means no asset/input. Required context and data must be supplied.
Do not use synthetic model assets in research/validated profiles.

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
  predictions: null
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
promoters:
  weights: provided
allocation:
  eta: 0
  missing_policy: fixed_gene_set
sequence:
  model_path: null
  max_n_fraction: 0.05
fusion:
  calibrator_path: null
  quality_stratum: default
genome:
  individual_id: null
  reference_path: null
  variant_path: null
  callability_path: null
  ploidy_path: null
  sample_id: null
  phase_policy: require_phase_or_single_variant_scenario
  unrecorded_site_policy: require_callable
  sv_assessed: false
multiomics:
  mode: annotate
  model_path: null
comparison:
  full_delta_requires_complete: true
  allow_conditional_intersection: true
  minimum_common_units: 2
output:
  format: tsv_gz
  retain_all_candidates: true
seed: 17
```

## Scientific decisions

| Parameter | Meaning and boundary |
|---|---|
| eta | 0 or 1 only; 0 is the default and skips B completely |
| panel | ATAC, DNase, H3K27ac, ATAC+H3K27ac or DNase+H3K27ac; no per-edge fallback |
| catalog profile | canonical_grid requires aligned equal-width nonoverlapping units; provided_regions is measured-only |
| include_promoter_units | checks that trusted TSS cells occur in the catalog; the small hand-calculation example explicitly disables it |
| promoters weights | provided pi must sum to one per complete gene; equal uses deduplicated physical TSSs |
| contact mode | observed, prior_only, or explicit linear shrinkage; no hidden quality heuristic |
| near_diagonal_bp | declare an additional near-distance region; observed same-bin pairs always use the explicit near-diagonal policy, even when this value is 0 |
| near_diagonal_policy | use a same-scale prior or leave unresolved; no invented diagonal correction |
| reliability | explicit [0,1] and reliability_source for shrinkage; endpoint weights do not require the unused source |
| minimum_callable_fraction | protocol-specific observation filter; it is not a universal quality score |
| replicate_aggregation | normalized technical means within biological replicate, then equal donor mean; counts must be pooled upstream |
| genome unrecorded policy | require callable evidence by default; assume_reference must be an explicit research assumption |
| sequence max_n_fraction | recorded model-input QC threshold; output is unavailable above it |
| comparison minimum | at least two common units and positive denominators; no full Delta through technical missingness |

The main score uses no arbitrary epsilon, residual unknown support, RNA multiplier or auxiliary
mark multiplier. Weights, units, windows and model scope are checked before evidence resolution.
All numerical and scientific defaults are engineering starting points, not species-specific
optimal settings. The supplied [model contract](model.md) explains their assumptions.
