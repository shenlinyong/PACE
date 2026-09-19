# Canonical run parameters

The defaults below are generated from `pace_livestock.config.DEFAULTS`. Input paths resolve
relative to the YAML file; `null` means no asset/input. Required context and data must be supplied.
Do not use synthetic model assets in research/validated profiles.

## Choose the small set of parameters first

Most users only need to decide four things before reading the full defaults: `regime` (`measured`, `hybrid` or `genome_only`), the fixed activity `panel`, the contact `mode`, and the output target level. The remaining fields describe files, provenance and optional calibration.

| First decision | Typical choice | What it controls |
|---|---|---|
| Evidence mode | `measured`, `hybrid`, `genome_only` | Which sources may resolve activity and contact |
| Activity panel | `ATAC`, `DNase`, `H3K27ac`, or one supported two-layer panel | The layers used by the geometric-mean activity formula |
| Contact mode | `observed`, `shrinkage`, `prior_only` | Whether contacts come from measurements, an explicit mixture, or a prior |
| Replicates | Any positive number | Sample rows and aggregation; there is no fixed requirement of three |

Additional RNA-seq, histone, CTCF and methylation inputs are declared in their own tables and keep their own `sample_id` and `assay`. They are annotations or separately validated ML features by default, not unlabelled extra score multipliers.

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

## Scientific decisions

| Parameter | Meaning and boundary |
|---|---|
| eta | `auto` (default) or a finite number in [0,1]; auto falls back to 0 without suitable functional calibration; 0 skips B completely |
| allocation.labels_path / calibrator_path | Mutually exclusive functional calibration TSV or frozen eta JSON; see [calibration](eta_calibration.md) |
| allocation.minimum_genes | At least 2; default 3 informative genes, including each CV training subset; an engineering guard, not a sample-size calculation |
| allocation.minimum_groups | Default 3 independent connected components; overlapping group/gene/element identities connect rows |
| allocation.validation_folds | Default 5; fold count is limited by independent components and information requirements |
| allocation.minimum_positive_fraction | Default 0.8; required fraction of held-out groups showing positive AP improvement |
| panel | ATAC, DNase, H3K27ac, ATAC+H3K27ac or DNase+H3K27ac; no per-edge fallback |
| catalog profile | canonical_grid requires aligned equal-width nonoverlapping units; provided_regions is measured-only |
| include_promoter_units | Require eligible trusted TSS cells in the catalog; complete reference-boundary cells are checked using chrom_sizes_path |
| catalog.chrom_sizes_path | Exported chromosome lengths, required to justify exclusion of incomplete terminal TSS cells |
| catalog.candidate_radius_bp | Declared candidate search radius, default 5,000,000 bp; match preparation, not a universal optimal distance |
| promoters weights | provided pi must sum to one per complete gene; equal uses deduplicated physical TSSs |
| contact.resolution / normalization_id / balancing / window_id | Optional declarations strengthen the measurement contract; actual observations/prior must agree. Null means unspecified, not interchangeable |
| contact mode | observed, prior_only, or explicit linear shrinkage; no hidden quality heuristic |
| near_diagonal_bp | declare an additional near-distance region; observed same-bin pairs always use the explicit near-diagonal policy, even when this value is 0 |
| near_diagonal_policy | use a same-scale prior or leave unresolved; no invented diagonal correction |
| reliability | explicit [0,1] and reliability_source for shrinkage; endpoint weights do not require the unused source |
| minimum_callable_fraction | protocol-specific observation filter; it is not a universal quality score |
| replicate_aggregation | normalized technical means within biological replicate, then equal donor mean; counts must be pooled upstream |
| genome unrecorded policy | require callable evidence by default; assume_reference must be an explicit research assumption |
| sequence max_n_fraction | recorded model-input QC threshold; output is unavailable above it |
| methylation.minimum_coverage | Minimum reads at a CpG, default 1; protocol-specific and frozen for comparisons/ML |
| methylation.promoter_upstream_bp / promoter_downstream_bp | Strand-aware promoter windows, defaults 2000/500 bp, include TSS base |
| methylation.reference_cpg_path | Optional reference counts by entity_type/entity_id/n_cpg; missing denominator keeps coverage NA |
| comparison minimum | at least two common units and positive denominators; no full Delta through technical missingness |

The main score uses no arbitrary epsilon, residual unknown support, RNA multiplier or auxiliary
mark multiplier. Weights, units, windows and model scope are checked before evidence resolution.
All numerical and scientific defaults are engineering starting points, not species-specific
optimal settings. The supplied [model contract](model.md) explains their assumptions.

## Required and optional inputs by mode

The block above is a defaults reference, not a runnable project: context and input
paths are intentionally null. [Three complete real-project templates](USER_GUIDE.zh-CN.md)
and [runnable fixture configurations](TUTORIAL.md) show how to fill them.

| Section | measured | hybrid | genome_only |
|---|---|---|---|
| Context, units, promoters, candidates | Required | Required | Required |
| Samples and sources | Required for observed data | Required for observed data | Only when importing actual annotations |
| Observed activity | Required or valid observed resolved evidence | Optional qualified observations | Forbidden |
| Observed contact | Required unless explicit prior policy | Optional under declared prior policy | Forbidden |
| Sequence model | Not required | Required when using predictions | Required |
| Fusion calibrator | Not used for primary activity | Required for actual source fusion | Not used |
| Contact prior | Optional explicit fallback/shrinkage | Optional/required according to contact mode | Required; prior_only |
| Genome reference | Not required for measured core | Required for local sequence prediction | Required for local sequence prediction |
| VCF, callable regions, ploidy | Not a substitute for a genomic model | Required for individual reconstruction | Required for individual reconstruction |
| RNA, auxiliary marks, methylation | Optional annotation/ML | Optional annotation/ML | Optional available annotation; measured activity/contact requires hybrid |

External quantitative tables require a matching manifest. Genomic inputs, when
supplied, remain binding constraints even with imported predictions. Validated
profile additionally checks applicable validation evidence; it cannot create it.
