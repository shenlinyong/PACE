# Practical experimental-data workflow

[中文完整手册](USER_GUIDE.zh-CN.md) · [Formula](FORMULA.md) · [CLI](cli.md)

These commands use your experimental files. Replace the example paths and biological context before running. They do not generate missing experiments.

## Initialize a project

```bash
pace init --species chicken --assembly GRCg7w --tissue liver \
  --panel H3K27ac --out chicken_liver
```

This creates a short `config.yaml`, empty input tables with correct headers, and a preparation checklist. `--catalog-dir prepared/catalog` imports existing tables and the exported catalog definition. Fill the missing tables before `pace validate`. Context identifiers must use the same spelling in all inputs. `--target-level population_mean` is required when explicitly pooling several donors.

The canonical catalog still comes from `pace prepare --config prepare_catalog.yaml --out prepared/catalog`, as documented in [input preparation](input_preparation.md). It includes promoter units and a mapping from original regions to grid cells. Requantify the resulting cells; do not assign a whole peak's read count to every overlapping cell.

## Generate Hi-C query pairs

```bash
pace prepare-pairs --catalog-dir prepared/catalog --out prepared/pairs
```

Use `prepared/pairs/pairs.tsv` as `pairs:` in the cooler preparation YAML. No hand-written table-joining script is needed. A physical element–promoter pair appears once even when shared by gene annotations.

## Merge measured replicates

```bash
pace merge-tables --table observed_activity \
  --inputs prepared/animal1_atac/observed_activity.tsv \
           prepared/animal1_h3k27ac/observed_activity.tsv \
  --out prepared/activity
```

The command rejects repeated element/sample/assay keys. The same command supports `observed_contacts`, `expression`, and `methylation`. Merging rows does not average donors; the run uses the specified sample hierarchy. Register sources and samples before scoring.

## Normalize raw counts

If your window table contains **raw filtered fragment counts**, prepare these TSV headers:

```text
# counts.tsv
element_id sample_id assay count
# library_sizes.tsv
sample_id library_size
```

Use actual tabs, not the spaces shown above. One count per element/sample/assay; one full-library total per sample. Missing counts are NA; measured zero counts are 0.

```bash
pace normalize-activity --counts counts.tsv --library-sizes library_sizes.tsv \
  --units prepared/catalog/units.tsv --out prepared/normalized
```

Output signal = count × 10^6 / library_size / window_length, in CPM_per_bp. The library total counts all retained fragments, not only peak-overlapping fragments. Paired-end data must use the same fragment-counting convention in numerator and denominator. Variable-width regions require `--window-id`. The default identifier is `CPM_density_v1`; a project can supply a more specific `--normalization-id`.

This command does not infer total library size from a bigWig, renormalize TPM, or normalize fold-enrichment tracks. Use the bigWig adapter for already normalized quantitative tracks and name their actual processing protocol. Changing a protocol string does not make incompatible assays comparable. No default cross-animal quantile normalization is imposed: it could remove real global biological differences.

## Fit a matching contact prior from cool/mcool

```bash
pace fit-prior --cooler data/animal1.mcool::/resolutions/5000 \
  --species chicken --assembly GRCg7w --tissue liver \
  --scale balanced_contact_protocol_1 --normalization-id hic_norm_protocol_1 \
  --test-chromosomes chr2 --out models/animal1_contact
```

Use chromosome names actually present in the file. `--test-chromosomes` is optional; omitting it reports that no holdout was performed. It is not correct to rename a liver fit as a brain-validated asset. The defaults fit distances from 5 kb to 5 Mb, in 30 logarithmic bins; resolution must not exceed the chosen minimum fit distance. For a 10 kb map, supply `--min-distance 10000`. `--no-balanced` uses raw counts and requires matching raw-count preparation. The command only accepts fixed-bin, symmetric-upper coolers.

The output includes `manifest.json`, `distance_bins.tsv`, and `fit_report.json`. All callable bin opportunities enter distance means, including zeros absent from sparse storage; invalid balancing weights are excluded. Fitting is streamed over sparse pixels. The optional chromosome holdout assesses decay fit, not enhancer function or generalization across tissues.

Use the same scale/normalization/balancing/resolution when preparing contact observations. For observations extracted from this map:

```yaml
contact:
  mode: observed
  scale: balanced_contact_protocol_1
  prior_path: models/animal1_contact
  near_diagonal_policy: prior_or_neighbor
  allow_prior_fallback: true
  pseudocount: auto
```

The matching prior now supports near-diagonal replacement, sparse-count regularization and explicitly requested missing-contact fallback. Rows retain their individual evidence sources. Multiple samples require genuinely comparable scales before pooling; fitting one sample's arbitrary raw-count scale does not calibrate another sample.

`pace fit-contact-prior --config fit.yaml --out models/prior` remains available for explicitly prepared distance/contact tables. `pace fit-prior --config fit.yaml ...` is its short alias.

## No Hi-C: an explicit baseline

```bash
pace init --species chicken --assembly GRCg7w --tissue liver \
  --panel H3K27ac --prior-preset abc_human --out chicken_distance_baseline
```

You must still fill the experimental activity and annotation tables. This chooses a human-derived contact **shape** on a relative scale, labeled `abc_human_default` and `unvalidated_for_target_context` in `run_manifest.json` under `asset_manifests.contact_prior`, alongside its source species and numerical parameters. It is not a livestock optimum and cannot be mixed with experimental contact tables or used with the validated profile. Prefer a suitable same-species fit where possible and compare both as prespecified sensitivity analyses.

## Promoter weights from a measured signal

Prepare `promoter_signal.tsv` with `promoter_id` and a nonnegative `signal` measured over the same specified promoter window. Supply one already aggregated value per physical promoter.

```bash
pace prepare-promoter-weights --promoters prepared/catalog/promoters.tsv \
  --signals promoter_signal.tsv --assay H3K4me3 \
  --normalization-id H3K4me3_CPM_TSS500 --out prepared/promoter_weights
```

Use the resulting promoters.tsv in the run. Weights are fixed before enhancer scoring. An unavailable TSS signal errors. All-zero genes error unless `--zero-policy equal` is explicitly chosen. ATAC, DNase and CAGE are also supported; these measurements are promoter-activity proxies and require their own biological validation.

## Score and inspect coverage

```bash
pace validate --config experiment.yaml
pace run --config experiment.yaml --out results/experiment
```

| Output | Read it as |
|---|---|
| pace_score | Full planned-background share; NA for incomplete genes by default |
| pace_score_conditional | Share within the available subset only |
| pace_score_lo / pace_score_hi | Sensitivity range under specified support bounds, not a confidence interval |
| region_scores.tsv | Sum of unique cells for each original region/gene, without another denominator |
| resolved_contacts.tsv | Raw, prior and pseudocount values, reasons and sources |
| qc_report.json | Coverage, unresolved near-diagonal contacts and actual regularization counts |

An optional support-bounds TSV has `element_id,gene_id,support_lower,support_upper,bound_source`. Bounds apply to the **final unnormalized support**, including the chosen allocation exponent, and only to unresolved edges. An NA upper bound is unbounded. Do not estimate bounds from held-out functional outcomes or reuse an old bound after changing activity/contact scales without reassessment.

Default outputs must be new directories. `--force` permits replacing a recognized PACE result only after preserving it as a sibling `NAME.backup-*` directory. A failed computation leaves the previous result intact. Arbitrary input directories and symlinks cannot be overwritten with this flag.

## Optional settings for sparse activity and alternative TSSs

中文逐项说明见[活性](FORMULA.zh-CN.md#二活性)和[多-TSS](FORMULA.zh-CN.md#四多-tss)。

These options are off by default. Add only the settings needed for the analysis:

```yaml
activity:
  panel: [ATAC, H3K27ac]
  pseudocounts:
    H3K27ac: 0.1  # example only; use the units of your normalized signal
promoters:
  weights: provided
  minimum_weight: 0.01
  missing_policy: drop_missing
  minimum_retained_weight: 0.9
```

The pseudocount is added after replicate aggregation; NA remains NA. Compare
results with zero offsets because positive offsets can increase background support.
There is no universally appropriate 0.1 offset across assays or normalizations.

TSS filtering is performed once per gene. Every candidate for that gene uses the
same retained set; a TSS missing any planned contact is removed globally under
`drop_missing`. At least 90% of the original weight must remain in this example.
The original candidate catalog is preserved. Inspect `promoter_weights.tsv`,
`tss_contact_scope` and `tss_retained_weight` before interpreting scores.
Different selected promoter sets do not define the same biological comparison.

## Using a prior fitted in another tissue

Fit a prior using the source tissue name. For 25 kb matrices the fitting lower
bound defaults to 25 kb:

```bash
pace fit-prior --cooler reference_liver.mcool::/resolutions/25000 \
  --species chicken --assembly GRCg7w --tissue liver \
  --out priors/liver
```

In the target muscle configuration, add:

```yaml
context:
  species: chicken
  assembly: GRCg7w
  context_id: muscle
contact:
  mode: prior_only
  prior_path: priors/liver
  scale: cooler_native
  allow_cross_context_prior: true
```

This block supplements the existing activity, sample and catalog inputs. Keep the
prior's source tissue name intact. The output records source and target contexts
and marks the transfer unvalidated; it cannot use the `validated` profile. Measured
muscle activity remains necessary. If combining a prior with measured target
contacts, first establish matching normalization, scale, resolution and windows;
a source-tissue amplitude is not automatically calibrated to a target library.
