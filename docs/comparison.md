# Measured-run comparisons and functional benchmarks

PACE preserves each original run and recomputes conditional support shares on the intersection
of measurable units. A simple join-and-subtract of existing normalized scores is not valid when
their denominators differ.

## Compare two runs

```yaml
left: results/animal_a
right: results/animal_b
minimum_common_units: 2
allow_eta_difference: false
allow_evidence_difference: false
```

```bash
pace-livestock compare --config comparison.yaml --out results/comparison
```

The comparison checks catalog, candidate and promoter hashes, target level, context, panel,
units, normalization and formula parameters. At least two common units and positive common
denominators are required. Complete Delta also requires complete planned backgrounds and
compatible structural status. Technical absence leaves complete Delta as NA, even when a
conditional Delta can be reported. Eta comparisons can be explicitly enabled, but then no
complete-background score difference is claimed. B retains its original candidate-gene set in each run.
Contact-policy and measurement-quality comparisons can similarly enable `allow_evidence_difference`; these
produce conditional comparisons, with full Delta withheld because the inference policy differs.

The table contains original/common scores, full/conditional Delta, delta A, delta support,
gene totals, common denominators, comparison IDs and reasons. All deltas are right minus left.
Do not infer enhancer activity or expression direction from the sign of a support-share change.

## Replicate stability

```yaml
replicates: replicate_runs.tsv
minimum_common_units: 2
```

The TSV contains `run_path,donor_id,replicate_type`, where paths resolve relative to the TSV and
replicate type is biological or technical. Every pair is compared on recomputed common
denominators. Pearson correlation is NA for constant or insufficient vectors. Distinct-donor
counts are reported; technical replicates are never counted as independent animals. No
confidence interval is manufactured from a small number of animals.

## Functional-label benchmark

When the run config supplies eta calibration labels or a fixed artifact, the benchmark
adds `PACE_calibrated` alongside the endpoint baselines. It rejects functional evaluation
labels sharing fitting genes, elements or groups. A fixed continuous exponent adds
`PACE_fixed_eta`. Freeze eta before between-animal comparisons; separate refits can
confound evidence changes with parameter changes. See [eta calibration](eta_calibration.md).

```yaml
run_config: examples/measured/config.yaml
labels: functional_labels.tsv
region_membership: region_membership.tsv
stratify: [gene_id]
thresholds: {}
# external_methods:
#   - name: gABC
#     path: gabc_scores.tsv
#     version: recorded_external_version
#     configuration: recorded_external_configuration
```

Decision thresholds are omitted by default. If a threshold has been independently fixed
from appropriate training/calibration data, map its method name to `value`, `source_split`
and `source_id`; omit it when no such evidence exists. The benchmark never selects
a deployment threshold using the test set. Labels use the [standard dictionary](data_dictionary.md)
and preserve every tested positive in the evaluation universe, including absent predictions.
The current implementation trains/evaluates only unambiguous one-to-one region mappings; ambiguous
regions are exported rather than copied to several tiles.

Built-ins are negative distance ranking, ABC-style single-TSS, and PACE eta 0/1. The single TSS
is chosen by smallest coordinate then promoter ID before evaluation, and saved in the report.
All configured formula methods are additionally normalized on their common scoreable sets. External
scores require version/configuration provenance and columns `element_id,gene_id,score`; absent
files become `not_available`. Without external raw support, their common-denominator status is
explicitly not assessed. PACE does not claim to reimplement the external gABC package.

Average precision is the tie-aware sum of recall increments times precision, not trapezoidal PR
area. AUROC is reported only with both classes. Coverage, missing tested positives, candidate
recall and end-to-end recall are separate. Missing scores are not filled with zero. Requested
strata are explicit label columns and include their sample counts. QTL association and input
Hi-C are not substituted for independent functional perturbation labels.
