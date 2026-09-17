# Inputs and migration

## Quantified-table CLI

Use `scripts/pace.py --pairs table.tsv --output scores.tsv`. Each row is a unique enhancer/gene/TSS. Process a single sample/assembly per invocation, or provide distinct `sample_id` values when using the core directly. Do not mix assemblies within a sample ID.

Required columns are `chr`, `start`, `end`, `TargetGeneEnsemblID`, `TargetGeneTSS`, `activity`, `contact_prior`. The activity can instead be computed using `--activity-config` and named signal columns. Coordinates are BED0. NA is a missing measurement; literal 0 means measured zero. All scores and quality fields reject invalid negative/infinite values.

Optional fields

| Column | Meaning |
|---|---|
| sample_id | Normalization scope, one sample and assembly |
| contact_observed | Actual contact H before prior replacement |
| contact_expected | Compatible distance expectation D, positive |
| contact_reliability | Local QC weight in [0,1], NaN if unknown |
| contact_source | matched / surrogate / distance_prior / unknown |
| tss_weight | Promoter-use weight, normalized over all unique gene TSSs before window filtering |
| activity_quality | Independently supplied activity measurement QC |
| tss_quality | Independently supplied TSS annotation QC |
| catalogue_quality | Independent candidate catalogue coverage QC |
| unassigned_mass | Independently estimated residual support in raw_support units, constant per gene |

No quality field is inferred to be 1 because an assay file exists. Supplying QC is a modelling decision that must be documented and validated. Missing quality does not suppress structural scoring, but it prevents a sufficient-input-evidence designation. Expected contact is not fitted by this CLI; label-derived or selected-edge averages are not suitable substitutes.

Activity configuration example

```json
{"signals": {
  "ATAC": {"weight": 1, "scale": 1, "quality_column": "ATAC_quality"},
  "H3K27ac": {"weight": 1, "scale": 1, "quality_column": "H3K27ac_quality"}
}}
```

The scale of 1 is suitable only for already normalized example values. A missing signal column becomes all-NA and remains part of the planned-modality QC denominator. A configured quality column absent from the table also becomes all-NA; its quality is not inferred from signal presence. If a modality is absent by design, declare the planned input regime explicitly. Changing that regime changes the interpretation of activity quality.

## Raw file adapters

The standalone signal calculator and workflow neighborhood adapter share activity aggregation and the prediction kernel. Their input quantification remains subject to assay-specific preprocessing. Use pre-normalized bigWig tracks for comparable units. Raw BAM/BED counts do not become cross-assay comparable merely because both enter the same aggregation function. The legacy normalization name `rpkm` denotes candidate-total scaling, not true library RPKM.

Both raw-read adapters use read counts per kb of candidate width and exact bigWig summaries. Candidate files may be BED3, BED4, BED6 or wider; coordinates retain their original columns and duplicate display names do not merge regions. The neighborhood adapter's former candidate-total normalization was corrected on 16 September 2026. Explicitly requested files must exist; omit the corresponding argument when an assay is unavailable.

`workflow/scripts/pace_predict.py` accepts `--contact_metadata qc.tsv`, keyed by chr/start/end/TargetGeneEnsemblID/TargetGeneTSS. `scripts/calculate_pace_score.py` accepts the corresponding per-sample `contact_metadata` file path. File-derived H is retained if not overridden. A compatible D, reliability and source are still required before measured contact affects the score.

Use `scripts/prepare_tss.py --gtf annotation.gtf.gz --output tss.tsv`, then provide this headered TSV as genes input. The historical gene-only reference builder is not a transcript-TSS discovery pipeline.

RNA annotation requires `gene_id` and `TPM`. No automatic symbol fallback or hard expression filtering is used in prediction. A missing RNA identifier remains unknown. Deprecated expression-weight options log that no multiplication is performed.

Outputs are gene-level; `TargetGeneTSSs` lists the contributing input alternatives and `contact_gene` is their weighted contact. `TargetGeneTSS`, `contact` and `distance` inherited from one representative record are diagnostics only, not a promoter-specific biological assignment. Full alternative weights and contacts remain in the input manifest; save it with the prediction file.

## Independent comparison and eQTL annotation

```bash
python scripts/benchmark_compare.py --labels tested_pairs.tsv \
  --abc abc.tsv --pace pace.tsv --abc-threshold 0.02 \
  --pace-threshold 0.02 --output benchmark
```

The thresholds above illustrate syntax only. Select each on an independent validation set, freeze it and document achieved test precision. `label` must be 0/1 or NA. Both classes must reflect actual experiments, not unsupported candidates. Supply unfiltered predictions so omitted candidates can be counted. Only ABC and PACE are supported as external method inputs.

```bash
python scripts/eqtl_pair_support.py --predictions pace.tsv \
  --eqtls independently_filtered_eqtls.tsv --output supported.tsv
```

The eQTL table requires chr, **pos0**, and gene_id. Convert ordinary 1-based VCF/eQTL positions to BED0 before use. Prefilter associations independently, use matching tissue/assembly, and account for LD in downstream inference. This command annotates association support only; it does not estimate precision or FDR from absent associations.
