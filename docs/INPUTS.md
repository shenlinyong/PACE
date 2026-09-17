# Input specification

PACE (Prediction of Activity-based regulatory Connections for Enhancers) needs a consistent species, assembly, candidate background and gene-ID namespace. Run one sample/assembly at a time unless explicitly supplying distinct `sample_id` values to the numerical-table interface.

## Coordinate and missing-value rules

- Enhancer intervals are zero-based, half-open: `start` is included and `end` is excluded.
- TSS coordinates are zero-based single positions. `prepare_tss.py` converts GTF starts/ends and strand correctly.
- Chromosome strings must match exactly across inputs. `1` and `chr1` are different identifiers.
- Use stable gene IDs for joins and normalization. Gene symbols are labels and need not be unique.
- Use blank/`NA` for unavailable measurements; literal `0` means measured zero. Negative or infinite scores/quality inputs are invalid.
- Do not restrict candidates to experimentally tested or high-scoring links before computing the score.

## Quantified candidate table

Entry point: `scripts/pace.py --pairs pairs.tsv --output scores.tsv`.

Each row represents **one enhancer–gene–TSS combination**. Use a tab-separated file with a header. The included [fixture](../example_quantified/candidates.tsv) is a complete example.

| Column | Required? | Definition |
| --- | --- | --- |
| `chr` | Yes | Chromosome or contig |
| `start`, `end` | Yes | Integer enhancer coordinates with `end > start` |
| `TargetGeneEnsemblID` | Yes | Stable target ID; the column name does not require the ID to come from Ensembl |
| `TargetGeneTSS` | Yes | Integer BED0 TSS |
| `contact_prior` | Yes | Finite positive prior for this enhancer–TSS pair |
| `activity` | Yes, unless using activity JSON | Nonnegative enhancer-level activity; may be `NA` |
| Named signal columns | With activity JSON | Nonnegative preprocessed assay values; `NA` allowed |
| `sample_id` | Optional | Separate sample/assembly normalization scope; defaults to `single_sample` |
| `tss_weight` | Optional | Prior TSS-use weight; constant for a gene/TSS across enhancers |
| `contact_observed`, `contact_expected`, `contact_reliability`, `contact_source` | Optional | Contact evidence; see below |
| `activity_quality`, `tss_quality`, `catalogue_quality` | Optional | Independent QC values in [0,1]; unknown by default |
| `unassigned_mass` | Optional | Nonnegative residual support, constant per gene, in raw-support units |

Activity and activity quality must be consistent for the same enhancer within a sample. TSS quality is gene/TSS-specific; catalogue quality and residual support are gene-specific. Conflicting duplicate rows are rejected. Exact duplicate rows do not create additional evidence.

The table CLI does not generate pairs or apply a distance cutoff. Assign TSS weights using the full available gene catalogue **before** any window-based row selection. If weights are omitted, the kernel uses uniform weights over TSSs visible in the supplied table; it cannot infer omitted alternatives.

## Activity JSON

```json
{
  "signals": {
    "ATAC": {"weight": 1.0, "scale": 1.0, "quality_column": "ATAC_quality"},
    "H3K27ac": {"weight": 1.0, "scale": 1.0, "quality_column": "H3K27ac_quality"}
  }
}
```

Pass this with `--activity-config activity.json`. It recomputes activity and its quality fields from the named signals, replacing any pre-existing activity fields. `scale: 1` is appropriate for the normalized synthetic example; it is not an instruction to mix arbitrary raw count units. Use fixed, positive, independently chosen scales for a biological run.

A named assay column absent from the table becomes missing. A named quality column absent from the table also remains unknown. A missing planned assay contributes no measurement quality; removing the assay from the JSON instead changes the planned input regime. Weights are design priors, not learned feature importances.

## Candidate BED and signal files

The direct adapters accept BED3, BED4, BED6 or wider candidate files without a header. Interval coordinates identify candidates; duplicate display names do not merge different intervals.

| Input | Format | Handling |
| --- | --- | --- |
| Accessibility | BAM, BED/tagAlign (including gzip), or bigWig | BAM/BED reads per kb via bedtools; exact bigWig regional summaries |
| H3K27ac | Same signal formats | Optional activating layer |
| Peaks for candidate construction | Ten-column narrowPeak | Uses summit offset and `signalValue`; missing/invalid summit definitions need upstream correction |
| Chromosome sizes | Two-column, headerless TSV | Exact contig names and lengths |
| Transcript annotation | GTF or `.gtf.gz` | Convert with `scripts/prepare_tss.py` |
| TSS catalogue | Headered TSV from `prepare_tss.py` | Preferred gene input for all distinct annotated TSSs |

Missing bigWig intervals remain missing. A zero read count in an available read file is an observed count of zero. Neither case alone establishes biological inactivity. Raw-file adapters do not automatically normalize library depths or provide local assay QC. Use the table route for controlled scaling and planned missing-modality QC.

## Direct sample sheet and YAML

The standalone signal calculator needs a `biosample` column and file columns for the assays used. A runnable example is [config/biosamples_direct.example.tsv](../config/biosamples_direct.example.tsv), selected by [config/config_direct.example.yaml](../config/config_direct.example.yaml).

```text
biosample    ATAC                H3K27ac
Pig_Liver    data/liver_ATAC.bw  data/liver_H3K27ac.bw
```

The display above uses spaces for readability; the actual file must use **tabs**. Optional columns include `DHS`, `RNA_seq`, `HiC_file`, `HiC_type`, `HiC_resolution` and `contact_metadata`. Leave unavailable files blank and ensure the chosen sample occurs exactly once. Paths are resolved from the working directory; run commands from the repository root or use absolute paths.

The Snakemake sample table has a fuller required header. Copy [config/config_biosamples.tsv](../config/config_biosamples.tsv) and follow the [workflow instructions](TUTORIAL.md#optional-snakemake-workflow); the minimal direct sample sheet is not a Snakemake replacement.

## Contact metadata

The direct predictor accepts `--contact_metadata contact_metadata.tsv`. The standalone calculator reads the path from the sample-sheet `contact_metadata` column. The join key is:

```text
chr  start  end  TargetGeneEnsemblID  TargetGeneTSS
```

Add these columns:

| Field | Meaning |
| --- | --- |
| `contact_observed` | Optional override of the file-derived measurement $H$ |
| `contact_expected` | Positive expected value $D$ with matching sample, normalization, resolution and distance |
| `contact_reliability` | Independently derived local measurement reliability in [0,1] |
| `contact_source` | `matched`, `surrogate`, `distance_prior` or `unknown` |

A qualified observation requires finite $H$, positive $D$, known reliability and a matched/surrogate source. Otherwise the score uses the prior and records the limitation. Do not estimate $D$ only from selected positive links, or define reliability directly from the edge's contact count. Surrogate tissue/species/assembly correspondence must be resolved outside PACE and recorded.

BEDPE contact input is headerless with seven columns:

```text
chr1  start1  end1  chr2  start2  end2  score
```

Use tabs and matching coordinate conventions. An explicit BEDPE zero is measured zero; an omitted edge remains missing. Binary `.hic`/`.cool` input requires the [optional readers](INSTALLATION.md#optional-contact-readers).

## RNA, eQTL and experimental labels

RNA context uses `gene_id` and `TPM` with stable IDs matching the TSS catalogue. No symbol-based fallback is performed by the primary predictor. Missing IDs remain unknown; RNA does not multiply structural scores.

`eqtl_pair_support.py` requires `chr`, zero-based `pos0` and `gene_id`. Convert conventional 1-based VCF/eQTL coordinates before use. Independent association filtering and tissue/assembly matching are prerequisites.

`benchmark_compare.py` joins on `chr`, `start`, `end`, `TargetGeneEnsemblID` and requires an experimental `label` (0, 1 or missing). Use the appropriate score column names shown by `--help`, unfiltered predictions and independently frozen thresholds. Lack of an association or lack of a prediction is not automatically a negative experimental label.
