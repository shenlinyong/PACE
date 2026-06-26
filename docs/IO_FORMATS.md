# PACE Input/Output Format Specification

A precise reference for every file PACE reads and writes. All tabular files
are **tab-separated**. Coordinates are 0-based, half-open (BED convention)
unless noted.

---

## Inputs

### Accessibility — ATAC-seq / DNase-seq (required)
- **tagAlign** (`.tagAlign[.gz]`): BED6 of reads
  `chrom  start  end  name  score  strand`.
- **BAM** (`.bam`): coordinate-sorted, indexed.
- **bigWig** (`.bw`): coverage track (read via pyBigWig).

### Histone ChIP-seq (optional)
`.bam` or `.bw`, one file per mark. Activating: `H3K27ac, H3K4me1, H3K4me3,
H3K36me3, H3K9ac`. Inhibitory: `H3K27me3, H3K9me3`.

### DNA methylation (optional, inhibitory)
Pre-process with `scripts/process_methylation.py`. Accepted raw formats:
- **Bismark coverage** (`bismark_cov`): `chrom start end pct count_M count_U`
- **Bismark cytosine report** (`bismark_report`):
  `chrom pos strand count_M count_U context tri-context`
- **bedGraph** (`bedgraph`): `chrom start end value` (value 0–1 or 0–100)

Output consumed by PACE: a bedGraph `chrom start end beta` with beta ∈ [0, 1].

### Hi-C (optional)
- `.hic` (Juicer; needs `hic-straw`), `.cool`/`.mcool` (needs `cooler`),
  or **BEDPE** `chr1 start1 end1 chr2 start2 end2 score`.
- `HiC_type` ∈ {`hic`, `cool`, `bedpe`, `avg`}; `HiC_resolution` in bp (5000).

### RNA-seq expression (optional)
TSV with a gene-id column and an expression column (default `TPM`):
```
gene_id   gene_name   TPM
ENSSSCG…  GENE1       42.7
```
Accepted id columns: `gene_id, gene_name, GeneID, Gene` (first match used).

### eQTL summary statistics (optional)
TSV with at least chromosome and position; gene and p-value recommended:
```
chr   pos       gene       pvalue   beta
1     10456321  ENSSSCG…   3.1e-08  0.42
```
Column names are configurable (`eqtl_validation.columns` in `config.yaml`;
chromosome/position synonyms auto-detected by `eqtl_enrichment.py`).

### TF ChIP-seq (optional)
`TF_binding`: comma-separated `.bam`/`.bw` files. `TF_names`: matching names.
`TF_modes`: matching `activator`/`repressor` (default activator).

### Gene annotation BED
BED6 **or** extended 8-column BED:
`chrom start end NAME[;ENSID] score strand [gene_id] [gene_type]`.
TSS BED: `chrom start end name score strand` (TSS ± window).

### Validation pairs (for ML / benchmarking)
```
enhancer        gene     validated
chr1:1000-1500  GENE1    1
chr1:5000-5500  GENE2    0
```
`enhancer` may be `chr:start-end` or `chr_start_end`; `validated`/`label` is
1 (positive) or 0 (negative).

---

## Outputs

### `Peaks/candidateRegions.bed`
`chrom start end name [isTSS]` — 500 bp candidate elements.

### `Neighborhoods/EnhancerList.txt` (has header)
`chr start end name <SIGNAL>_signal … activity`
One `<signal>_signal` column per quantified track; `activity` is the
aggregated, inhibition-adjusted enhancer activity.

### `Neighborhoods/GeneList.txt` (has header)
`chr start end name score strand [gene_id gene_type] gene_name TSS [Expression]`

### `Predictions/EnhancerPredictionsAllPutative.tsv.gz` (has header)
All scored E-G pairs within 5 Mb.

| Column | Description |
|--------|-------------|
| `chr, start, end` | enhancer coordinates |
| `name` | enhancer id |
| `TargetGene` | gene symbol |
| `TargetGeneEnsemblID` | Ensembl gene id |
| `TargetGeneTSS` | gene TSS |
| `TargetGeneStrand` | gene strand |
| `distance` | enhancer-midpoint → TSS distance (bp) |
| `class` | `promoter` / `proximal` / `distal` |
| `activity` | enhancer activity |
| `contact` | Hi-C or power-law contact |
| `ABC.Score` | normalized score (per-gene), 0–1 |
| `PACE.Score` | expression-weighted score (standalone calculator) |
| `Expression, ExpressionWeight, isExpressed` | if RNA-seq provided |
| `<signal>` | per-signal component scores (standalone calculator) |
| `eQTL_support, eQTL_pvalue, eQTL_beta` | if eQTL validation enabled |

### `Predictions/EnhancerPredictions.tsv` / `_Full.tsv`
Filtered predictions (score ≥ threshold): slim key columns / all columns.

### `Predictions/EnhancerPredictions_ML.tsv`
Adds `ML.Score` (ML probability), `ML.Decision` (`accept`/`reject`), and
`Combined.Score` — a gated selection (`ML.Score` if the learned weights agree
with the default priors, otherwise the formula-based `ABC.Score`; never an
average).

### `Metrics/QCSummary_*.tsv`, `Metrics/QCPlots_*.pdf`
Per-sample summary statistics and QC plots (score/distance distributions,
score-vs-distance, enhancer class, activity, enhancers per gene).

### Analysis-module outputs
| Tool | Output |
|------|--------|
| `eqtl_enrichment.py` | TSV: contingency counts, OR + 95% CI, Fisher p, score-stratified fold enrichment |
| `sensitivity_analysis.py` | TSV: per-weight factor, Jaccard/overlap of top-k, Spearman, AUROC |
| `benchmark_compare.py` | `benchmark_summary.tsv` (AUROC/AUPRC/sens@spec) + `benchmark_curves.png` |
| `pace_ml.py importance` | gain & permutation importance TSVs, default-vs-learned TSV, plot |
