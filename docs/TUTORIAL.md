# PACE Tutorial: From Raw Data to Validated E-P Predictions

This tutorial walks through a complete PACE analysis, documenting the **inputs
and outputs of every step**. For a formal specification of every file format,
see [IO_FORMATS.md](IO_FORMATS.md).

Contents
1. [Quick start (5 minutes)](#1-quick-start)
2. [Preparing a reference](#2-preparing-a-reference)
3. [Configuring samples](#3-configuring-samples)
4. [The five pipeline steps](#4-the-five-pipeline-steps)
5. [Running with Snakemake](#5-running-with-snakemake)
6. [Adding multi-omics layers](#6-adding-multi-omics-layers)
7. [Validation and benchmarking](#7-validation-and-benchmarking)
8. [Optional ML module](#8-optional-ml-module)
9. [Robustness analysis](#9-robustness-analysis)
10. [Troubleshooting](#10-troubleshooting)

---

## 1. Quick start

Run the bundled synthetic example end-to-end with **no Snakemake or MACS2**
required (only Python + bedtools):

```bash
bash example/run_example_direct.sh
```

This generates synthetic data (5 Mb genome, 50 genes, 240 peaks), runs all
five steps, and prints the top predictions. Outputs land in
`example/results/Example_Sample/`. Expected runtime: ~1–2 minutes.

To confirm your installation, run the test suite (45 tests):

```bash
python test_pace_complete.py
```

---

## 2. Preparing a reference

PACE needs, per species: chromosome sizes, a gene BED, and a TSS BED. Generate
them from an Ensembl GTF + FASTA:

```bash
python scripts/prepare_reference.py \
    --gtf  Sus_scrofa.Sscrofa11.1.111.gtf.gz \
    --fasta Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz \
    --output_dir reference/pig --species_name pig
```

| Output | Description |
|--------|-------------|
| `genes.bed` | Gene coordinates (`chr start end NAME;ENSID score strand gene_id gene_type`) |
| `genes.TSS500bp.bed` | TSS ± 250 bp promoter windows |
| `genome.chrom.sizes` | `chrom <TAB> length` |

> PACE parses both 6-column BED6 and the extended 8-column gene BED above.

---

## 3. Configuring samples

Edit `config/config_biosamples.tsv` (tab-separated). Minimum: a sample name,
an accessibility file, and the accessibility type.

```tsv
biosample    ATAC                           H3K27ac                 RNA_seq            default_accessibility_feature
Pig_Liver    data/pig_liver_ATAC.tagAlign.gz data/pig_liver_H3K27ac.bam data/pig_liver_rna.tsv ATAC
```

Optional columns: `H3K4me1 H3K4me3 H3K36me3 H3K9ac H3K27me3 H3K9me3
methylation HiC_file HiC_type HiC_resolution TF_binding TF_names TF_modes
eQTL_file alt_TSS alt_genes`. See [IO_FORMATS.md](IO_FORMATS.md) for each.

Edit `config/config.yaml` for pipeline parameters (activity method, weights,
thresholds, genome size for MACS2, etc.). Every default is explained in
[PARAMETERS.md](PARAMETERS.md).

---

## 4. The five pipeline steps

You can run the steps individually (useful for debugging and for understanding
the I/O). The `example/run_example_direct.sh` script does exactly this.

### Step 1 — Candidate regions
**In:** MACS2 narrowPeak + chrom.sizes → **Out:** `candidateRegions.bed`

```bash
python workflow/scripts/pace_candidate_regions.py \
    --narrowPeak peaks.narrowPeak --chrom_sizes genome.chrom.sizes \
    --output Peaks/candidateRegions.bed \
    --nStrongestPeaks 150000 --peakExtendFromSummit 250
```
Selects the strongest peaks, extends summits to 500 bp windows, removes
blacklist/TSS overlaps.

### Step 2 — Neighborhoods (activity)
**In:** candidate regions, gene BED, accessibility (+ optional histone/methylation),
optional RNA-seq → **Out:** `EnhancerList.txt`, `GeneList.txt`

```bash
python workflow/scripts/pace_neighborhoods.py \
    --candidate_regions Peaks/candidateRegions.bed \
    --genes reference/pig/genes.bed \
    --chrom_sizes reference/pig/genome.chrom.sizes \
    --output_dir Neighborhoods/ \
    --accessibility_file data/pig_liver_ATAC.tagAlign.gz \
    --H3K27ac data/pig_liver_H3K27ac.bam \
    --H3K27me3 data/pig_liver_H3K27me3.bam \
    --methylation data/pig_liver_meth.bedGraph \
    --activity_method weighted_geometric
```
`EnhancerList.txt` gains one column per quantified signal plus an `activity`
column. Repressive marks (`--H3K27me3`, `--H3K9me3`) and `--methylation` are
applied as **inhibitory** signals that down-weight activity.

### Step 3 — Predict
**In:** `EnhancerList.txt`, `GeneList.txt`, optional Hi-C/expression →
**Out:** `EnhancerPredictionsAllPutative.tsv.gz`

```bash
python workflow/scripts/pace_predict.py \
    --enhancers Neighborhoods/EnhancerList.txt \
    --genes Neighborhoods/GeneList.txt \
    --output Predictions/EnhancerPredictionsAllPutative.tsv.gz \
    --hic_file data/pig_liver.hic --hic_type hic \
    --expression data/pig_liver_rna.tsv --use_expression_weight \
    --threshold 0.02
```
Forms all E-G pairs within 5 Mb, estimates contact (Hi-C if given, else
power-law), computes the normalized PACE score.

### Step 4 — Filter
**In:** all-putative predictions → **Out:** `EnhancerPredictions.tsv`
(+ `_Full.tsv`)

```bash
python workflow/scripts/pace_filter.py \
    --predictions Predictions/EnhancerPredictionsAllPutative.tsv.gz \
    --output Predictions/EnhancerPredictions.tsv --threshold 0.02
```

### Step 5 — QC metrics
**In:** predictions → **Out:** `QCSummary_*.tsv`, `QCPlots_*.pdf`

```bash
python workflow/scripts/pace_metrics.py \
    --predictions Predictions/EnhancerPredictionsAllPutative.tsv.gz \
    --output_dir Metrics/ --sample_name Pig_Liver
```

---

## 5. Running with Snakemake

For real datasets (with MACS2 peak calling and parallelism):

```bash
snakemake --cores 8                       # full pipeline for all samples
snakemake --cores 8 peaks_only            # stop after peak calling
snakemake --cores 8 neighborhoods_only    # stop after activity
```

Outputs follow `results/{biosample}/{Peaks,Neighborhoods,Predictions,Metrics}/`.

---

## 6. Adding multi-omics layers

PACE scales from minimal to comprehensive inputs:

| Configuration | Columns to add | Effect |
|---------------|----------------|--------|
| Minimal | `ATAC` | accessibility-only activity |
| Standard | `+ H3K27ac` | active-enhancer signal |
| Expression-aware | `+ RNA_seq` | expression weighting/filter |
| Repression-aware | `+ methylation`, `H3K27me3`, `H3K9me3` | inhibitory down-weighting |
| 3D-aware | `+ HiC_file/HiC_type/HiC_resolution` | measured contacts |
| Mechanism | `+ TF_binding/TF_names/TF_modes` | activator/repressor TFs |

### Methylation (WGBS/RRBS) preprocessing

Convert per-cytosine calls to a methylation track before use:

```bash
python scripts/process_methylation.py \
    --input sample.bismark.cov.gz --format bismark_cov \
    --context CG --min_coverage 4 \
    --regions Peaks/candidateRegions.bed \
    --output data/sample_region_methylation.bedGraph
```

Then set the `methylation` column to this file.

---

## 7. Validation and benchmarking

### eQTL enrichment (with full contingency tables)

```bash
python scripts/eqtl_enrichment.py \
    --predictions Predictions/EnhancerPredictionsAllPutative.tsv.gz \
    --eqtl pig_liver_eqtl.tsv --chrom_sizes reference/pig/genome.chrom.sizes \
    --score_col ABC.Score --pvalue_col pvalue --label Pig_Liver \
    --output validation/pig_liver_eqtl_enrichment.tsv
```
Reports the 2×2 table counts, odds ratio + 95% CI, Fisher p-value, and
score-stratified fold enrichment.

### Comparing methods on a gold standard

```bash
python scripts/benchmark_compare.py \
    --labels crispri_validated.tsv \
    --method PACE:pace_pred.tsv:PACE.Score \
    --method ABC:abc_pred.tsv:ABC.Score \
    --method EPIPDLF:epipdlf.tsv:score \
    --method GATv2EPI:gatv2.tsv:score \
    --output benchmark/
```
Produces an AUROC/AUPRC table and overlaid ROC/PR curves on the common
evaluable set.

---

## 8. Optional ML module

Use ML **only when species-/context-matched validated E-P pairs are
available** (e.g. tissue-matched eQTL or CRISPRi). PACE does not transfer human
models blindly to livestock.

```bash
python scripts/pace_ml.py train \
    --predictions Predictions/EnhancerPredictionsAllPutative.tsv.gz \
    --validation matched_validated_pairs.tsv \
    --output models/pig_liver.pkl --balance_classes

python scripts/pace_ml.py predict \
    --predictions Predictions/EnhancerPredictionsAllPutative.tsv.gz \
    --model models/pig_liver.pkl --output Predictions/EnhancerPredictions_ML.tsv

python scripts/pace_ml.py importance \
    --model models/pig_liver.pkl --output models/pig_liver_importance \
    --predictions Predictions/EnhancerPredictionsAllPutative.tsv.gz \
    --validation matched_validated_pairs.tsv
```

---

## 9. Robustness analysis

Confirm that predictions are stable to the weight choices (±50%):

```bash
python scripts/sensitivity_analysis.py \
    --predictions predictions_with_signals.tsv \
    --output robustness/sensitivity.tsv --plot robustness/sensitivity.png \
    --validated_col validated
```

---

## 10. Troubleshooting

| Symptom | Cause / fix |
|---------|-------------|
| `bedtools: not found` | Install bedtools (`conda install -c bioconda bedtools`). |
| All-zero activity | Accessibility file does not overlap candidate regions; check chromosome naming (`chr1` vs `1`). |
| `pyBigWig` ImportError on bigWig input | `pip install pyBigWig`, or supply BED/tagAlign/BAM. |
| Empty predictions | Gene BED and accessibility on different chromosome names; verify both use the same assembly. |
| Hi-C ignored | Install `hic-straw` (.hic) or `cooler` (.cool); otherwise PACE falls back to power-law. |
| ML "Too few validated pairs" | Need ≥10 matched validated positives; otherwise use the formula-based score. |
