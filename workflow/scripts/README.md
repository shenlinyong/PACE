# PACE Workflow Scripts

This directory contains Python scripts for the PACE pipeline.

**Note**: These scripts are developed specifically for PACE, inspired by but not copied from the [ABC-Enhancer-Gene-Prediction](https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction) project.

## Core Modules

| Module | Description | Lines |
|--------|-------------|-------|
| `tools.py` | Utility functions (I/O, math, helpers) | ~470 |
| `peaks.py` | Peak processing and candidate region generation | ~400 |
| `neighborhoods.py` | Multi-omics neighborhood analysis | ~490 |
| `hic.py` | Hi-C data processing and contact estimation | ~420 |
| `predictor.py` | PACE score prediction engine | ~430 |
| `metrics.py` | QC metrics and visualization | ~420 |

## Command-Line Scripts

| Script | Description | Usage |
|--------|-------------|-------|
| `pace_candidate_regions.py` | Generate candidate enhancer regions | [See below](#candidate-regions) |
| `pace_neighborhoods.py` | Run neighborhood analysis | [See below](#neighborhoods) |
| `pace_predict.py` | Predict enhancer-gene interactions | [See below](#predictions) |
| `pace_filter.py` | Filter predictions by threshold | [See below](#filtering) |
| `pace_metrics.py` | Generate QC metrics and plots | [See below](#metrics) |

## Usage Examples

### Candidate Regions

```bash
python pace_candidate_regions.py \
    --narrowPeak peaks.narrowPeak \
    --chrom_sizes genome.chrom.sizes \
    --output candidates.bed \
    --blacklist blacklist.bed \
    --tss_regions tss.bed \
    --nStrongestPeaks 150000 \
    --peakExtendFromSummit 250
```

### Neighborhoods

```bash
python pace_neighborhoods.py \
    --candidate_regions candidates.bed \
    --genes genes.bed \
    --chrom_sizes genome.chrom.sizes \
    --output_dir output/Neighborhoods \
    --accessibility_file atac.tagAlign.gz \
    --accessibility_type ATAC \
    --H3K27ac h3k27ac.bam \
    --H3K4me1 h3k4me1.bam \
    --methylation methylation.bw \
    --expression expression.tsv \
    --activity_method geometric_mean
```

### Predictions

```bash
python pace_predict.py \
    --enhancers EnhancerList.txt \
    --genes GeneList.txt \
    --output predictions.tsv.gz \
    --hic_file sample.hic \
    --hic_type hic \
    --expression expression.tsv \
    --use_expression_weight \
    --threshold 0.02
```

### Filtering

```bash
python pace_filter.py \
    --predictions predictions.tsv.gz \
    --output filtered_predictions.tsv \
    --threshold 0.02 \
    --only_expressed
```

### Metrics

```bash
python pace_metrics.py \
    --predictions predictions.tsv.gz \
    --output_dir output/Metrics \
    --sample_name sample1
```

## Key Differences from ABC

PACE extends the ABC model with several enhancements:

1. **Multi-omics Activity Calculation**
   - ABC: `Activity = ATAC × H3K27ac`
   - PACE: `Activity = f(ATAC, H3K27ac, H3K4me1, H3K4me3, ...) × (1 - Inhibition)`

2. **Flexible Aggregation Methods**
   - `geometric_mean`: Simple geometric mean of signals
   - `weighted_geometric`: Weighted geometric mean with custom weights
   - `weighted_sum`: Weighted sum of normalized signals
   - `arithmetic_mean`: Simple arithmetic mean

3. **Inhibitory Signal Support**
   - DNA methylation
   - H3K27me3 (Polycomb repression)
   - H3K9me3 (Heterochromatin)

4. **Expression Integration**
   - Filter predictions to expressed genes
   - Weight predictions by expression level
   - Multiple weight methods: binary, linear, log

5. **Livestock-Focused Design**
   - Pre-configured for pig, cattle, sheep, chicken, horse
   - Optimized defaults for livestock genomes

## Dependencies

These scripts require:
- Python ≥ 3.8
- numpy, pandas, scipy
- pybigwig
- pybedtools
- matplotlib, seaborn (for metrics)

Install via conda:
```bash
conda activate pace-env
```

## Author

Linyong Shen (沈林泳)  
Northwest A&F University (西北农林科技大学)

## Citation

If you use PACE, please cite both:

1. PACE (this work)
2. Original ABC paper: Fulco et al. Nature Genetics 2019
