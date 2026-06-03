# PACE Example Tutorial

This tutorial walks you through running PACE with example data to help you understand the complete workflow.

## Overview

The example dataset includes:
- **Synthetic genome**: 5 Mb chromosome (chr1)
- **50 genes**: With simulated promoters and enhancers
- **ATAC-seq data**: ~60,000 simulated reads
- **H3K27ac data**: ~58,000 simulated reads

Running time: **~2-5 minutes** on a typical laptop

## Quick Run

```bash
# One command to run everything
bash example/run_example.sh
```

## Step-by-Step Guide

### Step 1: Generate Example Data

```bash
# Generate synthetic test data
python scripts/generate_example_data.py
```

This creates:
```
example/
├── data/
│   ├── example_ATAC.tagAlign.gz      # Simulated ATAC-seq reads
│   └── example_H3K27ac.tagAlign.gz   # Simulated H3K27ac reads
└── reference/
    ├── example_genome.fa.gz          # Reference genome (5 Mb)
    ├── example_annotation.gtf.gz     # Gene annotations (50 genes)
    └── example.chrom.sizes           # Chromosome sizes
```

### Step 2: Prepare Reference Files

```bash
# Generate all required reference files from GTF and FASTA
python scripts/prepare_reference.py \
    --gtf example/reference/example_annotation.gtf.gz \
    --fasta example/reference/example_genome.fa.gz \
    --output_dir example/reference \
    --species_name example
```

This creates:
```
example/reference/
├── example.chrom.sizes    # Chromosome sizes
├── genes.bed              # Gene annotation BED file
├── genes.TSS500bp.bed     # TSS regions (±500bp)
└── blacklist.bed          # Regions to exclude
```

### Step 3: Review Configuration

**Example config file** (`example/config.yaml`):

```yaml
# Output directory
results_dir: "example/results/"

# Sample configuration
biosamplesTable: "example/config_biosamples.tsv"

# Reference files
ref:
  chrom_sizes: "example/reference/example.chrom.sizes"
  genes: "example/reference/genes.bed"
  genome_tss: "example/reference/genes.TSS500bp.bed"
  regions_blocklist: "example/reference/blacklist.bed"

# MACS2 parameters (small genome)
params_macs:
  genome_size: 5.0e6  # 5 Mb for example
  
# Candidate regions (reduced for example)
params_candidate:
  nStrongestPeaks: 5000
```

**Sample configuration** (`example/config_biosamples.tsv`):

| biosample | ATAC | H3K27ac | default_accessibility_feature |
|-----------|------|---------|-------------------------------|
| Example_Sample | example/data/example_ATAC.tagAlign.gz | example/data/example_H3K27ac.tagAlign.gz | ATAC |

### Step 4: Run the Pipeline

```bash
# Dry run (check what will be executed)
snakemake --configfile example/config.yaml -n

# Run with 4 cores
snakemake --configfile example/config.yaml --cores 4
```

### Step 5: Explore Results

```bash
# List output files
ls -la example/results/Example_Sample/

# View directory structure
tree example/results/
```

**Output structure:**
```
example/results/Example_Sample/
├── Peaks/
│   ├── macs2_peaks.narrowPeak           # Raw peaks
│   └── *.candidateRegions.bed           # Candidate enhancer regions
├── Neighborhoods/
│   ├── EnhancerList.txt                 # Enhancer list with activity
│   └── GeneList.txt                     # Gene list with expression
├── Predictions/
│   ├── EnhancerPredictionsAllPutative.tsv.gz    # All predictions
│   ├── EnhancerPredictionsFull_*.tsv            # Filtered predictions
│   └── EnhancerPredictions_*.tsv                # Slim predictions
└── Metrics/
    ├── QCSummary_*.tsv                  # QC summary
    └── QCPlots_*.pdf                    # QC plots
```

### Step 6: Analyze Predictions

```bash
# View first few predictions
zcat example/results/Example_Sample/Predictions/EnhancerPredictionsAllPutative.tsv.gz | head -5

# Count predictions
zcat example/results/Example_Sample/Predictions/EnhancerPredictionsAllPutative.tsv.gz | wc -l

# View high-confidence predictions
cat example/results/Example_Sample/Predictions/EnhancerPredictions_*.tsv | head -10

# Get predictions for a specific gene
zcat example/results/Example_Sample/Predictions/EnhancerPredictionsAllPutative.tsv.gz | grep "GENE1"
```

**Key columns in output:**

| Column | Description |
|--------|-------------|
| chr | Enhancer chromosome |
| start | Enhancer start position |
| end | Enhancer end position |
| TargetGene | Predicted target gene |
| ABC.Score | ABC score (0-1, higher = more confident) |
| distance | Distance from enhancer to TSS |
| class | Enhancer class (promoter/genic/intergenic) |

## Understanding the Output

### ABC Score Interpretation

| ABC Score | Interpretation |
|-----------|----------------|
| > 0.1 | High confidence |
| 0.02 - 0.1 | Medium confidence |
| < 0.02 | Low confidence (filtered by default) |

### Enhancer Classes

| Class | Description |
|-------|-------------|
| promoter | Overlaps gene TSS region |
| genic | Within gene body |
| intergenic | Between genes |

## Common Issues

### Issue 1: "ABC core scripts not found"

```bash
# Run setup script to download required scripts
bash setup.sh
```

### Issue 2: Memory errors

```bash
# Reduce parallelization
snakemake --configfile example/config.yaml --cores 1
```

### Issue 3: Missing dependencies

```bash
# Reinstall conda environment
conda env remove -n pace-env
conda env create -f workflow/envs/pace-env.yml
conda activate pace-env
```

## Next Steps

After successfully running the example, you're ready to analyze your own data:

1. **Prepare your data**: BAM/tagAlign files for ATAC-seq and H3K27ac
2. **Download reference files**: See README.md for species-specific links
3. **Generate reference files**: Use `scripts/prepare_reference.py`
4. **Configure your samples**: Edit `config/config_biosamples.tsv`
5. **Run PACE**: `snakemake --cores 8`

See [README.md](../README.md) for detailed documentation.
