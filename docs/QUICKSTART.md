# PACE Quick Start Guide

## 1-Minute Quick Start (Example Data)

```bash
# Clone and setup
git clone https://github.com/shenlinyong/PACE.git
cd PACE
bash setup.sh
conda activate pace-env

# Run example pipeline (~2-5 minutes)
bash example/run_example.sh
```

Results will be in `example/results/Example_Sample/`

---

## 5-Minute Quick Start (Your Own Data)

### Step 1: Installation

```bash
git clone https://github.com/shenlinyong/PACE.git
cd PACE
bash setup.sh
conda activate pace-env
```

### Step 2: Prepare Reference Files

```bash
# Download reference (example: pig)
wget https://ftp.ensembl.org/pub/release-111/fasta/sus_scrofa/dna/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/release-111/gtf/sus_scrofa/Sus_scrofa.Sscrofa11.1.111.gtf.gz

# Generate reference files
python scripts/prepare_reference.py \
    --gtf Sus_scrofa.Sscrofa11.1.111.gtf.gz \
    --fasta Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz \
    --output_dir reference/pig \
    --species_name susScr11
```

### Step 3: Configure

1. Update `config/config.yaml`:
```yaml
ref:
  chrom_sizes: "reference/pig/susScr11.chrom.sizes"
  genes: "reference/pig/genes.bed"
  genome_tss: "reference/pig/genes.TSS500bp.bed"
  regions_blocklist: "reference/pig/blacklist.bed"

params_macs:
  genome_size: 2.5e9  # Pig genome size
```

2. Edit `config/config_biosamples.tsv`:
```tsv
biosample	ATAC	H3K27ac	default_accessibility_feature
Pig_Liver	data/pig_liver_ATAC.tagAlign.gz	data/pig_liver_H3K27ac.bam	ATAC
```

### Step 4: Run

```bash
# Dry run (check configuration)
snakemake -n

# Run with 8 cores
snakemake --cores 8
```

### Step 5: View Results

```bash
# View predictions
zcat results/Pig_Liver/Predictions/EnhancerPredictionsAllPutative.tsv.gz | head
```

---

## Configuration Levels

| Level | Data Required | Configuration |
|-------|---------------|---------------|
| **Minimal** | ATAC-seq only | Default config |
| **Standard** | + H3K27ac | Enable H3K27ac |
| **Enhanced** | + More histones | Enable H3K4me1, H3K4me3 |
| **Expression-aware** | + RNA-seq | Enable expression filter |
| **Comprehensive** | + Methylation, TFs | Enable all options |

---

## Output Files

| File | Description |
|------|-------------|
| `EnhancerPredictionsAllPutative.tsv.gz` | All E-G predictions |
| `EnhancerPredictions_*.tsv` | Filtered high-confidence |
| `QCPlots_*.pdf` | Quality control plots |

---

## Common Issues

**"ABC core scripts not found"**
```bash
bash setup.sh  # Downloads required scripts
```

**Memory errors**
```bash
snakemake --cores 2  # Reduce parallelization
```

**Missing dependencies**
```bash
conda env remove -n pace-env
conda env create -f workflow/envs/pace-env.yml
```

---

See [README.md](../README.md) for full documentation.
