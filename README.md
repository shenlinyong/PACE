# PACE: Predicting Activity-Contact for Enhancer-promoter

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![GitHub Stars](https://img.shields.io/github/stars/shenlinyong/PACE?style=social)](https://github.com/shenlinyong/PACE)

```
╔══════════════════════════════════════════════════════════╗
║                                                          ║
║   ██████╗  █████╗  ██████╗███████╗                       ║
║   ██╔══██╗██╔══██╗██╔════╝██╔════╝                       ║
║   ██████╔╝███████║██║     █████╗                         ║
║   ██╔═══╝ ██╔══██║██║     ██╔══╝                         ║
║   ██║     ██║  ██║╚██████╗███████╗                       ║
║   ╚═╝     ╚═╝  ╚═╝ ╚═════╝╚══════╝                       ║
║                                                          ║
║   Predicting Activity-Contact for Enhancer-promoter      ║
║                                                          ║
╚══════════════════════════════════════════════════════════╝
```

**Author**: Linyong Shen @ Northwest A&F University (西北农林科技大学)

## Overview

**PACE** is a flexible multi-omics computational framework for predicting enhancer-promoter regulatory interactions in livestock species, based on the [Activity-by-Contact (ABC) Model](https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction).

### Core Innovation

PACE extends the ABC model with flexible multi-omics integration:

```
                    A(E) × C(E,G) × W_expr(G)
PACE Score(E,G) = ─────────────────────────────
                   Σ[A(e) × C(e,G) × W_expr(G)]
```

Where Activity can integrate multiple data types:
```
Activity = f(Accessibility, Histones, TFs, Methylation, ...)
```

### Key Features

- 🐄 **Universal for Livestock**: Pig, cattle, sheep, chicken, horse, and other domestic animals
- 📊 **Flexible Multi-omics Integration**: From minimal (ATAC-only) to comprehensive (all data types)
- 🧬 **RNA-seq Integration**: Filter predictions to expressed genes
- 🔬 **Methylation Support**: DNA methylation as inhibitory signal
- 🎯 **TF Binding**: Incorporate transcription factor ChIP-seq
- ✅ **eQTL Validation**: Validate predictions with genetic evidence
- 🚀 **Automated Pipeline**: Snakemake workflow for one-command execution

## Documentation

| Document | Description |
|----------|-------------|
| 📖 [Quick Start Guide](docs/QUICKSTART.md) | 5-minute guide to get started |
| 📊 [Methods & Algorithm](docs/METHODS.md) | Detailed methods and formulas |
| 🤖 [ML Integration](docs/ML_INTEGRATION.md) | Machine learning module guide |

## Quick Start

### Installation

```bash
# Clone the repository
git clone https://github.com/shenlinyong/PACE.git
cd PACE

# Run setup script
bash setup.sh

# Activate environment
conda activate pace-env
```

### System Requirements

| Component | Version | Required |
|-----------|---------|----------|
| Python | 3.8 - 3.12 | ✓ |
| Conda/Mamba | Latest | ✓ |
| Memory | 8GB+ RAM | ✓ |
| Storage | 5GB+ | ✓ |

### Dependencies

#### Core (Required)
| Package | Version | Purpose |
|---------|---------|---------|
| numpy | ≥1.20 | Numerical computing |
| pandas | ≥1.3 | Data manipulation |
| scipy | ≥1.7 | Scientific computing |
| matplotlib | ≥3.4 | Visualization |
| seaborn | ≥0.11 | Statistical plots |

#### Bioinformatics (Required)
| Package | Version | Purpose |
|---------|---------|---------|
| bedtools | ≥2.30 | BED file operations |
| samtools | ≥1.15 | BAM processing |
| macs2 | ≥2.2 | Peak calling |
| pybigwig | ≥0.3 | BigWig I/O |
| pybedtools | ≥0.9 | BED operations |

#### Optional
| Package | Version | Purpose |
|---------|---------|---------|
| scikit-learn | ≥1.0 | Machine learning |
| hic-straw | ≥1.3 | Hi-C data |
| snakemake | 7.x | Workflow |

#### Alternative Installation (pip)
```bash
# Install system tools first (via conda or apt)
conda install -c bioconda bedtools samtools macs2

# Then install Python packages
pip install -r requirements.txt
```

### Verify Test Environment Setup

To confirm that the testing environment is correctly configured, run the full test suite:

```bash
# Run comprehensive test suite (37 tests)
python test_pace_complete.py
```

This command evaluates the following components:
- Module imports (7 tests)
- Core functions (9 tests)
- CLI scripts (6 tests)
- Complete pipeline (5 tests)
- Machine learning module (10 tests)

Successful completion of all tests indicates that the environment has been properly set up and that all major components of the PACE framework are functioning as expected.

### Try with Example Data (Recommended for First-Time Users)

We provide a complete example dataset to help you get started quickly:

```bash
# Run the example pipeline (takes ~2-5 minutes)
bash example/run_example.sh
```

This will:
1. Generate synthetic test data (5 Mb genome, 50 genes)
2. Prepare reference files
3. Run the complete PACE pipeline
4. Show prediction results

See [example/README.md](example/README.md) for a detailed step-by-step tutorial.

### Basic Usage

```bash
# 1. Prepare reference files
python scripts/prepare_reference.py \
    --gtf annotation.gtf.gz \
    --fasta genome.fa.gz \
    --output_dir reference/my_species

# 2. Configure your analysis:
#    - Edit config/config.yaml for pipeline parameters
#    - Edit config/config_biosamples.tsv for sample information

# 3. Run the pipeline
snakemake --cores 8
```

## Supported Data Types

### Required (at least one)

| Data Type | Format | Description |
|-----------|--------|-------------|
| **ATAC-seq** | BAM/tagAlign | Chromatin accessibility |
| **DNase-seq** | BAM | Chromatin accessibility |

### Recommended

| Data Type | Format | Description | Effect |
|-----------|--------|-------------|--------|
| **H3K27ac** | BAM/bigWig | Active enhancers | Activating |
| **RNA-seq** | TSV (TPM) | Gene expression | Filter/Weight |

### Optional (Enhanced Predictions)

| Data Type | Format | Description | Effect |
|-----------|--------|-------------|--------|
| **H3K4me1** | BAM/bigWig | Enhancer mark | Activating |
| **H3K4me3** | BAM/bigWig | Promoter mark | Activating |
| **H3K36me3** | BAM/bigWig | Transcribed regions | Activating |
| **H3K27me3** | BAM/bigWig | Polycomb repression | **Inhibitory** |
| **H3K9me3** | BAM/bigWig | Heterochromatin | **Inhibitory** |
| **DNA Methylation** | bedGraph/bigWig | CpG methylation | **Inhibitory** |
| **TF ChIP-seq** | BAM/bigWig | TF binding | Activating |
| **Hi-C** | .hic/bedpe | 3D chromatin contact | Contact |
| **eQTL** | TSV | Genetic evidence | Validation |

## Configuration

### Activity Calculation Methods

PACE supports multiple methods for combining signals:

| Method | Formula | Use Case |
|--------|---------|----------|
| `geometric_mean` | ∏(Sᵢ)^(1/n) | Simple, like original ABC |
| `weighted_geometric` | ∏(Sᵢ^wᵢ) | Prioritize certain signals |
| `weighted_sum` | Σ(wᵢ × Sᵢ) | Linear combination |

### Signal Weights (Defaults)

| Signal | Type | Default Weight |
|--------|------|----------------|
| ATAC/DNase | Accessibility | 1.5 |
| H3K27ac | Active enhancer | 1.0 |
| H3K4me1 | Enhancer mark | 0.8 |
| H3K4me3 | Promoter mark | 0.5 |
| H3K27me3 | Repressive | 0.5 (inhibitory) |
| Methylation | CpG | 0.5 (inhibitory) |
| TF binding | Transcription factor | 0.3 |

### Example Configurations

The following examples show different configuration levels for `config/config.yaml`:

#### Minimal (ATAC only)

Edit `config/config.yaml`:
```yaml
activity_method: "geometric_mean"
histone_marks:
  H3K27ac:
    enabled: false
```

#### Standard (Recommended)

Edit `config/config.yaml`:
```yaml
activity_method: "geometric_mean"
histone_marks:
  H3K27ac:
    enabled: true
    weight: 1.0
```

#### Multi-omics Enhanced

Edit `config/config.yaml`:
```yaml
activity_method: "weighted_geometric"
accessibility:
  weight: 1.5
histone_marks:
  H3K27ac: {enabled: true, weight: 1.0}
  H3K4me1: {enabled: true, weight: 0.8}
  H3K4me3: {enabled: true, weight: 0.5}
methylation:
  enabled: true
  weight: 0.5
expression:
  enabled: true
  min_expression: 1.0
  weight_method: "log"
```

## Sample Configuration

Edit `config/config_biosamples.tsv` to add your samples. The file uses tab-separated format with the following columns:

### Required Columns
| Column | Description |
|--------|-------------|
| `biosample` | Sample name |
| `default_accessibility_feature` | `ATAC` or `DHS` |

### Accessibility (at least one required)
| Column | Description |
|--------|-------------|
| `ATAC` | ATAC-seq tagAlign/BAM file |
| `DHS` | DNase-seq BAM file |

### Optional Columns
| Column | Description |
|--------|-------------|
| `H3K27ac`, `H3K4me1`, `H3K4me3`, etc. | Histone ChIP-seq files |
| `RNA_seq` | Gene expression file (TPM) |
| `methylation` | DNA methylation file |
| `HiC_file`, `HiC_type`, `HiC_resolution` | Hi-C data |
| `TF_binding`, `TF_names` | TF ChIP-seq files |

### Example Entry
```tsv
biosample	ATAC	H3K27ac	RNA_seq	default_accessibility_feature
Pig_Liver	data/atac.tagAlign.gz	data/h3k27ac.bam	data/rnaseq.tsv	ATAC
```

See `config/config_biosamples.tsv` for complete column list and more examples.

## Output Files

```
results/{biosample}/
├── Peaks/
│   ├── macs2_peaks.narrowPeak              # MACS2 peak calls
│   └── macs2_peaks.narrowPeak.sorted.candidateRegions.bed
├── Neighborhoods/
│   ├── EnhancerList.txt                    # Enhancers with activity scores
│   └── GeneList.txt                        # Genes with expression
├── Predictions/
│   ├── EnhancerPredictionsAllPutative.tsv.gz    # All E-G predictions
│   ├── EnhancerPredictionsFull_*.tsv            # Filtered (all columns)
│   └── EnhancerPredictions_*.tsv                # Filtered (key columns)
├── Metrics/
│   ├── QCSummary_*.tsv                     # QC metrics summary
│   └── QCPlots_*.pdf                       # QC visualization
└── logs/
    └── *.log                               # Pipeline logs
```

### Key Output Columns

| Column | Description |
|--------|-------------|
| chr, start, end | Enhancer coordinates |
| TargetGene | Target gene symbol |
| ABC.Score | Prediction score (0-1) |
| distance | Enhancer-TSS distance |
| class | Enhancer class (promoter/proximal/distal) |

## Species-Specific References

### Pig (Sus scrofa)
```bash
wget https://ftp.ensembl.org/pub/release-111/fasta/sus_scrofa/dna/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/release-111/gtf/sus_scrofa/Sus_scrofa.Sscrofa11.1.111.gtf.gz
```

### Cattle (Bos taurus)
```bash
wget https://ftp.ensembl.org/pub/release-111/fasta/bos_taurus/dna/Bos_taurus.ARS-UCD1.3.dna.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/release-111/gtf/bos_taurus/Bos_taurus.ARS-UCD1.3.111.gtf.gz
```

### Sheep (Ovis aries)
```bash
wget https://ftp.ensembl.org/pub/release-111/fasta/ovis_aries_rambouillet/dna/Ovis_aries_rambouillet.Oar_rambouillet_v1.0.dna.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/release-111/gtf/ovis_aries_rambouillet/Ovis_aries_rambouillet.Oar_rambouillet_v1.0.111.gtf.gz
```

### Chicken (Gallus gallus)
```bash
wget https://ftp.ensembl.org/pub/release-111/fasta/gallus_gallus/dna/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/release-111/gtf/gallus_gallus/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.111.gtf.gz
```

### Horse (Equus caballus)
```bash
wget https://ftp.ensembl.org/pub/release-111/fasta/equus_caballus/dna/Equus_caballus.EquCab3.0.dna.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/release-111/gtf/equus_caballus/Equus_caballus.EquCab3.0.111.gtf.gz
```

### Genome Sizes for MACS2

| Species | Genome Size |
|---------|-------------|
| Pig | 2.5e9 |
| Cattle | 2.7e9 |
| Sheep | 2.6e9 |
| Chicken | 1.0e9 |
| Horse | 2.5e9 |

## FAQ

**Q: Can I run PACE with only ATAC-seq data?**  
A: Yes, PACE works with minimal data (ATAC-seq only) and scales up with more data types.

**Q: Can I run without Hi-C data?**  
A: Yes, PACE uses a power-law function to estimate 3D contact when Hi-C is unavailable.

**Q: How do I choose the score threshold?**  
A: Default 0.02 works well. Increase to 0.03-0.05 for higher precision, decrease to 0.015 for higher recall.

**Q: How does methylation affect predictions?**  
A: DNA methylation is treated as an inhibitory signal - high methylation reduces enhancer activity scores.

## Machine Learning Integration (Experimental)

PACE includes an optional ML module that can learn optimal feature combinations from validated E-G pairs (e.g., eQTL data).

### Model Selection

**Default**: Gradient Boosting Classifier
- Provides feature importance for interpretability
- Handles non-linear feature interactions
- Robust to class imbalance (few validated pairs vs many predictions)
- Widely used in genomics applications

**Alternative**: Random Forest (`--model_type random_forest`)

See [docs/ML_INTEGRATION.md](docs/ML_INTEGRATION.md) for detailed model comparison.

### When to Use ML

- You have validated E-G pairs (eQTL, CRISPRi, etc.)
- You want to optimize predictions for your specific tissue/species
- You want to understand which features are most predictive

### Usage

#### Step 1: Prepare Validation Data

Create a TSV file with validated E-G pairs:
```tsv
enhancer	gene	validated
chr1:1000-1500	Gene1	1
chr1:5000-5500	Gene2	1
chr2:3000-3500	Gene3	0
```

#### Step 2: Train Model

```bash
python scripts/pace_ml.py train \
    --predictions results/Sample/Predictions/EnhancerPredictionsAllPutative.tsv.gz \
    --validation data/eqtl_validated_pairs.tsv \
    --output models/my_model.pkl \
    --balance_classes
```

#### Step 3: Apply Model

```bash
python scripts/pace_ml.py predict \
    --predictions results/NewSample/Predictions/EnhancerPredictionsAllPutative.tsv.gz \
    --model models/my_model.pkl \
    --output results/NewSample/Predictions/EnhancerPredictions_ML.tsv \
    --abc_weight 0.5
```

### ML Output

The ML module adds:
- `ML.Score`: ML-based prediction probability
- `Combined.Score`: Weighted combination of ABC.Score and ML.Score

### Configuration

In `config/config.yaml`:
```yaml
ml_integration:
  enabled: false          # Set to true to enable
  model_type: "gradient_boosting"  # or "random_forest"
  use_pretrained: false   # Use pre-trained model
```

**Note**: ML is experimental and requires scikit-learn (`pip install scikit-learn`).

## Citation

If you use PACE in your research, please cite:

```bibtex
@software{PACE,
  author = {Shen, Linyong},
  title = {PACE: Predicting Activity-Contact for Enhancer-promoter},
  year = {2026},
  publisher = {GitHub},
  url = {https://github.com/shenlinyong/PACE}
}
```

And the original ABC model papers:

```bibtex
@article{fulco2019activity,
  title={Activity-by-contact model of enhancer--promoter regulation from thousands of CRISPR perturbations},
  author={Fulco, Charles P and Nasser, Joseph and others},
  journal={Nature Genetics},
  volume={51},
  pages={1664--1669},
  year={2019}
}

@article{nasser2021genome,
  title={Genome-wide enhancer maps link risk variants to disease genes},
  author={Nasser, Joseph and Bergman, Drew T and others},
  journal={Nature},
  volume={593},
  pages={238--243},
  year={2021}
}
```

## Contact

- **Author**: Linyong Shen (申林用)
- **Institution**: Northwest A&F University (西北农林科技大学)
- **Email**: [shenlinyong@nwafu.edu.cn]
- **Issues**: [GitHub Issues](https://github.com/shenlinyong/PACE/issues)

## License

MIT License - see [LICENSE](LICENSE) file for details.

## Acknowledgments

PACE is built upon the [ABC-Enhancer-Gene-Prediction](https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction) framework developed by the Broad Institute and Engreitz Lab.
