#!/bin/bash
# PACE Example Runner
# Run the complete pipeline with example data
#
# Usage: bash example/run_example.sh
#
# Author: Linyong Shen @ Northwest A&F University

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo ""
echo -e "${BLUE}╔══════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║                                                          ║${NC}"
echo -e "${BLUE}║   ██████╗  █████╗  ██████╗███████╗                      ║${NC}"
echo -e "${BLUE}║   ██╔══██╗██╔══██╗██╔════╝██╔════╝                      ║${NC}"
echo -e "${BLUE}║   ██████╔╝███████║██║     █████╗                        ║${NC}"
echo -e "${BLUE}║   ██╔═══╝ ██╔══██║██║     ██╔══╝                        ║${NC}"
echo -e "${BLUE}║   ██║     ██║  ██║╚██████╗███████╗                      ║${NC}"
echo -e "${BLUE}║   ╚═╝     ╚═╝  ╚═╝ ╚═════╝╚══════╝                      ║${NC}"
echo -e "${BLUE}║                                                          ║${NC}"
echo -e "${BLUE}║   Example Pipeline Runner                                ║${NC}"
echo -e "${BLUE}║                                                          ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════╝${NC}"
echo ""

# Check if we're in the PACE directory
if [ ! -f "workflow/Snakefile" ]; then
    echo -e "${RED}Error: Please run this script from the PACE root directory${NC}"
    echo "Usage: bash example/run_example.sh"
    exit 1
fi

# Check if example data exists
if [ ! -f "example/data/example_ATAC.tagAlign.gz" ]; then
    echo -e "${YELLOW}[Step 1/4] Generating example data...${NC}"
    python scripts/generate_example_data.py
else
    echo -e "${GREEN}[Step 1/4] Example data already exists, skipping generation${NC}"
fi

# Check if reference files exist
if [ ! -f "example/reference/genes.bed" ]; then
    echo -e "${YELLOW}[Step 2/4] Preparing reference files...${NC}"
    python scripts/prepare_reference.py \
        --gtf example/reference/example_annotation.gtf.gz \
        --fasta example/reference/example_genome.fa.gz \
        --output_dir example/reference \
        --species_name example
else
    echo -e "${GREEN}[Step 2/4] Reference files already exist, skipping preparation${NC}"
fi

# Check if ABC scripts exist
echo ""
echo -e "${YELLOW}[Step 3/4] Checking ABC core scripts...${NC}"
if [ ! -f "workflow/scripts/predict.py" ]; then
    echo -e "${RED}Error: ABC core scripts not found!${NC}"
    echo ""
    echo "Please run setup.sh first to download the required scripts:"
    echo "  bash setup.sh"
    echo ""
    echo "Or manually download from:"
    echo "  https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction"
    echo ""
    exit 1
else
    echo -e "${GREEN}  ABC core scripts found${NC}"
fi

# Run the pipeline
echo ""
echo -e "${YELLOW}[Step 4/4] Running PACE pipeline...${NC}"
echo ""
echo "Configuration:"
echo "  Config file: example/config.yaml"
echo "  Sample: Example_Sample"
echo "  ATAC-seq: example/data/example_ATAC.tagAlign.gz"
echo "  H3K27ac: example/data/example_H3K27ac.tagAlign.gz"
echo ""

# Determine number of cores
CORES=${1:-4}
echo "Using ${CORES} cores (change with: bash example/run_example.sh <num_cores>)"
echo ""

# Run snakemake
snakemake \
    --configfile example/config.yaml \
    --cores ${CORES} \
    --printshellcmds \
    --reason \
    --keep-going

# Check results
echo ""
echo -e "${GREEN}╔══════════════════════════════════════════════════════════╗${NC}"
echo -e "${GREEN}║                   Pipeline Complete!                     ║${NC}"
echo -e "${GREEN}╚══════════════════════════════════════════════════════════╝${NC}"
echo ""
echo "Results saved to: example/results/Example_Sample/"
echo ""
echo "Key output files:"
echo ""

if [ -f "example/results/Example_Sample/Predictions/EnhancerPredictionsAllPutative.tsv.gz" ]; then
    echo -e "${GREEN}✓${NC} EnhancerPredictionsAllPutative.tsv.gz"
    NUM_PREDICTIONS=$(zcat example/results/Example_Sample/Predictions/EnhancerPredictionsAllPutative.tsv.gz | wc -l)
    echo "  Total predictions: $((NUM_PREDICTIONS - 1))"
else
    echo -e "${RED}✗${NC} EnhancerPredictionsAllPutative.tsv.gz (not found)"
fi

# Show filtered predictions if exists
for f in example/results/Example_Sample/Predictions/EnhancerPredictions_*.tsv; do
    if [ -f "$f" ]; then
        echo -e "${GREEN}✓${NC} $(basename $f)"
        NUM_FILTERED=$(wc -l < "$f")
        echo "  Filtered predictions: $((NUM_FILTERED - 1))"
    fi
done

# Show QC files if exist
for f in example/results/Example_Sample/Metrics/QC*.tsv; do
    if [ -f "$f" ]; then
        echo -e "${GREEN}✓${NC} $(basename $f)"
    fi
done

echo ""
echo "View predictions:"
echo "  zcat example/results/Example_Sample/Predictions/EnhancerPredictionsAllPutative.tsv.gz | head -20"
echo ""
