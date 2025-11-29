#!/bin/bash
# PACE Setup Script
# Predicting Activity-Contact for Enhancer-promoter
#
# Author: Linyong Shen @ Northwest A&F University

set -e

echo "========================================"
echo "PACE Setup Script"
echo "========================================"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Get script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Check for conda
check_conda() {
    if command -v conda &> /dev/null; then
        echo -e "${GREEN}✓ Conda found${NC}"
        return 0
    else
        echo -e "${RED}✗ Conda not found${NC}"
        echo "Please install Miniconda or Anaconda first:"
        echo "https://docs.conda.io/en/latest/miniconda.html"
        return 1
    fi
}

# Create conda environment
create_env() {
    ENV_FILE="${SCRIPT_DIR}/workflow/envs/pace-env.yml"
    
    if [ -f "$ENV_FILE" ]; then
        echo "Creating conda environment from ${ENV_FILE}..."
        
        # Check if environment already exists
        if conda env list | grep -q "^pace-env "; then
            echo -e "${YELLOW}Environment 'pace-env' already exists.${NC}"
            read -p "Do you want to update it? (y/n) " -n 1 -r
            echo
            if [[ $REPLY =~ ^[Yy]$ ]]; then
                conda env update -f "$ENV_FILE" -n pace-env
                echo -e "${GREEN}✓ Environment updated${NC}"
            fi
        else
            conda env create -f "$ENV_FILE"
            echo -e "${GREEN}✓ Environment created${NC}"
        fi
    else
        echo -e "${RED}Environment file not found: ${ENV_FILE}${NC}"
        echo "Creating environment manually..."
        
        conda create -n pace-env -y python=3.9
        conda activate pace-env
        conda install -y -c bioconda -c conda-forge \
            snakemake>=7.0 \
            pandas>=1.3 \
            numpy>=1.20 \
            scipy>=1.7 \
            pybigwig>=0.3 \
            pybedtools>=0.9 \
            bedtools>=2.30 \
            samtools>=1.15 \
            macs2>=2.2 \
            matplotlib>=3.5 \
            seaborn>=0.11
        
        echo -e "${GREEN}✓ Environment created manually${NC}"
    fi
}

# Verify PACE scripts
verify_scripts() {
    echo "Verifying PACE scripts..."
    
    SCRIPTS_DIR="${SCRIPT_DIR}/workflow/scripts"
    
    required_scripts=(
        "tools.py"
        "peaks.py"
        "neighborhoods.py"
        "hic.py"
        "predictor.py"
        "metrics.py"
        "pace_candidate_regions.py"
        "pace_neighborhoods.py"
        "pace_predict.py"
        "pace_filter.py"
        "pace_metrics.py"
    )
    
    all_found=true
    for script in "${required_scripts[@]}"; do
        if [ -f "${SCRIPTS_DIR}/${script}" ]; then
            echo -e "  ${GREEN}✓${NC} ${script}"
        else
            echo -e "  ${RED}✗${NC} ${script} - NOT FOUND"
            all_found=false
        fi
    done
    
    if [ "$all_found" = true ]; then
        echo -e "${GREEN}✓ All PACE scripts found${NC}"
        return 0
    else
        echo -e "${RED}✗ Some scripts are missing${NC}"
        return 1
    fi
}

# Make scripts executable
make_executable() {
    echo "Making scripts executable..."
    chmod +x "${SCRIPT_DIR}"/workflow/scripts/*.py 2>/dev/null || true
    chmod +x "${SCRIPT_DIR}"/scripts/*.py 2>/dev/null || true
    echo -e "${GREEN}✓ Scripts made executable${NC}"
}

# Create necessary directories
create_dirs() {
    echo "Creating directories..."
    mkdir -p "${SCRIPT_DIR}/results"
    mkdir -p "${SCRIPT_DIR}/reference"
    mkdir -p "${SCRIPT_DIR}/data"
    echo -e "${GREEN}✓ Directories created${NC}"
}

# Verify example data
verify_example() {
    echo "Verifying example data..."
    
    EXAMPLE_DIR="${SCRIPT_DIR}/example"
    
    if [ -f "${EXAMPLE_DIR}/data/example_ATAC.tagAlign.gz" ] && \
       [ -f "${EXAMPLE_DIR}/reference/genes.bed" ]; then
        echo -e "${GREEN}✓ Example data found${NC}"
        return 0
    else
        echo -e "${YELLOW}! Example data not found or incomplete${NC}"
        echo "You can generate example data by running:"
        echo "  python scripts/generate_example_data.py"
        return 0  # Don't fail on missing example data
    fi
}

# Print usage instructions
print_usage() {
    echo ""
    echo "========================================"
    echo "Setup Complete!"
    echo "========================================"
    echo ""
    echo "To activate the environment:"
    echo -e "  ${GREEN}conda activate pace-env${NC}"
    echo ""
    echo "To run the example pipeline:"
    echo -e "  ${GREEN}bash example/run_example.sh${NC}"
    echo ""
    echo "To run with your own data:"
    echo "  1. Prepare reference files for your species"
    echo "  2. Edit config/config.yaml"
    echo "  3. Edit config/config_biosamples.tsv"
    echo -e "  4. Run: ${GREEN}snakemake --cores 8${NC}"
    echo ""
    echo "For more information, see README.md"
    echo "========================================"
}

# Main function
main() {
    echo ""
    
    # Step 1: Check conda
    check_conda || exit 1
    
    # Step 2: Create environment
    create_env
    
    # Step 3: Verify scripts
    verify_scripts || exit 1
    
    # Step 4: Make scripts executable
    make_executable
    
    # Step 5: Create directories
    create_dirs
    
    # Step 6: Verify example data
    verify_example
    
    # Print usage
    print_usage
}

# Run main
main
