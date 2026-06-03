#!/bin/bash
# PACE Setup Script
# Predicting Activity-Contact for Enhancer-promoter
#
# Author: Linyong Shen @ Northwest A&F University
#
# This script prepares a working PACE environment. It supports two paths:
#   1. Conda  (recommended for the full Snakemake/MACS2 workflow)
#   2. pip + system tools  (dependency-light path: Python + bedtools only,
#      matching example/run_example_direct.sh)
# When conda is not available it automatically falls back to path 2.

set -e

echo "========================================"
echo "PACE Setup Script"
echo "========================================"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

ok()   { echo -e "  ${GREEN}✓${NC} $1"; }
warn() { echo -e "  ${YELLOW}!${NC} $1"; }
err()  { echo -e "  ${RED}✗${NC} $1"; }

# Get script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Python interpreter to use throughout
PYTHON="$(command -v python3 || command -v python || true)"

# Set by check_conda(); selects the environment-setup path.
USE_CONDA=false

# Check for conda (no longer fatal: we fall back to pip if absent)
check_conda() {
    if command -v conda &> /dev/null; then
        echo -e "${GREEN}✓ Conda found${NC} - using conda environment path"
        USE_CONDA=true
    else
        echo -e "${YELLOW}! Conda not found${NC} - using pip + system tools (dependency-light path)"
        echo "  (To use the full Snakemake/MACS2 workflow, install Miniconda:"
        echo "   https://docs.conda.io/en/latest/miniconda.html )"
        USE_CONDA=false
    fi
}

# Create conda environment (conda path)
create_env_conda() {
    ENV_FILE="${SCRIPT_DIR}/workflow/envs/pace-env.yml"

    if [ -f "$ENV_FILE" ]; then
        echo "Creating conda environment from ${ENV_FILE}..."

        # Check if environment already exists
        if conda env list | grep -q "^pace-env "; then
            echo -e "${YELLOW}Environment 'pace-env' already exists.${NC}"
            # Default to "update" when running non-interactively (e.g. CI).
            if [ -t 0 ]; then
                read -p "Do you want to update it? (y/n) " -n 1 -r
                echo
            else
                REPLY="y"
            fi
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
        conda install -n pace-env -y -c bioconda -c conda-forge \
            "snakemake>=7.0" \
            "pandas>=1.3" \
            "numpy>=1.20" \
            "scipy>=1.7" \
            "pybigwig>=0.3" \
            "pybedtools>=0.9" \
            "bedtools>=2.30" \
            "samtools>=1.15" \
            "macs2>=2.2" \
            "matplotlib>=3.5" \
            "seaborn>=0.11"

        echo -e "${GREEN}✓ Environment created manually${NC}"
    fi
}

# Install Python dependencies with pip (no-conda path). Best-effort: the
# authoritative pass/fail check is verify_python_deps() below, so a network
# failure here does not abort setup when the deps are already present.
setup_pip() {
    REQ="${SCRIPT_DIR}/requirements.txt"
    PIP="$(command -v pip3 || command -v pip || true)"

    if [ -z "$PIP" ]; then
        warn "pip not found - skipping automatic install; will verify what is present"
        return 0
    fi
    if [ ! -f "$REQ" ]; then
        warn "requirements.txt not found - skipping automatic install"
        return 0
    fi

    echo "Installing Python dependencies from requirements.txt (best-effort)..."
    LOG=/tmp/pace_pip_install.log
    if "$PIP" install -r "$REQ" >"$LOG" 2>&1 \
       || "$PIP" install --user -r "$REQ" >>"$LOG" 2>&1 \
       || "$PIP" install --break-system-packages -r "$REQ" >>"$LOG" 2>&1; then
        echo -e "${GREEN}✓ Python dependencies installed${NC}"
    else
        warn "Could not auto-install all packages (offline?). See $LOG"
        warn "Continuing - already-installed packages will still be verified."
    fi
}

# Verify Python dependencies are importable.
verify_python_deps() {
    echo "Verifying Python dependencies..."
    if [ -z "$PYTHON" ]; then
        err "Python interpreter not found"
        return 1
    fi

    # module:display-name. Core deps are required for the prediction pipeline.
    local core="yaml:PyYAML numpy:numpy pandas:pandas scipy:scipy"
    local optional="sklearn:scikit-learn pyBigWig:pyBigWig pysam:pysam matplotlib:matplotlib seaborn:seaborn"
    local missing_core=""

    local pair mod name
    for pair in $core; do
        mod="${pair%%:*}"; name="${pair##*:}"
        if "$PYTHON" -c "import $mod" 2>/dev/null; then
            ok "$name"
        else
            err "$name (required)"
            missing_core="$missing_core $name"
        fi
    done
    for pair in $optional; do
        mod="${pair%%:*}"; name="${pair##*:}"
        if "$PYTHON" -c "import $mod" 2>/dev/null; then
            ok "$name"
        else
            warn "$name (optional)"
        fi
    done

    if [ -z "$missing_core" ]; then
        echo -e "${GREEN}✓ All required Python dependencies present${NC}"
        return 0
    fi
    err "Missing required Python dependencies:$missing_core"
    if [ "$USE_CONDA" = false ]; then
        echo "  Install them with: pip install -r requirements.txt"
    fi
    return 1
}

# Verify external command-line tools.
verify_external_tools() {
    echo "Verifying external tools..."
    local hard_fail=0

    # bedtools is required for the dependency-light (no-conda) pipeline.
    if command -v bedtools &> /dev/null; then
        ok "bedtools ($(bedtools --version 2>/dev/null | head -1))"
    elif [ "$USE_CONDA" = true ]; then
        warn "bedtools (will be provided by the conda environment)"
    else
        err "bedtools - required for the direct pipeline"
        echo "    Install via: apt-get install bedtools  (or: conda install -c bioconda bedtools)"
        hard_fail=1
    fi

    # Optional tools: only needed for the full Snakemake/MACS2 workflow.
    local t
    for t in samtools macs2 snakemake tabix bgzip; do
        if command -v "$t" &> /dev/null; then
            ok "$t"
        else
            warn "$t (optional - needed for the full Snakemake/MACS2 workflow)"
        fi
    done

    [ "$hard_fail" -eq 0 ] && return 0 || return 1
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
            ok "${script}"
        else
            err "${script} - NOT FOUND"
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
    chmod +x "${SCRIPT_DIR}"/example/*.sh 2>/dev/null || true
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
        echo "  It will be generated automatically on first run, or manually with:"
        echo "    python scripts/generate_example_data.py"
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
    if [ "$USE_CONDA" = true ]; then
        echo "To activate the environment:"
        echo -e "  ${GREEN}conda activate pace-env${NC}"
        echo ""
    fi
    echo "To run the bundled example:"
    echo -e "  ${GREEN}bash example/run_example_direct.sh${NC}   # Python + bedtools only"
    if [ "$USE_CONDA" = true ]; then
        echo -e "  ${GREEN}bash example/run_example.sh${NC}          # full Snakemake/MACS2 workflow"
    fi
    echo ""
    echo "To verify the installation:"
    echo -e "  ${GREEN}python test_pace_complete.py${NC}"
    echo ""
    echo "To run with your own data:"
    echo "  1. Prepare reference files for your species"
    echo "  2. Edit config/config.yaml"
    echo "  3. Edit config/config_biosamples.tsv"
    if [ "$USE_CONDA" = true ]; then
        echo -e "  4. Run: ${GREEN}snakemake --cores 8${NC}"
    else
        echo "  4. Run the step scripts directly (see docs/TUTORIAL.md), or install"
        echo "     conda + snakemake for the full automated workflow."
    fi
    echo ""
    echo "For more information, see README.md and docs/TUTORIAL.md"
    echo "========================================"
}

# Main function
main() {
    echo ""

    # Step 1: Detect environment manager (conda vs pip)
    check_conda

    # Step 2: Set up dependencies
    if [ "$USE_CONDA" = true ]; then
        create_env_conda
    else
        setup_pip
    fi

    # Step 3: Verify Python dependencies (conda installs into pace-env, which is
    # not active in this shell, so only verify directly on the pip path).
    if [ "$USE_CONDA" = false ]; then
        verify_python_deps || exit 1
        verify_external_tools || exit 1
    fi

    # Step 4: Verify scripts
    verify_scripts || exit 1

    # Step 5: Make scripts executable
    make_executable

    # Step 6: Create directories
    create_dirs

    # Step 7: Verify example data
    verify_example

    # Print usage
    print_usage
}

# Run main
main
