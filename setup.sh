#!/usr/bin/env bash
# Author: shenlinyong — Linyong Shen, Northwest A&F University.
# Create a Conda environment for PACE.
set -euo pipefail

PACE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACE_RECIPE="$PACE_ROOT/environment.yml"
PACE_ENV_NAME=pace
PACE_PREFIX=

while [ "$#" -gt 0 ]; do
    case "$1" in
        --workflow)
            PACE_RECIPE="$PACE_ROOT/workflow/envs/pace-env.yml"
            PACE_ENV_NAME=pace-workflow
            shift ;;
        --prefix)
            if [ "$#" -lt 2 ] || [ -z "$2" ]; then
                echo '--prefix requires an environment path' >&2
                exit 2
            fi
            PACE_PREFIX="$2"
            shift 2 ;;
        -h|--help)
            echo 'Usage: bash setup.sh [--workflow] [--prefix /path/to/environment]'
            echo 'Default: create pace from environment.yml using Conda.'
            exit 0 ;;
        *)
            echo "Unknown option: $1" >&2
            exit 2 ;;
    esac
done

if ! command -v conda >/dev/null 2>&1; then
    echo 'Install Conda first; see docs/INSTALLATION.md.' >&2
    exit 1
fi

PACE_ARGS=(env create --file "$PACE_RECIPE")
if [ -n "$PACE_PREFIX" ]; then
    PACE_ARGS+=(--prefix "$PACE_PREFIX")
    PACE_ENV_NAME="$PACE_PREFIX"
fi
CONDA_CHANNEL_PRIORITY=strict conda "${PACE_ARGS[@]}"
printf '\nActivate with: conda activate %s\n' "$PACE_ENV_NAME"
echo 'Then run: python scripts/pace.py --help'
echo 'To verify: python -m pytest tests -q'
