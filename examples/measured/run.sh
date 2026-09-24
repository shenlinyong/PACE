#!/usr/bin/env bash
# Score the synthetic example tables in this folder.  Run from this folder: bash run.sh
set -euo pipefail

pace validate -d . --species synthetic --assembly toy_assembly --tissue toy_tissue \
  --profile demonstration --no-include-promoters
pace run -d . --species synthetic --assembly toy_assembly --tissue toy_tissue \
  --profile demonstration --no-include-promoters -o results
