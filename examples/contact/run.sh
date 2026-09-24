#!/usr/bin/env bash
# Sparse-contact workflow on synthetic tables; every input here is synthetic.
# Run from this folder:  bash run.sh
set -euo pipefail

CONTEXT="--species synthetic --assembly toy_assembly --tissue toy_tissue"
TABLES="--units units.tsv --promoters promoters.tsv --candidates candidates.tsv --samples samples.tsv --sources sources.tsv --evidence evidence.tsv --activity observed_activity.tsv --no-include-promoters --profile demonstration"

# 1. Candidate CTCF boundaries from motif hits
pace boundaries --motifs motifs.tsv -o boundary_asset

# 2a. A prior from known parameters ...
pace prior --a 1 --gamma 1 --beta 0 --kappa 2 --d-ref 1000 --d-min 1000 --boundaries boundary_asset/boundaries.tsv $CONTEXT --synthetic --model-id synthetic_boundary_prior --scale toy_contact --resolution 500 --normalization-id toy_mean --balancing unbalanced -o prior_asset

# 2b. ... or fitted on Hi-C counts, holding out chromosome c4
pace fit-hic --contacts observed_contacts.tsv --gamma-grid 0.5 1.0 1.5 --beta-grid 0 0.5 1 --test-chromosomes c4 --d-ref 1000 --d-min 1000 --boundaries boundary_asset/boundaries.tsv $CONTEXT --synthetic --model-id synthetic_boundary_prior -o fitted_prior

# 3. Score with per-pair Gamma-Poisson shrinkage of the sparse contacts
pace fuse $TABLES --contacts observed_contacts.tsv --contact-prior prior_asset $CONTEXT -o fused

# 4. Distance-only baseline run, then tune gamma/beta/eta on eQTL PIPs (c4 held out)
pace run $TABLES --contact-mode prior_only --contact-prior prior_asset $CONTEXT -o baseline
pace fit-labels --run baseline --eqtl eqtl.tsv --gamma-grid 1.0 --beta-grid 0 0.5 --eta-grid 0 0.5 1 --test-chromosomes c4 -o weak_fit

# 5. Apply the eQTL-tuned model
pace run $TABLES --contact-mode prior_only --contact-prior prior_asset --weak-model weak_fit/weak_model.json $CONTEXT -o calibrated
