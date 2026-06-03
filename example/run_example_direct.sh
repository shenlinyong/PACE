#!/bin/bash
# PACE Example Runner (direct, no Snakemake/MACS2 required)
#
# Runs the complete PACE pipeline on the bundled synthetic dataset by calling
# the five step scripts directly. This is the recommended way to try PACE in
# a minimal environment (Python + bedtools only). For production runs on real
# data, use the Snakemake workflow (see README / docs/TUTORIAL.md).
#
# Usage: bash example/run_example_direct.sh
#
# Author: Linyong Shen @ Northwest A&F University

set -euo pipefail

# Resolve repository root (parent of this script's directory).
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$SCRIPT_DIR")"
cd "$ROOT"

S=workflow/scripts
DATA=example/data
REF=example/reference
OUT=example/results/Example_Sample
mkdir -p "$OUT"/Peaks "$OUT"/Neighborhoods "$OUT"/Predictions "$OUT"/Metrics

echo "============================================================"
echo "PACE example pipeline (direct mode)"
echo "============================================================"

# Step 0: generate synthetic data + reference if missing.
if [ ! -f "$DATA/example_peaks.narrowPeak" ] || [ ! -f "$DATA/example_ATAC.tagAlign.gz" ]; then
    echo "[0/5] Generating synthetic example data ..."
    python scripts/generate_example_data.py
fi
if [ ! -f "$REF/genes.bed" ]; then
    echo "[0/5] Preparing reference files ..."
    python scripts/prepare_reference.py \
        --gtf "$REF/example_annotation.gtf.gz" \
        --fasta "$REF/example_genome.fa.gz" \
        --output_dir "$REF" --species_name example
fi

PEAKS="$DATA/example_peaks.narrowPeak"
[ -f "$PEAKS" ] || PEAKS="$OUT/Peaks/macs2_peaks.narrowPeak"

echo "[1/5] Candidate regions ..."
python $S/pace_candidate_regions.py \
    --narrowPeak "$PEAKS" \
    --chrom_sizes "$REF/example.chrom.sizes" \
    --output "$OUT/Peaks/candidateRegions.bed" \
    --nStrongestPeaks 150000 --peakExtendFromSummit 250

echo "[2/5] Neighborhood analysis (ATAC + H3K27ac, weighted geometric) ..."
python $S/pace_neighborhoods.py \
    --candidate_regions "$OUT/Peaks/candidateRegions.bed" \
    --genes "$REF/genes.bed" \
    --chrom_sizes "$REF/example.chrom.sizes" \
    --output_dir "$OUT/Neighborhoods" \
    --accessibility_file "$DATA/example_ATAC.tagAlign.gz" \
    --H3K27ac "$DATA/example_H3K27ac.tagAlign.gz" \
    --activity_method weighted_geometric

echo "[3/5] Predicting enhancer-gene interactions ..."
python $S/pace_predict.py \
    --enhancers "$OUT/Neighborhoods/EnhancerList.txt" \
    --genes "$OUT/Neighborhoods/GeneList.txt" \
    --output "$OUT/Predictions/EnhancerPredictionsAllPutative.tsv.gz" \
    --threshold 0.02

echo "[4/5] Filtering predictions ..."
python $S/pace_filter.py \
    --predictions "$OUT/Predictions/EnhancerPredictionsAllPutative.tsv.gz" \
    --output "$OUT/Predictions/EnhancerPredictions.tsv" --threshold 0.02

echo "[5/5] Quality-control metrics ..."
python $S/pace_metrics.py \
    --predictions "$OUT/Predictions/EnhancerPredictionsAllPutative.tsv.gz" \
    --output_dir "$OUT/Metrics" --sample_name Example_Sample

echo ""
echo "============================================================"
echo "Done. Results in $OUT/"
echo "------------------------------------------------------------"
N=$(zcat "$OUT/Predictions/EnhancerPredictionsAllPutative.tsv.gz" | tail -n +2 | wc -l)
echo "Total E-G predictions: $N"
echo "Top predictions:"
# awk (rather than head) reads the whole stream, avoiding a spurious SIGPIPE
# non-zero exit under `set -o pipefail`.
zcat "$OUT/Predictions/EnhancerPredictionsAllPutative.tsv.gz" \
    | awk -F'\t' 'NR<=4{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}'
echo "============================================================"
