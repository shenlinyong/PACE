#!/usr/bin/env python3
"""
PACE: DNA Methylation (WGBS/RRBS) Processing

Convert per-cytosine methylation calls into the region-level fractional
methylation track that PACE uses as an inhibitory signal.

Supported input formats
-----------------------
1. ``bismark_cov``   Bismark coverage file (``*.bismark.cov[.gz]``):
                     chrom  start  end  methylation_percentage  count_M  count_U
2. ``bismark_report``Bismark genome-wide cytosine report
                     (``*CpG_report.txt[.gz]``):
                     chrom  pos  strand  count_M  count_U  context  tri-context
3. ``bedgraph``      Generic bedGraph: chrom  start  end  value (0-1 or 0-100)

Processing steps
----------------
* keep cytosines in the requested context (default CpG / "CG");
* require a minimum read coverage per cytosine (default >= 4x), discarding
  low-confidence calls;
* compute fractional methylation beta = M / (M + U) in [0, 1];
* optionally aggregate to regions (mean beta over each region), which is the
  representation consumed by ``scripts/multiomics_activity.py`` and the
  Snakemake ``methylation`` column.

Outputs a bedGraph (chrom start end beta) of per-CpG methylation, and — when
``--regions`` is supplied — a per-region mean-methylation bedGraph.

Author: Linyong Shen @ Northwest A&F University
"""

import argparse
import gzip
import os
import sys

import numpy as np
import pandas as pd


def _open(path, mode="rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def load_bismark_cov(path: str) -> pd.DataFrame:
    """Load a Bismark coverage file into a per-cytosine DataFrame."""
    df = pd.read_csv(path, sep="\t", header=None,
                     names=["chr", "start", "end", "pct", "M", "U"])
    df["pos"] = df["start"]
    df["coverage"] = df["M"] + df["U"]
    df["context"] = "CG"  # coverage files are typically already CpG-only
    return df[["chr", "pos", "M", "U", "coverage", "context"]]


def load_bismark_report(path: str) -> pd.DataFrame:
    """Load a Bismark genome-wide cytosine report."""
    df = pd.read_csv(path, sep="\t", header=None,
                     names=["chr", "pos", "strand", "M", "U",
                            "context", "tricontext"])
    df["coverage"] = df["M"] + df["U"]
    return df[["chr", "pos", "M", "U", "coverage", "context"]]


def load_bedgraph(path: str) -> pd.DataFrame:
    """Load a generic methylation bedGraph (chrom start end value)."""
    df = pd.read_csv(path, sep="\t", header=None,
                     names=["chr", "start", "end", "value"])
    df["pos"] = df["start"]
    # Coverage is unknown for bedGraph input; treat each record as one call.
    df["coverage"] = np.nan
    val = df["value"].astype(float)
    if val.max() > 1.5:  # percentage -> fraction
        val = val / 100.0
    df["beta"] = val.clip(0, 1)
    df["context"] = "CG"
    return df[["chr", "pos", "beta", "coverage", "context"]]


def per_cytosine_beta(df: pd.DataFrame, context: str,
                      min_coverage: int) -> pd.DataFrame:
    """Filter by context/coverage and compute beta = M / (M + U)."""
    if "beta" in df.columns:  # bedGraph already has beta
        out = df.copy()
    else:
        ctx = context.upper()
        if "context" in df.columns and ctx not in ("ALL", "ANY"):
            df = df[df["context"].str.upper() == ctx]
        df = df[df["coverage"] >= min_coverage].copy()
        df["beta"] = df["M"] / df["coverage"].replace(0, np.nan)
        out = df.dropna(subset=["beta"])
    return out[["chr", "pos", "beta"]].sort_values(["chr", "pos"])


def aggregate_over_regions(cpg: pd.DataFrame, regions: pd.DataFrame) -> pd.DataFrame:
    """Mean beta over each region using a simple interval join."""
    regions = regions.copy().reset_index(drop=True)
    if "name" not in regions.columns:
        regions["name"] = [f"region_{i}" for i in range(len(regions))]
    out = []
    for chrom, grp in cpg.groupby("chr"):
        rsub = regions[regions["chr"] == chrom]
        if len(rsub) == 0:
            continue
        positions = grp["pos"].to_numpy()
        betas = grp["beta"].to_numpy()
        order = positions.argsort()
        positions, betas = positions[order], betas[order]
        for _, r in rsub.iterrows():
            lo = np.searchsorted(positions, r["start"], side="left")
            hi = np.searchsorted(positions, r["end"], side="right")
            mean_beta = betas[lo:hi].mean() if hi > lo else 0.0
            out.append({"chr": chrom, "start": int(r["start"]),
                        "end": int(r["end"]), "name": r["name"],
                        "methylation": float(mean_beta),
                        "n_cpg": int(hi - lo)})
    return pd.DataFrame(out)


def main():
    parser = argparse.ArgumentParser(
        description="Process WGBS/RRBS methylation calls for PACE",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Bismark coverage -> per-CpG methylation bedGraph
  python scripts/process_methylation.py \\
      --input sample.bismark.cov.gz --format bismark_cov \\
      --output sample.CpG.methylation.bedGraph

  # Aggregate to candidate enhancer regions (for PACE inhibitory signal)
  python scripts/process_methylation.py \\
      --input sample.CpG_report.txt.gz --format bismark_report \\
      --regions results/Sample/Peaks/candidateRegions.bed \\
      --output sample.region_methylation.bedGraph
        """)
    parser.add_argument("--input", required=True, help="Methylation input file")
    parser.add_argument("--format", default="bismark_cov",
                        choices=["bismark_cov", "bismark_report", "bedgraph"],
                        help="Input format")
    parser.add_argument("--output", required=True, help="Output bedGraph")
    parser.add_argument("--context", default="CG",
                        help="Cytosine context to keep (CG/CHG/CHH/ALL)")
    parser.add_argument("--min_coverage", type=int, default=4,
                        help="Minimum read coverage per cytosine (default: 4)")
    parser.add_argument("--regions", default=None,
                        help="Optional BED of regions to aggregate over")
    args = parser.parse_args()

    if not os.path.exists(args.input):
        sys.exit(f"[PACE] Input not found: {args.input}")

    loaders = {"bismark_cov": load_bismark_cov,
               "bismark_report": load_bismark_report,
               "bedgraph": load_bedgraph}
    print(f"[PACE] Loading {args.format} from {args.input}")
    raw = loaders[args.format](args.input)

    cpg = per_cytosine_beta(raw, args.context, args.min_coverage)
    print(f"[PACE] Retained {len(cpg):,} cytosines "
          f"(context={args.context}, min_cov={args.min_coverage})")

    os.makedirs(os.path.dirname(os.path.abspath(args.output)) or ".",
                exist_ok=True)

    if args.regions:
        regions = pd.read_csv(args.regions, sep="\t", header=None,
                              usecols=[0, 1, 2, 3],
                              names=["chr", "start", "end", "name"])
        out = aggregate_over_regions(cpg, regions)
        out.to_csv(args.output, sep="\t", index=False, header=False,
                   columns=["chr", "start", "end", "methylation"])
        print(f"[PACE] Wrote region methylation for {len(out):,} regions "
              f"to {args.output}")
    else:
        cpg_out = cpg.copy()
        cpg_out["end"] = cpg_out["pos"] + 1
        cpg_out[["chr", "pos", "end", "beta"]].to_csv(
            args.output, sep="\t", index=False, header=False)
        print(f"[PACE] Wrote per-CpG methylation to {args.output}")


if __name__ == "__main__":
    main()
