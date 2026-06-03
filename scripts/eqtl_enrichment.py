#!/usr/bin/env python3
"""
PACE: eQTL Enrichment Analysis with Full Contingency Tables

Quantify the enrichment of significant cis-eQTL variants in PACE-predicted
enhancer regions relative to size-matched control regions, using Fisher's
exact test. For every comparison the script reports the odds ratio, its 95%
confidence interval, the p-value AND the four cell counts of the underlying
2x2 contingency table, so the results can be tabulated directly.

It additionally stratifies predictions by PACE-score percentile (top 10% /
top 25% / above-median / below-median) and reports the fold enrichment in
each stratum, giving a score dose-response analysis.

Inputs
------
--predictions  PACE predictions (chr, start, end, score column, optionally
               TargetGene). Plain TSV or .gz.
--eqtl         Significant cis-eQTL variants (chr, pos[, gene, pvalue]).
--chrom_sizes  Chromosome sizes file (for size-matched control regions).

Author: Linyong Shen @ Northwest A&F University
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd


def load_table(path, **kw):
    comp = "gzip" if str(path).endswith(".gz") else None
    return pd.read_csv(path, sep="\t", compression=comp, **kw)


def odds_ratio_ci(a, b, c, d, alpha=0.05):
    """Odds ratio and Wald 95% CI with Haldane-Anscombe correction."""
    from scipy.stats import norm
    a_, b_, c_, d_ = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    or_ = (a_ * d_) / (b_ * c_)
    se = np.sqrt(1 / a_ + 1 / b_ + 1 / c_ + 1 / d_)
    z = norm.ppf(1 - alpha / 2)
    return or_, or_ * np.exp(-z * se), or_ * np.exp(z * se)


def variants_in_regions(regions: pd.DataFrame,
                        variants: pd.DataFrame) -> np.ndarray:
    """Boolean mask: does each region contain >=1 variant?

    Uses a per-chromosome sorted search; O((R+V) log V).
    """
    has = np.zeros(len(regions), dtype=bool)
    var_by_chr = {c: np.sort(g["pos"].to_numpy())
                  for c, g in variants.groupby("chr")}
    regions = regions.reset_index(drop=True)
    for i, r in regions.iterrows():
        positions = var_by_chr.get(r["chr"])
        if positions is None:
            continue
        lo = np.searchsorted(positions, r["start"], side="left")
        hi = np.searchsorted(positions, r["end"], side="right")
        has[i] = hi > lo
    return has


def make_control_regions(regions: pd.DataFrame, chrom_sizes: dict,
                        seed: int = 42) -> pd.DataFrame:
    """Size-matched control regions: same chromosome and width, random start."""
    rng = np.random.default_rng(seed)
    rows = []
    for _, r in regions.iterrows():
        size = chrom_sizes.get(r["chr"])
        width = int(r["end"] - r["start"])
        if not size or width <= 0 or size <= width:
            continue
        start = int(rng.integers(0, size - width))
        rows.append({"chr": r["chr"], "start": start, "end": start + width})
    return pd.DataFrame(rows)


def fisher_table(pred_has, ctrl_has):
    """Build the 2x2 table and run Fisher's exact test (greater)."""
    from scipy.stats import fisher_exact
    a = int(pred_has.sum())                 # predicted, with eQTL
    b = int((~pred_has).sum())              # predicted, without
    c = int(ctrl_has.sum())                 # control, with
    d = int((~ctrl_has).sum())              # control, without
    or_fisher, pval = fisher_exact([[a, b], [c, d]], alternative="greater")
    or_wald, lo, hi = odds_ratio_ci(a, b, c, d)
    return {
        "pred_with_eqtl": a, "pred_without_eqtl": b,
        "ctrl_with_eqtl": c, "ctrl_without_eqtl": d,
        "odds_ratio": or_wald, "or_ci_low": lo, "or_ci_high": hi,
        "fisher_or": or_fisher, "pvalue": pval,
        "pred_frac": a / (a + b) if (a + b) else 0,
        "ctrl_frac": c / (c + d) if (c + d) else 0,
    }


def main():
    parser = argparse.ArgumentParser(
        description="eQTL enrichment with full contingency tables")
    parser.add_argument("--predictions", required=True)
    parser.add_argument("--eqtl", required=True,
                        help="Significant eQTL variants (chr, pos[, pvalue])")
    parser.add_argument("--chrom_sizes", required=True)
    parser.add_argument("--output", required=True, help="Output TSV (tables)")
    parser.add_argument("--score_col", default="ABC.Score")
    parser.add_argument("--threshold", type=float, default=0.02,
                        help="Score threshold for the primary comparison")
    parser.add_argument("--pvalue_col", default=None,
                        help="eQTL p-value column (filter to significant)")
    parser.add_argument("--pvalue_threshold", type=float, default=1e-5)
    parser.add_argument("--label", default="all",
                        help="Label (e.g. species_tissue) for the output rows")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    for p in (args.predictions, args.eqtl, args.chrom_sizes):
        if not os.path.exists(p):
            sys.exit(f"[PACE] File not found: {p}")

    preds = load_table(args.predictions)
    if args.score_col not in preds.columns:
        sys.exit(f"[PACE] Score column '{args.score_col}' not found.")

    eqtl = load_table(args.eqtl)
    # Normalize eQTL column names.
    rename = {}
    for c in eqtl.columns:
        cl = c.lower()
        if cl in ("chr", "chrom", "chromosome", "#chr"):
            rename[c] = "chr"
        elif cl in ("pos", "position", "bp", "snp_pos", "variant_pos"):
            rename[c] = "pos"
    eqtl = eqtl.rename(columns=rename)
    if "chr" not in eqtl.columns or "pos" not in eqtl.columns:
        sys.exit("[PACE] eQTL file must contain chromosome and position columns.")
    if args.pvalue_col and args.pvalue_col in eqtl.columns:
        eqtl = eqtl[eqtl[args.pvalue_col] <= args.pvalue_threshold]
    eqtl["chr"] = eqtl["chr"].astype(str)
    preds["chr"] = preds["chr"].astype(str)

    chrom_sizes = {}
    with open(args.chrom_sizes) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) >= 2:
                chrom_sizes[parts[0]] = int(parts[1])

    print(f"[PACE] {len(preds):,} predictions, {len(eqtl):,} eQTL variants")

    # --- Primary comparison: predicted enhancers vs size-matched controls ----
    sig_preds = preds[preds[args.score_col] >= args.threshold]
    enh = sig_preds[["chr", "start", "end"]].drop_duplicates()
    ctrl = make_control_regions(enh, chrom_sizes, seed=args.seed)

    pred_has = variants_in_regions(enh, eqtl)
    ctrl_has = variants_in_regions(ctrl, eqtl)
    primary = fisher_table(pred_has, ctrl_has)
    primary.update({"label": args.label, "comparison": "enhancers_vs_control",
                    "score_threshold": args.threshold,
                    "n_enhancers": len(enh)})

    results = [primary]
    print(f"[PACE] Enhancers vs control: OR={primary['odds_ratio']:.2f} "
          f"(95% CI {primary['or_ci_low']:.2f}-{primary['or_ci_high']:.2f}), "
          f"P={primary['pvalue']:.2e}")

    # --- Score-stratified fold enrichment (Fig. 5B) -------------------------
    all_enh = preds[["chr", "start", "end", args.score_col]].drop_duplicates(
        subset=["chr", "start", "end"])
    all_has = variants_in_regions(all_enh[["chr", "start", "end"]], eqtl)
    all_enh = all_enh.assign(has_eqtl=all_has)
    ctrl_all = make_control_regions(all_enh, chrom_sizes, seed=args.seed)
    ctrl_frac = variants_in_regions(ctrl_all, eqtl).mean()
    ctrl_frac = ctrl_frac if ctrl_frac > 0 else 1e-9

    scores = all_enh[args.score_col]
    strata = {
        "top_10pct": scores >= scores.quantile(0.90),
        "top_25pct": scores >= scores.quantile(0.75),
        "above_median": scores >= scores.quantile(0.50),
        "below_median": scores < scores.quantile(0.50),
    }
    for name, mask in strata.items():
        frac = all_enh.loc[mask, "has_eqtl"].mean() if mask.any() else 0.0
        results.append({
            "label": args.label, "comparison": f"stratum_{name}",
            "n_enhancers": int(mask.sum()),
            "pred_frac": frac, "ctrl_frac": ctrl_frac,
            "fold_enrichment": frac / ctrl_frac,
        })
        print(f"[PACE] {name}: fold enrichment = {frac / ctrl_frac:.2f}")

    out = pd.DataFrame(results)
    os.makedirs(os.path.dirname(os.path.abspath(args.output)) or ".",
                exist_ok=True)
    out.to_csv(args.output, sep="\t", index=False)
    print(f"[PACE] Wrote enrichment tables to {args.output}")


if __name__ == "__main__":
    main()
