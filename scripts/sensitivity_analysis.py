#!/usr/bin/env python3
"""
PACE: Parameter Sensitivity and Top-Prediction Stability Analysis

Systematically perturb PACE signal weights (activating AND inhibitory, i.e.
DNA methylation / H3K27me3 / H3K9me3) and quantify how much the predictions
change. For each weight we report, in addition to global AUROC:

  * Jaccard similarity of the top-ranked E-P predictions (default top 10%);
  * overlap percentage of the top-ranked predictions;
  * Spearman correlation of the full score vectors;

between the perturbed and the default configuration. Inhibitory weights are
included in the sensitivity analysis, demonstrating that the biological
conclusions (the set of top-ranked pairs) — not only the global AUROC curve —
remain stable when weights are varied by +/- 50%.

The analysis recomputes the PACE score directly from per-enhancer signal
columns, so it requires a predictions/enhancer table that contains the
individual signal tracks (e.g. ATAC_signal, H3K27ac_signal, ...), the
contact value, and the target-gene/distance columns. The standard PACE
``EnhancerPredictionsAllPutative.tsv.gz`` (joined with EnhancerList signals)
or any table with these columns can be used.

Author: Linyong Shen @ Northwest A&F University
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd


# Default weights matching config/config.yaml.
DEFAULT_ACTIVATING = {
    "ATAC": 1.5, "H3K27ac": 1.0, "H3K4me1": 0.8, "H3K4me3": 0.5,
}
DEFAULT_INHIBITORY = {
    "methylation": 0.5, "H3K27me3": 0.5, "H3K9me3": 0.5,
}


def _find_col(df: pd.DataFrame, signal: str):
    """Find the column for a signal (e.g. 'ATAC' -> 'ATAC_signal')."""
    for cand in (signal, f"{signal}_signal", signal.lower(),
                 f"{signal.lower()}_signal"):
        if cand in df.columns:
            return cand
    return None


def compute_activity(df: pd.DataFrame,
                     act_weights: dict,
                     inh_weights: dict,
                     act_cols: dict,
                     inh_cols: dict) -> np.ndarray:
    """Normalized weighted-geometric activity with inhibitory down-weighting."""
    n = len(df)
    log_sum = np.zeros(n)
    w_sum = 0.0
    for sig, w in act_weights.items():
        col = act_cols.get(sig)
        if col is None:
            continue
        vals = np.clip(df[col].to_numpy(dtype=float), 1e-10, None)
        log_sum += w * np.log(vals)
        w_sum += w
    w_sum = w_sum if w_sum > 0 else 1.0
    activity = np.exp(log_sum / w_sum)

    inhibitory = np.zeros(n)
    for sig, w in inh_weights.items():
        col = inh_cols.get(sig)
        if col is None:
            continue
        v = df[col].to_numpy(dtype=float)
        rng = v.max() - v.min()
        v = (v - v.min()) / rng if rng > 0 else np.zeros_like(v)
        inhibitory += w * v
    inhibitory = np.clip(inhibitory, 0, 1)
    return activity * (1 - inhibitory)


def compute_pace_score(df: pd.DataFrame, activity: np.ndarray,
                       contact_col: str, gene_col: str,
                       expr_col=None) -> np.ndarray:
    """ABC-style normalized score: A*C*W / sum_per_gene(A*C*W)."""
    axc = activity * df[contact_col].to_numpy(dtype=float)
    if expr_col and expr_col in df.columns:
        axc = axc * df[expr_col].to_numpy(dtype=float)
    tmp = pd.DataFrame({"gene": df[gene_col].values, "axc": axc})
    gene_sum = tmp.groupby("gene")["axc"].transform("sum").to_numpy()
    with np.errstate(divide="ignore", invalid="ignore"):
        score = np.where(gene_sum > 0, axc / gene_sum, 0.0)
    return score


def top_k_set(scores: np.ndarray, frac: float) -> set:
    """Indices of the top ``frac`` fraction of predictions by score."""
    k = max(1, int(np.ceil(len(scores) * frac)))
    return set(np.argsort(-scores)[:k].tolist())


def jaccard(a: set, b: set) -> float:
    if not a and not b:
        return 1.0
    return len(a & b) / len(a | b)


def run_sensitivity(df, act_weights, inh_weights, act_cols, inh_cols,
                    contact_col, gene_col, expr_col, top_frac,
                    factors, validated_col=None):
    """Perturb each weight by the given factors and measure stability."""
    from scipy.stats import spearmanr

    base_activity = compute_activity(df, act_weights, inh_weights,
                                     act_cols, inh_cols)
    base_score = compute_pace_score(df, base_activity, contact_col,
                                    gene_col, expr_col)
    base_top = top_k_set(base_score, top_frac)

    auroc_fn = None
    if validated_col and validated_col in df.columns:
        from sklearn.metrics import roc_auc_score
        y = df[validated_col].to_numpy()
        if len(np.unique(y)) == 2:
            auroc_fn = lambda s: roc_auc_score(y, s)

    rows = []
    all_weights = [(s, w, "activating") for s, w in act_weights.items()] + \
                  [(s, w, "inhibitory") for s, w in inh_weights.items()]

    for sig, base_w, kind in all_weights:
        col = (act_cols if kind == "activating" else inh_cols).get(sig)
        if col is None:
            continue  # signal not present in this dataset
        for factor in factors:
            aw = dict(act_weights)
            iw = dict(inh_weights)
            if kind == "activating":
                aw[sig] = base_w * factor
            else:
                iw[sig] = base_w * factor
            activity = compute_activity(df, aw, iw, act_cols, inh_cols)
            score = compute_pace_score(df, activity, contact_col,
                                      gene_col, expr_col)
            top = top_k_set(score, top_frac)
            rho = spearmanr(base_score, score).correlation
            row = {
                "parameter": sig, "type": kind,
                "factor": factor, "weight": base_w * factor,
                "jaccard_top": jaccard(base_top, top),
                "overlap_pct_top": 100.0 * len(base_top & top) / max(1, len(base_top)),
                "spearman_full": rho,
            }
            if auroc_fn is not None:
                row["auroc"] = auroc_fn(score)
            rows.append(row)
    return pd.DataFrame(rows), base_score


def main():
    parser = argparse.ArgumentParser(
        description="PACE parameter sensitivity & top-prediction stability",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--predictions", required=True,
                        help="Predictions/enhancer table with signal columns")
    parser.add_argument("--output", required=True, help="Output TSV")
    parser.add_argument("--gene_col", default="TargetGene")
    parser.add_argument("--contact_col", default="contact")
    parser.add_argument("--expr_col", default="ExpressionWeight",
                        help="Expression weight column (optional)")
    parser.add_argument("--validated_col", default=None,
                        help="Binary validation column for AUROC (optional)")
    parser.add_argument("--top_frac", type=float, default=0.10,
                        help="Top fraction for Jaccard/overlap (default 0.10)")
    parser.add_argument("--range", type=float, default=0.5,
                        help="Perturbation range as fraction (default 0.5 = +/-50%%)")
    parser.add_argument("--n_steps", type=int, default=5,
                        help="Number of factor steps across the range")
    parser.add_argument("--plot", default=None, help="Optional output plot")
    args = parser.parse_args()

    if not os.path.exists(args.predictions):
        sys.exit(f"[PACE] Predictions not found: {args.predictions}")

    df = pd.read_csv(args.predictions, sep="\t",
                     compression="gzip" if args.predictions.endswith(".gz") else None)

    act_cols = {s: _find_col(df, s) for s in DEFAULT_ACTIVATING}
    inh_cols = {s: _find_col(df, s) for s in DEFAULT_INHIBITORY}
    present_act = {s: w for s, w in DEFAULT_ACTIVATING.items() if act_cols[s]}
    present_inh = {s: w for s, w in DEFAULT_INHIBITORY.items() if inh_cols[s]}

    if not present_act:
        sys.exit("[PACE] No activating signal columns found "
                 "(expected e.g. ATAC_signal, H3K27ac_signal).")
    if args.contact_col not in df.columns:
        sys.exit(f"[PACE] Contact column '{args.contact_col}' not found.")

    print(f"[PACE] Activating signals: {list(present_act)}")
    print(f"[PACE] Inhibitory signals: {list(present_inh)}")

    factors = np.linspace(1 - args.range, 1 + args.range, args.n_steps)
    expr_col = args.expr_col if args.expr_col in df.columns else None

    results, _ = run_sensitivity(
        df, present_act, present_inh, act_cols, inh_cols,
        args.contact_col, args.gene_col, expr_col, args.top_frac,
        factors, validated_col=args.validated_col)

    os.makedirs(os.path.dirname(os.path.abspath(args.output)) or ".",
                exist_ok=True)
    results.to_csv(args.output, sep="\t", index=False)
    print(f"[PACE] Wrote sensitivity results to {args.output}")

    # Stability summary over the +/-range perturbation.
    print("\nStability summary (top-{:.0%} predictions):".format(args.top_frac))
    summary = (results.groupby("parameter")
               .agg(min_jaccard=("jaccard_top", "min"),
                    mean_jaccard=("jaccard_top", "mean"),
                    min_spearman=("spearman_full", "min"))
               .sort_values("min_jaccard"))
    print(summary.to_string())

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 5))
        for sig, grp in results.groupby("parameter"):
            grp = grp.sort_values("factor")
            ax.plot(grp["factor"], grp["jaccard_top"], marker="o", label=sig)
        ax.axvspan(1 - args.range, 1 + args.range, color="grey", alpha=0.08)
        ax.set_xlabel("Weight multiplier")
        ax.set_ylabel(f"Jaccard of top {args.top_frac:.0%} predictions")
        ax.set_title("PACE top-prediction stability under weight perturbation")
        ax.set_ylim(0, 1.02)
        ax.legend(fontsize=8, ncol=2)
        plt.tight_layout()
        plt.savefig(args.plot, dpi=200, bbox_inches="tight")
        print(f"[PACE] Wrote plot to {args.plot}")


if __name__ == "__main__":
    main()
