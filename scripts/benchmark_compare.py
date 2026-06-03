#!/usr/bin/env python3
"""
PACE: Benchmarking Harness for Comparing E-P Prediction Methods

Compare PACE against any number of other enhancer-promoter prediction methods
on a common set of experimentally validated E-P pairs (e.g. CRISPRi), using a
single, transparent evaluation protocol. The harness is fully reproducible and
lets users add independent algorithms such as EPIPDLF and GATv2EPI simply by
supplying their score files.

Each method is provided as ``NAME:FILE:SCORE_COLUMN``. The file must contain
``enhancer`` and ``gene`` identifier columns (or chr/start/end + TargetGene),
plus the named score column. All methods are intersected to the same set of
evaluable pairs so the comparison is apples-to-apples; for each method the
script reports AUROC, AUPRC, and sensitivity at fixed specificity, and draws
overlaid ROC and precision-recall curves.

Usage
-----
python scripts/benchmark_compare.py \
    --labels crispri_validated.tsv \
    --method PACE:pace_pred.tsv:PACE.Score \
    --method ABC:abc_pred.tsv:ABC.Score \
    --method EPIPDLF:epipdlf_scores.tsv:score \
    --method GATv2EPI:gatv2_scores.tsv:score \
    --output benchmark/

Author: Linyong Shen @ Northwest A&F University
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd


def _pair_id(df: pd.DataFrame) -> pd.Series:
    """Build a canonical enhancer||gene identifier for joining methods."""
    if "enhancer" in df.columns and "gene" in df.columns:
        enh, gene = df["enhancer"].astype(str), df["gene"].astype(str)
    else:
        gene_col = "gene" if "gene" in df.columns else (
            "TargetGene" if "TargetGene" in df.columns else None)
        if gene_col is None or not {"chr", "start", "end"}.issubset(df.columns):
            raise ValueError(
                "Need 'enhancer'+'gene' or 'chr/start/end'+'TargetGene' columns")
        enh = (df["chr"].astype(str) + ":" + df["start"].astype(str)
               + "-" + df["end"].astype(str))
        gene = df[gene_col].astype(str)
    return enh + "||" + gene


def load_method(spec: str) -> pd.DataFrame:
    """Parse ``NAME:FILE:SCORE_COLUMN`` and return [pair, NAME]."""
    parts = spec.split(":")
    if len(parts) < 3:
        raise ValueError(f"Bad --method spec '{spec}'. "
                         "Use NAME:FILE:SCORE_COLUMN")
    name, path, score_col = parts[0], ":".join(parts[1:-1]), parts[-1]
    comp = "gzip" if path.endswith(".gz") else None
    df = pd.read_csv(path, sep="\t", compression=comp)
    if score_col not in df.columns:
        raise ValueError(f"Score column '{score_col}' not in {path}")
    out = pd.DataFrame({"pair": _pair_id(df), name: df[score_col]})
    # Keep the strongest score per pair (in case of duplicates).
    return out.groupby("pair", as_index=False)[name].max()


def evaluate(y, scores):
    """Return AUROC, AUPRC, sensitivity@80%spec, sensitivity@90%spec."""
    from sklearn.metrics import roc_auc_score, average_precision_score, roc_curve
    res = {
        "auroc": roc_auc_score(y, scores),
        "auprc": average_precision_score(y, scores),
        "n_pairs": len(y), "n_positive": int(np.sum(y)),
    }
    fpr, tpr, _ = roc_curve(y, scores)
    for spec in (0.80, 0.90):
        idx = np.where(fpr <= (1 - spec))[0]
        res[f"sens_at_{int(spec*100)}spec"] = float(tpr[idx[-1]]) if len(idx) else 0.0
    return res


def main():
    parser = argparse.ArgumentParser(
        description="Benchmark E-P prediction methods on a common gold set")
    parser.add_argument("--labels", required=True,
                        help="Validated pairs: enhancer,gene,label (1/0)")
    parser.add_argument("--method", action="append", required=True,
                        help="NAME:FILE:SCORE_COLUMN (repeatable)")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--label_col", default="label")
    args = parser.parse_args()

    if not os.path.exists(args.labels):
        sys.exit(f"[PACE] Labels file not found: {args.labels}")
    os.makedirs(args.output, exist_ok=True)

    labels = pd.read_csv(args.labels, sep="\t")
    if args.label_col not in labels.columns:
        sys.exit(f"[PACE] Label column '{args.label_col}' not found.")
    labels = labels.assign(pair=_pair_id(labels))[["pair", args.label_col]]
    labels = labels.rename(columns={args.label_col: "label"})
    labels = labels.groupby("pair", as_index=False)["label"].max()

    merged = labels
    method_names = []
    for spec in args.method:
        mdf = load_method(spec)
        name = mdf.columns[1]
        method_names.append(name)
        merged = merged.merge(mdf, on="pair", how="left")

    # Restrict to pairs scored by all methods so the comparison is fair.
    before = len(merged)
    merged = merged.dropna(subset=method_names)
    print(f"[PACE] Evaluable pairs scored by all methods: {len(merged)} "
          f"(of {before} labelled)")
    if len(merged) == 0 or merged["label"].nunique() < 2:
        sys.exit("[PACE] Not enough overlapping/2-class pairs to benchmark.")

    y = merged["label"].to_numpy()
    rows = []
    for name in method_names:
        r = evaluate(y, merged[name].to_numpy())
        r["method"] = name
        rows.append(r)
    summary = pd.DataFrame(rows).set_index("method").sort_values(
        "auroc", ascending=False)
    summary_file = os.path.join(args.output, "benchmark_summary.tsv")
    summary.to_csv(summary_file, sep="\t")
    print(f"[PACE] Wrote summary to {summary_file}")
    print(summary[["auroc", "auprc", "sens_at_80spec"]].to_string())

    # ROC and PR plots.
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from sklearn.metrics import roc_curve, precision_recall_curve

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        for name in method_names:
            s = merged[name].to_numpy()
            fpr, tpr, _ = roc_curve(y, s)
            ax1.plot(fpr, tpr, label=f"{name} (AUROC={summary.loc[name,'auroc']:.2f})")
            prec, rec, _ = precision_recall_curve(y, s)
            ax2.plot(rec, prec, label=f"{name} (AUPRC={summary.loc[name,'auprc']:.2f})")
        ax1.plot([0, 1], [0, 1], "k--", alpha=0.4)
        ax1.set_xlabel("False positive rate"); ax1.set_ylabel("True positive rate")
        ax1.set_title("ROC"); ax1.legend(fontsize=8)
        ax2.set_xlabel("Recall"); ax2.set_ylabel("Precision")
        ax2.set_title("Precision-Recall"); ax2.legend(fontsize=8)
        plt.tight_layout()
        plot_file = os.path.join(args.output, "benchmark_curves.png")
        plt.savefig(plot_file, dpi=200, bbox_inches="tight")
        print(f"[PACE] Wrote curves to {plot_file}")
    except Exception as exc:  # pragma: no cover
        print(f"[PACE] Plotting skipped: {exc}")


if __name__ == "__main__":
    main()
