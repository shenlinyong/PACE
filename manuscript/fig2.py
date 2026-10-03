"""Figure 2: PACE on K562 CRISPR perturbation data."""
import sys, pickle; sys.path.insert(0, ".")
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from pacefig import figstyle as fs
fs.setup(); C = fs.COL
R = pickle.load(open("results/k562_benchmark.pkl", "rb"))
T = pd.read_csv("results/k562_benchmark_metrics.csv"); D = pd.read_csv("results/k562_benchmark_paired_differences.csv")
SENS = pd.read_csv("results/k562_sensitivity.csv") if __import__("os").path.exists("results/k562_sensitivity.csv") else None
show = [("PACE", C["pace"], "-", "PACE"), ("ABC (ENCODE K562 predictions)", C["abc"], "-", "ABC"),
        ("ABC power-law contact (ENCODE)", C["abc2"], "--", "ABC, power-law contact"),
        ("GBM (TargetFinder-style)", C["gbm"], "-", "GBM (TargetFinder-style)"),
        ("Inverse distance", C["dist"], "-", "Inverse distance"), ("Distance prior P_s(d)", "#BDBDBD", "--", r"Distance prior $P_s(d)$")]
fig = plt.figure(figsize=(183 * fs.MM, 125 * fs.MM))
gs = fig.add_gridspec(2, 12, hspace=0.62, wspace=3.2, left=0.06, right=0.995, top=0.95, bottom=0.12)
def pr(ax, sname, letter, title):
    for name, col, ls, lab in show:
        rc, p = R["curves"][(sname, name)]
        t = T[(T.set == sname) & (T.method == name)].iloc[0]
        ax.plot(rc, p, color=col, ls=ls, lw=1.3 if name == "PACE" else 0.9, label=f"{lab} ({t.auprc:.2f})", zorder=5 if name == "PACE" else 2)
    n = T[(T.set == sname)].iloc[0]
    ax.axhline(n.n_pos / n.n, color="k", lw=0.5, ls=":")
    ax.set_xlabel("Recall"); ax.set_ylabel("Precision"); ax.set_xlim(0, 1); ax.set_ylim(0, 1.02)
    ax.set_title(title, fontsize=6.3, pad=3); fs.panel(ax, letter, -0.2, 1.05)
    ax.legend(loc="lower left", bbox_to_anchor=(0.0, 0.05), fontsize=5, handlelength=1.6, labelspacing=0.2, borderaxespad=0.2)
ax = fig.add_subplot(gs[0, 0:4]); pr(ax, "main", "a", "Combined benchmark (10,356 pairs, 471 positives)")
ax = fig.add_subplot(gs[0, 4:8]); pr(ax, "heldout", "b", "Held-out set (1,918 pairs, 118 positives)")
# c: paired differences on held-out
ax = fig.add_subplot(gs[0, 8:12]); fs.panel(ax, "c", -0.7, 1.05)
d = D[D.set == "heldout"].copy()
lab = {"ABC (ENCODE K562 predictions)": "ABC", "ABC power-law contact (ENCODE)": "ABC, power-law", "Inverse distance": "Inverse distance",
       "Distance prior P_s(d)": r"Distance prior", "GBM (TargetFinder-style)": "GBM"}
d = d.iloc[::-1]
y = np.arange(len(d))
ax.errorbar(d.delta, y, xerr=[d.delta - d.lo, d.hi - d.delta], fmt="o", ms=3, color=C["pace"], ecolor=C["pace"], elinewidth=1, capsize=1.5)
ax.axvline(0, color="k", lw=0.6)
ax.set_yticks(y); ax.set_yticklabels([lab[m] for m in d.method]); ax.set_xlabel("ΔAUPRC (PACE − method)")
for yy, (_, r) in zip(y, d.iterrows()): ax.text(0.43, yy, f"n={int(r.n_common)}", fontsize=4.8, va="center", color=C["mid"])
ax.set_xlim(-0.2, 0.55); ax.set_title("Held-out, pairs scored by both", fontsize=6.3, pad=3)
# d: coverage
ax = fig.add_subplot(gs[1, 0:4]); fs.panel(ax, "d", -0.2, 1.05)
meths = [s[0] for s in show]; labs = [s[3] for s in show]
for sn, off, mk in [("main", -0.16, "o"), ("heldout", 0.16, "s")]:
    t = T[T.set == sn].set_index("method").loc[meths]
    yy = np.arange(len(meths)) + off
    ax.scatter(t.coverage, yy, marker=mk, s=14, color=C["ink"], label=("all pairs, " + ("combined" if sn == "main" else "held-out")), zorder=3)
    ax.scatter(t.cov_pos, yy, marker=mk, s=14, facecolor="white", edgecolor=C["pace"], lw=0.9, label=("positives, " + ("combined" if sn == "main" else "held-out")), zorder=3)
ax.set_yticks(range(len(meths))); ax.set_yticklabels(labs); ax.invert_yaxis(); ax.set_xlim(0.87, 1.005); ax.set_xlabel("Fraction of pairs scored")
ax.grid(axis="y", lw=0.3, color="#DDDDDD"); ax.legend(fontsize=4.8, loc="upper center", bbox_to_anchor=(0.45, -0.28), ncol=2, handletextpad=0.2, columnspacing=0.8)
# e, f: component analysis
comp = [("PACE", "Full model"), ("PACE, accessibility only", "Accessibility only"), ("PACE, H3K27ac only", "H3K27ac only"),
        ("PACE, single canonical TSS", "Single canonical TSS"), ("ABC (ENCODE K562 predictions)", "ABC (ENCODE)")]
for j, (sn, letter, ttl) in enumerate([("main", "e", "Combined benchmark"), ("heldout", "f", "Held-out set")]):
    ax = fig.add_subplot(gs[1, 4 + 4 * j:8 + 4 * j]); fs.panel(ax, letter, -0.2, 1.05)
    t = T[T.set == sn].set_index("method").loc[[c[0] for c in comp]]
    cols = [C["pace"], "#E69F00", "#CC79A7", "#F0A070", C["abc"]]
    x = np.arange(len(comp))
    ax.bar(x, t.auprc, 0.62, color=cols)
    ax.errorbar(x, t.auprc, yerr=[t.auprc - t.auprc_lo, t.auprc_hi - t.auprc], fmt="none", ecolor="k", elinewidth=0.7, capsize=1.5)
    for xi, v in zip(x, t.auprc): ax.text(xi, 0.02, f"{v:.2f}", ha="center", color="white", fontsize=5.5, fontweight="bold")
    ax.axhline(t.iloc[0].n_pos / t.iloc[0].n, color="k", lw=0.5, ls=":")
    ax.set_xticks(x); ax.set_xticklabels([c[1] for c in comp], rotation=35, ha="right"); ax.set_ylabel("AUPRC"); ax.set_ylim(0, 0.8)
    ax.set_title(ttl, fontsize=6.3, pad=3)
fs.save(fig, "Fig2_K562_CRISPR_benchmark")
