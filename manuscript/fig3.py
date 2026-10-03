"""Figure 3a-c: controlled degradation of complete K562 data (livestock panels d,e require livestock eQTL data)."""
import sys; sys.path.insert(0, ".")
import numpy as np, pandas as pd, pickle, matplotlib.pyplot as plt
from pacefig import figstyle as fs
fs.setup(); C = fs.COL
D = pd.read_csv("results/k562_degradation.csv"); M = D[D.set == "main"].set_index("condition")
fig = plt.figure(figsize=(183 * fs.MM, 62 * fs.MM))
gs = fig.add_gridspec(1, 4, width_ratios=[1.1, 1.2, 0.95, 0.95], wspace=0.62, left=0.065, right=0.995, top=0.86, bottom=0.25)
# a
ax = fig.add_subplot(gs[0]); fs.panel(ax, "a", -0.32, 1.12)
conds = [("A", "A\nmeasured\nHi-C"), ("B_same", "B\nsame-\ncell"), ("B_other", "B\nother-\ncell"), ("C", "C\nhuman\nref.")]
for i, (k, lab) in enumerate(conds):
    r = M.loc[k]; ax.bar(i, r.auprc, 0.62, color=fs.CASE[k])
    ax.errorbar(i, r.auprc, yerr=[[r.auprc - r.lo], [r.hi - r.auprc]], fmt="none", ecolor="k", elinewidth=0.7, capsize=1.5)
    ax.text(i, 0.02, f"{r.auprc:.3f}", ha="center", color="white", fontsize=5.3, fontweight="bold")
ax.set_xticks(range(4)); ax.set_xticklabels([c[1] for c in conds], fontsize=5); ax.set_ylim(0, 0.8); ax.set_ylabel("AUPRC")
ax.set_title("Contact case", fontsize=6.3, pad=3)
d1 = M.loc["diff_A_minus_Bsame"]; d2 = M.loc["diff_Bsame_minus_Bother"]
ax.text(0.02, 0.97, f"A − same-cell: {d1.auprc:+.3f}\n[{d1.lo:+.3f}, {d1.hi:+.3f}]\nsame − other: {d2.auprc:+.4f}\n[{d2.lo:+.4f}, {d2.hi:+.4f}]", transform=ax.transAxes, fontsize=4.6, va="top", color=C["mid"])
ax.set_ylim(0, 0.95)
# b
ax = fig.add_subplot(gs[1]); fs.panel(ax, "b", -0.28, 1.12)
fr = [100, 50, 25, 10, 5]; xs = np.arange(len(fr))
mean = [M.loc["depth100_ice"].auprc] + [M.loc[f"depth{d}_mean"].auprc for d in fr[1:]]
lo = [M.loc["depth100_ice"].lo] + [M.loc[f"depth{d}_mean"].lo for d in fr[1:]]
hi = [M.loc["depth100_ice"].hi] + [M.loc[f"depth{d}_mean"].hi for d in fr[1:]]
ax.axhline(M.loc["B_other"].auprc, color=fs.CASE["B_other"], ls="--", lw=0.9); ax.text(-0.3, 0.46, "dashed: other-cell prior (95% CI shaded)", fontsize=4.8, ha="left", color=fs.CASE["B_other"], va="bottom")
ax.axhspan(M.loc["B_other"].lo, M.loc["B_other"].hi, color=fs.CASE["B_other"], alpha=0.12, lw=0)
for d_, x in zip(fr[1:], xs[1:]):
    for s in [1, 2, 3]: ax.plot(x, M.loc[f"depth{d_}_s{s}"].auprc, ".", color=C["mid"], ms=2.5, zorder=3)
ax.errorbar(xs, mean, yerr=[np.array(mean) - lo, np.array(hi) - mean], fmt="o-", color=fs.CASE["A"], ms=3.5, lw=1, capsize=1.5, elinewidth=0.8, zorder=4)
ax.set_xticks(xs); ax.set_xticklabels([f"{f}%" for f in fr]); ax.set_xlabel("Hi-C read pairs retained"); ax.set_ylabel("AUPRC"); ax.set_ylim(0.45, 0.72)
ax.set_title("Hi-C depth (prior refitted)", fontsize=6.3, pad=3)
# c: left AUPRC, right coverage
ax = fig.add_subplot(gs[2]); fs.panel(ax, "c", -0.36, 1.12)
fr2 = [5, 10, 20]; x = np.arange(3)
cond = [M.loc[f"mask{f}"] for f in fr2]
ax.errorbar(x - 0.08, [r.auprc for r in cond], yerr=[[r.auprc - r.lo for r in cond], [r.hi - r.auprc for r in cond]], fmt="o", color=fs.CASE["A"], ms=3.5, capsize=1.5, elinewidth=0.8, label="conditional policy")
ax.errorbar(x + 0.08, [r.ref_auprc for r in cond], yerr=[[r.ref_auprc - r.ref_lo for r in cond], [r.ref_hi - r.ref_auprc for r in cond]], fmt="s", mfc="white", color=C["ink"], ms=3.3, capsize=1.5, elinewidth=0.8, label="complete data, same pairs")
ax.set_xticks(x); ax.set_xticklabels([f"{f}%" for f in fr2]); ax.set_xlabel("Genome with activity removed"); ax.set_ylabel("AUPRC (affected genes)"); ax.set_ylim(0.4, 0.8)
ax.legend(loc="lower left", fontsize=4.6, handletextpad=0.2); ax.set_title("Missing activity", fontsize=6.3, pad=3)
ax = fig.add_subplot(gs[3])
w = 0.34
ax.bar(x - w / 2, [r.cov_default for r in cond], w, yerr=[r.cov_default_sd for r in cond], color="#9A9A9A", error_kw=dict(lw=0.7, capsize=1.2), label="default (withhold)")
ax.bar(x + w / 2, [r.cov_cond for r in cond], w, yerr=[r.cov_cond_sd for r in cond], color=fs.CASE["A"], error_kw=dict(lw=0.7, capsize=1.2), label="conditional")
ax.set_xticks(x); ax.set_xticklabels([f"{f}%" for f in fr2]); ax.set_xlabel("Genome with activity removed"); ax.set_ylabel("Coverage (fraction of pairs scored)"); ax.set_ylim(0, 1.18)
ax.legend(loc="upper right", fontsize=4.6, handletextpad=0.3); ax.set_title("Coverage", fontsize=6.3, pad=3)
fig.text(0.005, 0.01, "Panels d, e (livestock eQTL proxy) are not shown: they require FAANG signal tracks, liver Hi-C and FarmGTEx fine-mapping that were not processed.", fontsize=4.8, color=C["mid"])
fs.save(fig, "Fig3_abc_K562_degradation")
