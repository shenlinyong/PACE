"""Supplementary Figs S2, S4, S5, S6, S7 (human K562/GM12878 parts)."""
import sys; sys.path.insert(0, ".")
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from pacefig import figstyle as fs, k562, engine
fs.setup(); C = fs.COL
# ------------------------------------------------------------------ S2 distance strata
S = pd.read_csv("results/k562_benchmark_distance_strata.csv")
meths = [("PACE", C["pace"], "PACE"), ("ABC (ENCODE K562 predictions)", C["abc"], "ABC"), ("GBM (TargetFinder-style)", C["gbm"], "GBM"), ("Inverse distance", C["dist"], "Inverse distance")]
strata = ["<50 kb", "50-200 kb", ">200 kb", "distance-matched negatives"]
fig, axs = plt.subplots(1, 2, figsize=(183 * fs.MM, 62 * fs.MM), sharey=True)
for ax, sn, ttl, L in zip(axs, ["main", "heldout"], ["Combined benchmark", "Held-out set"], "ab"):
    fs.panel(ax, L, -0.14, 1.06)
    for j, (m, col, lab) in enumerate(meths):
        d = S[(S.set == sn) & (S.method == m)].set_index("stratum").loc[strata]
        x = np.arange(4) + (j - 1.5) * 0.19
        ax.bar(x, d.auprc, 0.17, color=col, label=lab)
        ax.errorbar(x, d.auprc, yerr=[d.auprc - d.lo, d.hi - d.auprc], fmt="none", ecolor="k", elinewidth=0.5, capsize=0.8)
    prev = S[(S.set == sn) & (S.method == "PACE")].set_index("stratum").loc[strata].prevalence
    ax.set_xticks(range(4)); ax.set_xticklabels([f"{s}\n(prev. {p:.3f})" for s, p in zip(["<50 kb", "50–200 kb", ">200 kb", "distance-matched\nnegatives"], prev)], fontsize=5)
    ax.set_title(ttl, fontsize=6.3, pad=3)
axs[0].set_ylabel("AUPRC"); axs[0].legend(fontsize=5, ncol=2, loc="upper right")
fs.save(fig, "FigS2_distance_strata")
# ------------------------------------------------------------------ S4 decay fits
bsK = k562.BandSet.load(k562.W + "band/r10000", 10000); gm = k562.BandSet.load("/home/user/work/gm12878/band/r10000", 10000)
fig, ax = plt.subplots(figsize=(89 * fs.MM, 66 * fs.MM)); fs.panel(ax, "a", -0.2, 1.05)
for name, bs, col in [("K562", bsK, fs.CASE["A"]), ("GM12878", gm, fs.CASE["B_same"])]:
    p, rows = engine.fit_prior(bs.bands, bs.weights, 10000)
    r = pd.DataFrame(rows); r = r[r.fit_status == "used"]
    ax.plot(r.geometric_distance, r.mean_contact, "o", ms=2.6, color=col)
    xx = np.logspace(4, np.log10(5e6), 50); ax.plot(xx, p(xx) * (1 if True else 1), color=col, label=f"{name}: γ = {p.gamma:.3f}")
hum = engine.human_reference_prior(); xx = np.logspace(4, np.log10(5e6), 50)
p0, _ = engine.fit_prior(bsK.bands, bsK.weights, 10000)
ax.plot(xx, p0(1e4) * (xx / 1e4) ** (-hum.gamma), color=fs.CASE["C"], ls="--", label="human reference: γ = 1.024\n(scaled to K562 at 10 kb)")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("Distance (bp)"); ax.set_ylabel("Mean balanced contact"); ax.legend(fontsize=5)
fs.save(fig, "FigS4_distance_decay_human")
# ------------------------------------------------------------------ S5 sensitivity
R = pd.read_csv("results/k562_sensitivity.csv"); R = R[R.parameter != "default"].copy(); R["value"] = R.value.astype(float)
params = ["candidate window (Mb)", "expression threshold (TPM)", "pseudocount distance d0 (kb)", "element width (bp)"]
fig, axs = plt.subplots(2, 4, figsize=(183 * fs.MM, 80 * fs.MM), gridspec_kw=dict(hspace=0.75, wspace=0.5))
for j, p in enumerate(params):
    for i, sn in enumerate(["main", "heldout"]):
        ax = axs[i, j]; d = R[(R.parameter == p) & (R.set == sn)].sort_values("value")
        x = np.arange(len(d)); ax.errorbar(x, d.auprc, yerr=[d.auprc - d.lo, d.hi - d.auprc], fmt="o-", color=C["pace"], ms=3, lw=0.8, capsize=1.3, elinewidth=0.7)
        ax.set_xticks(x); ax.set_xticklabels([f"{v:g}" for v in d.value]); ax.set_ylim(0.35, 0.75)
        ax.set_title(("Combined: " if sn == "main" else "Held-out: ") + p.split(" (")[0], fontsize=5.5, pad=2)
        ax.set_xlabel(p.split("(")[1].rstrip(")") if "(" in p else "", fontsize=5.5)
        if j == 0: ax.set_ylabel("AUPRC")
fs.panel(axs[0, 0], "a", -0.4, 1.15); fs.panel(axs[1, 0], "b", -0.4, 1.15)
fs.save(fig, "FigS5_sensitivity")
# ------------------------------------------------------------------ S6 strata of contact cases + replicates
D = pd.read_csv("results/k562_degradation.csv"); M = D[D.set == "main"].set_index("condition")
fig, axs = plt.subplots(1, 2, figsize=(183 * fs.MM, 58 * fs.MM), gridspec_kw=dict(width_ratios=[2.2, 1], wspace=0.35))
ax = axs[0]; fs.panel(ax, "a", -0.1, 1.06)
labs = ["<50 kb", "50-200 kb", ">200 kb"]
for j, k in enumerate(["A", "B_same", "B_other", "C"]):
    d = [M.loc[f"strat_{l}_{k}"] for l in labs]; x = np.arange(3) + (j - 1.5) * 0.2
    ax.bar(x, [r.auprc for r in d], 0.18, color=fs.CASE[k], label={"A": "A measured", "B_same": "B same-cell", "B_other": "B other-cell", "C": "C human"}[k])
    ax.errorbar(x, [r.auprc for r in d], yerr=[[r.auprc - r.lo for r in d], [r.hi - r.auprc for r in d]], fmt="none", ecolor="k", elinewidth=0.5, capsize=0.8)
ax.set_xticks(range(3)); ax.set_xticklabels(["<50 kb", "50–200 kb", ">200 kb"]); ax.set_ylabel("AUPRC"); ax.legend(fontsize=5, ncol=2); ax.set_title("Contact cases by element–TSS distance", fontsize=6.3, pad=3)
ax = axs[1]; fs.panel(ax, "b", -0.3, 1.06)
for i, (k, l) in enumerate([("A", "2 replicates"), ("rep1", "1 replicate")]):
    r = M.loc[k]; ax.bar(i, r.auprc, 0.6, color=fs.CASE["A"] if i == 0 else "#E8A57A"); ax.errorbar(i, r.auprc, yerr=[[r.auprc - r.lo], [r.hi - r.auprc]], fmt="none", ecolor="k", elinewidth=0.7, capsize=1.5)
ax.set_xticks([0, 1]); ax.set_xticklabels(["2 replicates", "1 replicate"]); ax.set_ylim(0, 0.8); ax.set_title("Biological replicates", fontsize=6.3, pad=3)
fs.save(fig, "FigS6_distance_strata_and_replicates")
# ------------------------------------------------------------------ S7 resolution
fig, axs = plt.subplots(1, 2, figsize=(120 * fs.MM, 55 * fs.MM), gridspec_kw=dict(wspace=0.4))
for ax, sn, L in zip(axs, ["main", "heldout"], "ab"):
    m = D[D.set == sn].set_index("condition"); fs.panel(ax, L, -0.25, 1.06)
    r = [m.loc[f"res{x}"] for x in [5000, 10000, 25000, 40000]]; x = np.arange(4)
    ax.errorbar(x, [a.auprc for a in r], yerr=[[a.auprc - a.lo for a in r], [a.hi - a.auprc for a in r]], fmt="o-", color=fs.CASE["A"], ms=3.5, capsize=1.5, elinewidth=0.8)
    ax.set_xticks(x); ax.set_xticklabels(["5", "10", "25", "40"]); ax.set_xlabel("Hi-C resolution (kb)"); ax.set_ylabel("AUPRC"); ax.set_ylim(0.4, 0.75)
    ax.set_title("Combined benchmark" if sn == "main" else "Held-out set", fontsize=6.3, pad=3)
fs.save(fig, "FigS7_resolution")
