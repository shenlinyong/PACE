"""Figure 1: the PACE model, one illustrative locus; every number follows from Equation 1."""
import sys; sys.path.insert(0, ".")
import numpy as np, matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, FancyArrowPatch, Circle
from pacefig import figstyle as fs
fs.setup()
C = fs.COL
# ---------------------------------------------------------------- worked example (illustrative values)
atac = dict(E1=8, E2=6, E3=5); k27 = dict(E1=2, E2=6, E3=0)
A = {e: float(np.sqrt(atac[e] * k27[e])) for e in atac}
# contact of each element with the TSSs of G1 (two TSSs) and with G2 (one TSS); None = unresolved
c = {("E1", "G1"): [2, 4], ("E2", "G1"): [1, 2], ("E3", "G1"): [3, 3],
     ("E1", "G2"): [5], ("E2", "G2"): [1], ("E3", "G2"): [None]}
cbar = {k: (None if any(x is None for x in v) else float(np.mean(v))) for k, v in c.items()}
supp = {k: (None if v is None else A[k[0]] * v) for k, v in cbar.items()}
score = {}
for g in ["G1", "G2"]:
    ks = [k for k in supp if k[1] == g]
    if any(supp[k] is None for k in ks):
        for k in ks: score[k] = None
    else:
        tot = sum(supp[k] for k in ks)
        for k in ks: score[k] = supp[k] / tot
assert abs(sum(score[k] for k in score if k[1] == "G1") - 1) < 1e-12
EC = {"E1": "#D55E00", "E2": "#0072B2", "E3": "#009E73"}

fig = plt.figure(figsize=(183 * fs.MM, 170 * fs.MM))
gs = fig.add_gridspec(3, 3, height_ratios=[1.05, 1.0, 1.05], hspace=0.5, wspace=0.42, left=0.06, right=0.995, top=0.965, bottom=0.04)

# ---- a: locus
ax = fig.add_subplot(gs[0, :2]); fs.panel(ax, "a", -0.07, 1.04)
ax.set_xlim(0, 100); ax.set_ylim(-3, 19); ax.axis("off")
ax.plot([2, 98], [0, 0], color=C["ink"], lw=0.7)
pos = {"E1": 18, "E2": 40, "E3": 72}
tss = {"G1a": 52, "G1b": 58, "G2": 90}
for e, x in pos.items():
    ax.add_patch(Rectangle((x - 3, 0), 6, 1.4, fc=EC[e], ec="none"))
    ax.text(x, -1.6, e, ha="center", va="center", color=EC[e], fontweight="bold")
    # signal tracks
    ax.bar([x - 1.2], [atac[e] * 0.45], width=2.2, bottom=3, color="#7F7F7F"); ax.bar([x + 1.2], [k27[e] * 0.45], width=2.2, bottom=3, color="#E69F00")
ax.text(1, 3.2, "ATAC", color="#7F7F7F", fontsize=5.5, va="bottom"); ax.text(1, 5.4, "H3K27ac", color="#E69F00", fontsize=5.5, va="bottom")
for name, x in tss.items():
    ax.plot([x, x], [0, 2.6], color=C["ink"], lw=0.7); ax.annotate("", xy=(x + 3, 2.6), xytext=(x, 2.6), arrowprops=dict(arrowstyle="->", lw=0.7, color=C["ink"]))
ax.add_patch(Rectangle((50, 0), 10, 1.0, fc="#B8B8B8", ec="none")); ax.text(55, -1.6, r"$G_1$", ha="center", fontweight="bold")
ax.add_patch(Rectangle((88, 0), 5, 1.0, fc="#B8B8B8", ec="none")); ax.text(90.5, -1.6, r"$G_2$", ha="center", fontweight="bold")
ax.text(52, 3.1, r"$t_1$", ha="center", fontsize=5.5); ax.text(58.5, 3.1, r"$t_2$", ha="center", fontsize=5.5)
def arc(x0, x1, h, ls, col, lw=0.9):
    xs = np.linspace(x0, x1, 60); ys = 6.2 + h * np.sin(np.pi * (xs - x0) / (x1 - x0))
    ax.plot(xs, ys + 3.2, color=col, lw=lw, ls=ls)
arc(pos["E1"], 52, 2.6, "-", EC["E1"]); arc(pos["E1"], 58, 3.4, "-", EC["E1"]); arc(pos["E2"], 52, 1.4, "-", EC["E2"]); arc(pos["E2"], 58, 2.0, "-", EC["E2"])
arc(pos["E3"], 52, 1.6, "-", EC["E3"]); arc(pos["E3"], 58, 1.0, "-", EC["E3"])
arc(pos["E1"], 90, 4.2, "-", EC["E1"], 0.6); arc(pos["E2"], 90, 3.0, "-", EC["E2"], 0.6); arc(pos["E3"], 90, 0.8, (0, (2, 1.5)), EC["E3"])
ax.text(81, 13.2, "contact unresolved (NA)", fontsize=5.3, color=EC["E3"], ha="center")
ax.text(40, 17.2, r"$G_1$: all candidates resolved, scored", fontsize=5.5, color=C["mid"], ha="center")
ax.text(99, 16.0, r"$G_2$: one NA, all scores withheld", fontsize=5.5, color=C["mid"], ha="right")
ax.set_title("Candidate element–gene graph (fixed before scoring)", fontsize=6.5, loc="left", pad=2)

# ---- b: activity
ax = fig.add_subplot(gs[0, 2]); fs.panel(ax, "b", -0.22, 1.04)
x = np.arange(3); w = 0.26
ax.bar(x - w, [atac[e] for e in EC], w, color="#7F7F7F", label="ATAC"); ax.bar(x, [k27[e] for e in EC], w, color="#E69F00", label="H3K27ac")
ax.bar(x + w, [A[e] for e in EC], w, color=[EC[e] for e in EC], label=r"$A_\star$")
for i, e in enumerate(EC):
    ax.text(i + w, A[e] + 0.25, f"{A[e]:.1f}", ha="center", fontsize=5.5, color=EC[e])
ax.text(2, 0.35, "0 (measured)", ha="center", fontsize=5, rotation=90, va="bottom", color="#B07A00")
ax.set_xticks(x); ax.set_xticklabels(list(EC)); ax.set_ylabel("signal (a.u.)"); ax.set_ylim(0, 9.5)
ax.legend(ncol=1, loc="upper right", handlelength=1, borderaxespad=0)
ax.set_title(r"$A_\star=\sqrt{x_{ATAC}\,x_{H3K27ac}}$", fontsize=6.5, pad=2)

# ---- c: three contact cases
axc = fig.add_subplot(gs[1, :2]); fs.panel(axc, "c", -0.07, 1.04); axc.axis("off")
sub = gs[1, :2].subgridspec(1, 3, wspace=0.4)
rng = np.random.default_rng(3)
d = np.logspace(4, 6.7, 100)
for i, (lab, ttl, sub_t) in enumerate([("A", "sample Hi-C", r"$\tilde{C}=C_{obs}+p_s(d)$"), ("B", "same-species Hi-C", r"$\tilde{C}=P_s(d)$"), ("C", "no species Hi-C", r"$\tilde{C}=P_h(d),\ \gamma_h\approx1.024$")]):
    a2 = fig.add_subplot(sub[i]); a2.set_xscale("log"); a2.set_yscale("log")
    if lab == "A":
        base = 4e3 * d ** -1.13; obs = base * np.exp(rng.normal(0, 0.45, len(d)))
        obs[rng.choice(100, 14, replace=False)] *= 4
        a2.plot(d, obs, ".", ms=2.2, color=C["pace"]); a2.plot(d, base, color=C["ink"], lw=0.8, ls="--")
        a2.text(0.04, 0.05, r"$p_s(d)$ pseudocount", transform=a2.transAxes, fontsize=5, color=C["ink"])
    elif lab == "B":
        a2.plot(d, 4e3 * d ** -1.13, color=CASE_B if False else "#0072B2", lw=1.4)
        a2.text(0.5, 0.78, r"$\gamma$ fitted" + "\n" + "on species Hi-C", transform=a2.transAxes, fontsize=5, color="#0072B2")
    else:
        a2.plot(d, 4e3 * d ** -1.024, color="#8C8C8C", lw=1.4)
        a2.text(0.5, 0.78, "fixed human\nreference", transform=a2.transAxes, fontsize=5, color="#6B6B6B")
    a2.set_xlabel("distance (bp)"); a2.set_ylabel("contact" if i == 0 else ""); a2.set_xticks([1e4, 1e5, 1e6]); a2.set_yticks([])
    a2.set_title(f"{lab}: {ttl}\n{sub_t}", fontsize=6, pad=3)

# ---- d: TSS averaging
ax = fig.add_subplot(gs[1, 2]); fs.panel(ax, "d", -0.22, 1.04); ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 10)
ax.set_title(r"TSS-averaged contact, $E_1$–$G_1$", fontsize=6.5, pad=2)
ax.text(5, 8.6, r"$\overline{C}=\sum_t\pi(t|G)\,\tilde{C}(E,t)$", ha="center", fontsize=6.5)
for j, (lab, val, xx) in enumerate([(r"$t_1$", 2, 2.6), (r"$t_2$", 4, 7.4)]):
    ax.add_patch(Rectangle((xx - 1.1, 3.0), 2.2, val * 0.9, fc=EC["E1"], ec="none", alpha=0.85))
    ax.text(xx, 2.2, lab, ha="center"); ax.text(xx, 3.2 + val * 0.9, f"{val}", ha="center", fontsize=6)
ax.plot([0.6, 9.4], [3.0 + 3 * 0.9] * 2, color=C["ink"], lw=0.7, ls=":")
ax.text(5, 6.9, r"$\frac{1}{2}\times2+\frac{1}{2}\times4=3$", ha="center", fontsize=6.3)
ax.text(5, 1.0, "equal weights by default", ha="center", fontsize=5.5, color=C["mid"])

# ---- e: score
ax = fig.add_subplot(gs[2, 0]); fs.panel(ax, "e", -0.28, 1.04)
ks = [("E1", "G1"), ("E2", "G1"), ("E3", "G1")]
su = [supp[k] for k in ks]; sc = [score[k] for k in ks]
b = ax.bar(range(3), su, color=[EC[k[0]] for k in ks], width=0.6)
for i, k in enumerate(ks):
    ax.text(i, su[i] + 0.4, f"{su[i]:.0f}" + ("" if su[i] else ""), ha="center", fontsize=5.8)
    ax.text(i, -3.3, f"{sc[i]:.2f}", ha="center", fontsize=6.5, fontweight="bold", color=EC[k[0]], clip_on=False)
ax.set_xticks(range(3)); ax.set_xticklabels([k[0] for k in ks]); ax.set_ylabel(r"$A_\star\times\overline{C}$"); ax.set_ylim(0, 16)
ax.text(1, -5.6, "PACE score (sums to 1)", ha="center", fontsize=5.5, color=C["mid"], clip_on=False)
ax.set_title(r"Normalisation within $G_1$", fontsize=6.5, pad=2)

# ---- f: output table
ax = fig.add_subplot(gs[2, 1:]); fs.panel(ax, "f", -0.045, 1.04); ax.axis("off"); ax.set_xlim(0, 100); ax.set_ylim(0, 8.4)
ax.set_title("Worked output for all six candidate links", fontsize=6.5, loc="left", pad=2)
cols = [("Element", 1), ("Gene", 12), (r"$A_\star$", 22), (r"$\overline{C}$", 34), (r"$A_\star\overline{C}$", 46), ("PACE score", 62), ("Status", 80)]
for t, x0 in cols: ax.text(x0, 7.5, t, fontweight="bold", fontsize=6)
ax.plot([0, 100], [7.0, 7.0], color=C["ink"], lw=0.6)
for r, k in enumerate([("E1", "G1"), ("E2", "G1"), ("E3", "G1"), ("E1", "G2"), ("E2", "G2"), ("E3", "G2")]):
    yy = 6.1 - r * 1.0
    ax.add_patch(Rectangle((0, yy - 0.35), 100, 0.9, fc="#F4F4F4" if r % 2 == 0 else "white", ec="none", zorder=0))
    cb = cbar[k]; sp = supp[k]; s_ = score[k]
    vals = [k[0], k[1], f"{A[k[0]]:.2f}".rstrip("0").rstrip("."), "NA" if cb is None else f"{cb:.1f}".rstrip("0").rstrip("."),
            "NA" if sp is None else f"{sp:.1f}".rstrip("0").rstrip("."), "withheld" if s_ is None else f"{s_:.3f}",
            "scored" if k[1] == "G1" else ("withheld (NA)" )]
    for (t, x0), v in zip(cols, vals):
        ax.text(x0, yy, v, fontsize=6, color=EC[k[0]] if t == "Element" else C["ink"], fontweight="bold" if t == "Element" else "normal")
ax.text(0, -0.2, "A measured zero (E3) is a valid observation; a single unresolved contact (NA) withholds every score of that gene. Values are illustrative.", fontsize=5.3, color=C["mid"])
fs.save(fig, "Fig1_PACE_model")
print({k: (round(v, 3) if v is not None else None) for k, v in score.items()})
