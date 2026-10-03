"""Controlled degradation of the complete K562 data (Fig. 3a-c; Supplementary Figs S4, S6, S7; Table S14)."""
import sys, time, pickle
sys.path.insert(0, ".")
import numpy as np, pandas as pd
from pacefig import k562, bench, engine

W = k562.W
t0 = time.time()
K = k562.K562()
P = K.pairs
y = P.y.to_numpy()
dist = P.distanceToTSS.to_numpy().astype(float)
sets = {"main": (P.set == "main").to_numpy(), "heldout": (P.set == "heldout").to_numpy()}
bs = k562.BandSet.load(W + "band/r10000", 10000)
p_k562 = bs.fit()
gm = k562.BandSet.load("/home/user/work/gm12878/band/r10000", 10000)
p_gm = gm.fit()
p_h = engine.human_reference_prior()
print("gamma K562 %.4f, GM12878 %.4f, human reference %.4f" % (p_k562.gamma, p_gm.gamma, p_h.gamma), flush=True)

S = {}      # condition -> pair scores
meta = {}
def score(name, **kw):
    S[name] = K.pair_scores(K.run(**kw))

# --- contact cases
score("A", mode="A", prior=p_k562, bands=bs)
score("B_same", mode="B", prior=p_k562)
score("B_other", mode="B", prior=p_gm)
score("C", mode="C", prior=p_h)
print("contact cases %.0fs" % (time.time() - t0), flush=True)

# --- Hi-C depth: binomial thinning, balancing and prior refit on each thinned map
bs100 = bs.rebalanced()
score("depth100_ice", mode="A", prior=bs100.fit(), bands=bs100)
for frac in [0.5, 0.25, 0.1, 0.05]:
    for seed in [1, 2, 3]:
        tb = bs.thinned(frac, seed)
        pr = tb.fit()
        meta[f"depth{int(frac*100)}_s{seed}"] = pr.gamma
        score(f"depth{int(frac*100)}_s{seed}", mode="A", prior=pr, bands=tb)
        del tb
    print("depth", frac, "%.0fs" % (time.time() - t0), flush=True)

# --- resolution (Supplementary Fig. S7)
for res in [5000, 10000, 25000, 40000]:
    if res == 10000:
        S["res10000"] = S["A"]; meta["res10000"] = p_k562.gamma; continue
    b = k562.BandSet.load(W + f"band/r{res}", res)
    if res == 40000:
        b = k562.BandSet(b.bands, {c: engine.ice_balance(x, maxd_bins=5_000_000 // res + 2) for c, x in b.bands.items()}, res)
    pr = b.fit(); meta[f"res{res}"] = pr.gamma
    score(f"res{res}", mode="A", prior=pr, bands=b)
    print("res", res, "%.0fs" % (time.time() - t0), flush=True)

# --- one biological replicate per assay
score("rep1", mode="A", prior=p_k562, bands=bs, A=K.activity("both_rep1"))
# prior-only with a single replicate too
score("B_same_rep1", mode="B", prior=p_k562, A=K.activity("both_rep1"))
score("B_same_rep2both", mode="B", prior=p_k562)
print("rep %.0fs" % (time.time() - t0), flush=True)

# --- missing activity: random 1 Mb blocks covering 5/10/20% of the genome
chrom_sizes = pd.read_csv("/home/user/data/ref/hg38.chrom.sizes", sep="\t", header=None, index_col=0)[1].to_dict()
blocks_all = [(c, i) for c in k562.CHROMS for i in range(chrom_sizes[c] // 1_000_000)]
A0 = K.activity("both")
units = K.cat.units
gene_pair = P.gene_id.to_numpy()
mask_res = {}
cands = K.cands()
for frac in [0.05, 0.10, 0.20]:
    for seed in [1, 2, 3]:
        rng = np.random.default_rng(100 * seed + int(frac * 100))
        n = int(round(frac * len(blocks_all)))
        pick = rng.choice(len(blocks_all), n, replace=False)
        masked = set(blocks_all[i] for i in pick)
        key = list(zip(units.chrom.to_numpy(), units.start.to_numpy() // 1_000_000))
        m = np.fromiter(((c, b) in masked for c, b in key), bool, len(units))
        A = A0.copy(); A[m] = np.nan
        tab = K.run("A", p_k562, A=A, bands=bs)
        s_w = K.pair_scores(tab, "score"); s_c = K.pair_scores(tab, "cond")
        # genes with at least one masked candidate
        aff = set(tab.gene_id[np.isnan(A[tab.unit_idx.to_numpy()])])
        mask_res[(frac, seed)] = dict(default=s_w, cond=s_c, affected=np.array([g in aff for g in gene_pair]),
                                      masked_fraction=float(m.mean()))
    print("mask", frac, "%.0fs" % (time.time() - t0), flush=True)

res = dict(S=S, meta=meta, mask=mask_res, gamma=dict(k562=p_k562.gamma, gm12878=p_gm.gamma, human=p_h.gamma))
pickle.dump(res, open("results/k562_degradation_scores.pkl", "wb"))

# -------------------------------------------------------------------------- evaluation
rows = []
boots = {k: bench.Boot(P.gene_id.to_numpy()[m]) for k, m in sets.items()}
contact = ["A", "B_same", "B_other", "C"]
depth = ["depth100_ice"] + [f"depth{d}_s{s}" for d in [50, 25, 10, 5] for s in [1, 2, 3]]
others = ["res5000", "res10000", "res25000", "res40000", "rep1"]
def common(sn, names):
    m = sets[sn].copy()
    for n in names:
        m &= np.isfinite(S[n])
    return m

for sn, base in sets.items():
    B = boots[sn]
    cm = common(sn, contact + depth + others)
    wsub = cm[base].astype(float)
    print(sn, "pairs scored under all conditions:", int(cm.sum()), "positives", int(y[cm].sum()))
    ys = y[base]
    def boot_ap(name, mask=None):
        s = S[name][base]
        w = wsub if mask is None else mask[base].astype(float)
        return B.stat(ys, s, w)
    bt = {n: boot_ap(n) for n in contact + depth + others + ["B_same_rep1"]}
    for n, v in bt.items():
        point = bench.ap(y[cm], S[n][cm])
        rows.append(dict(set=sn, condition=n, auprc=point, lo=np.nanpercentile(v, 2.5), hi=np.nanpercentile(v, 97.5),
                         n_common=int(cm.sum()), pos_common=int(y[cm].sum())))
    # depth: mean over seeds (paired bootstrap)
    for d in [50, 25, 10, 5]:
        names = [f"depth{d}_s{s}" for s in [1, 2, 3]]
        v = np.mean([bt[n] for n in names], axis=0)
        point = np.mean([bench.ap(y[cm], S[n][cm]) for n in names])
        sd = np.std([bench.ap(y[cm], S[n][cm]) for n in names], ddof=1)
        rows.append(dict(set=sn, condition=f"depth{d}_mean", auprc=point, lo=np.nanpercentile(v, 2.5), hi=np.nanpercentile(v, 97.5),
                         n_common=int(cm.sum()), pos_common=int(y[cm].sum()), seed_sd=sd,
                         gamma_mean=np.mean([meta[n] for n in names])))
    # paired differences
    def dif(a, b_, label):
        v = bt[a] - bt[b_] if isinstance(a, str) else a - bt[b_]
        p = bench.ap(y[cm], S[a][cm]) - bench.ap(y[cm], S[b_][cm])
        rows.append(dict(set=sn, condition=label, auprc=p, lo=np.nanpercentile(v, 2.5), hi=np.nanpercentile(v, 97.5),
                         n_common=int(cm.sum()), pos_common=int(y[cm].sum())))
    dif("A", "B_same", "diff_A_minus_Bsame")
    dif("B_same", "B_other", "diff_Bsame_minus_Bother")
    dif("B_other", "C", "diff_Bother_minus_C")
    dif("A", "B_other", "diff_A_minus_Bother")
    dif("A", "rep1", "diff_2rep_minus_1rep")
    dif("B_same", "B_same_rep1", "diff_Bsame_2rep_minus_1rep") if False else None
    # does thinned Hi-C stay above the other-cell prior?
    for d in [50, 25, 10, 5]:
        names = [f"depth{d}_s{s}" for s in [1, 2, 3]]
        v = np.mean([bt[n] for n in names], axis=0) - bt["B_other"]
        p = np.mean([bench.ap(y[cm], S[n][cm]) for n in names]) - bench.ap(y[cm], S["B_other"][cm])
        rows.append(dict(set=sn, condition=f"diff_depth{d}_minus_Bother", auprc=p, lo=np.nanpercentile(v, 2.5),
                         hi=np.nanpercentile(v, 97.5), n_common=int(cm.sum()), pos_common=int(y[cm].sum())))
    # distance-stratified comparison of contact cases (Supplementary Fig. S6)
    for lab, lo_, hi_ in [("<50 kb", 0, 50e3), ("50-200 kb", 50e3, 200e3), (">200 kb", 200e3, 1e12)]:
        sm = cm & (dist >= lo_) & (dist < hi_)
        if y[sm].sum() < 3:
            continue
        for n in contact:
            rows.append(dict(set=sn, condition=f"strat_{lab}_{n}", auprc=bench.ap(y[sm], S[n][sm]), n_common=int(sm.sum()),
                             pos_common=int(y[sm].sum()),
                             lo=np.nanpercentile(B.stat(ys, S[n][base], sm[base].astype(float)), 2.5),
                             hi=np.nanpercentile(B.stat(ys, S[n][base], sm[base].astype(float)), 97.5)))
    # missing measurements
    ref = S["A"]
    for frac in [0.05, 0.10, 0.20]:
        cond_ap, cond_cov, def_cov, ref_ap, boot_c, boot_r = [], [], [], [], [], []
        for seed in [1, 2, 3]:
            r = mask_res[(frac, seed)]
            use = sets[sn] & r["affected"] & np.isfinite(r["cond"]) & np.isfinite(ref)
            if y[use].sum() < 3:
                continue
            cond_ap.append(bench.ap(y[use], r["cond"][use])); ref_ap.append(bench.ap(y[use], ref[use]))
            cond_cov.append(np.isfinite(r["cond"][sets[sn]]).mean()); def_cov.append(np.isfinite(r["default"][sets[sn]]).mean())
            w = use[base].astype(float)
            boot_c.append(B.stat(ys, r["cond"][base], w)); boot_r.append(B.stat(ys, ref[base], w))
        vc = np.mean(boot_c, axis=0); vr = np.mean(boot_r, axis=0)
        rows.append(dict(set=sn, condition=f"mask{int(frac*100)}", auprc=np.mean(cond_ap), lo=np.nanpercentile(vc, 2.5),
                         hi=np.nanpercentile(vc, 97.5), ref_auprc=np.mean(ref_ap), ref_lo=np.nanpercentile(vr, 2.5),
                         ref_hi=np.nanpercentile(vr, 97.5), cov_default=np.mean(def_cov), cov_cond=np.mean(cond_cov),
                         cov_default_sd=np.std(def_cov, ddof=1), cov_cond_sd=np.std(cond_cov, ddof=1),
                         n_common=int(use.sum()), pos_common=int(y[use].sum())))
R = pd.DataFrame(rows)
R.to_csv("results/k562_degradation.csv", index=False)
pd.set_option("display.width", 250)
print(R.round(3).to_string())
print("gamma by condition:", {k: round(v, 3) for k, v in meta.items()})
print("done %.0fs" % (time.time() - t0))
