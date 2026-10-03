"""Figure 2 and Supplementary Fig. S2: PACE and comparison methods on the K562 CRISPR benchmark."""
import sys, pickle, time
sys.path.insert(0, ".")
import numpy as np, pandas as pd
from sklearn.metrics import precision_recall_curve
from pacefig import k562, bench, external, engine

t0 = time.time()
K = k562.K562()
P = K.pairs.copy()
y = P.y.to_numpy()
bs = k562.BandSet.load(k562.W + "band/r10000", 10000)
prior = bs.fit()
print("K562 prior from banded Hi-C: gamma=%.4f a=%.4f" % (prior.gamma, prior.a))

methods = {}
def pace(tag, **kw):
    tab = K.run(**kw)
    methods[tag] = K.pair_scores(tab)
    return tab

tab_A = pace("PACE", mode="A", prior=prior, bands=bs)
cbar_pairs = K.pair_scores(tab_A, "cbar")
# component alternatives (Fig 2e, f)
pace("PACE, accessibility only", mode="A", prior=prior, bands=bs, A=K.activity("dnase"))
pace("PACE, H3K27ac only", mode="A", prior=prior, bands=bs, A=K.activity("k27ac"))
pace("PACE, single canonical TSS", mode="A", prior=prior, bands=bs, tss="canonical")

# external comparisons
abc = external.load_abc()
m = external.map_intervals(P, abc, ("ABC.Score", "powerlaw.Score"))
methods["ABC (ENCODE K562 predictions)"] = m["ABC.Score"]
methods["ABC power-law contact (ENCODE)"] = m["powerlaw.Score"]
dist = P.distanceToTSS.to_numpy().astype(float)
methods["Inverse distance"] = 1.0 / dist
methods["Distance prior P_s(d)"] = prior(np.round(dist / 10000) * 10000)
# supervised: TargetFinder-style gradient boosting (re-implementation), trained on the benchmark
X = external.gbm_features(P, {"pace_contact": cbar_pairs})
methods["GBM (TargetFinder-style)"] = external.gbm_scores(P, X)
print("scored in %.0fs" % (time.time() - t0))

sets = {"main": (P.set == "main").to_numpy(), "heldout": (P.set == "heldout").to_numpy()}
res = {"methods": methods, "pairs": P.drop(columns=[c for c in P.columns if c.startswith("direct_")]), "prior_gamma": prior.gamma}
boots = {k: bench.Boot(P.gene_id.to_numpy()[m]) for k, m in sets.items()}
table = []
curves = {}
for sname, mask in sets.items():
    B = boots[sname]
    ys = y[mask]
    for name, s in methods.items():
        ss = s[mask]
        mt = bench.metrics(ys, ss)
        bt = B.stat(ys, ss)
        lo, hi = B.ci(bt)
        # AUROC bootstrap
        from sklearn.metrics import roc_auc_score
        s2 = np.where(np.isfinite(ss), ss, np.nanmin(ss) - 1)
        ar = np.array([roc_auc_score(ys, s2, sample_weight=w) for w in B.W[:300]])
        table.append(dict(set=sname, method=name, auprc=mt["auprc"], auprc_lo=lo, auprc_hi=hi, auroc=mt["auroc"],
                          auroc_lo=np.percentile(ar, 2.5), auroc_hi=np.percentile(ar, 97.5), coverage=mt["coverage"],
                          cov_pos=mt["cov_pos"], n=int(mask.sum()), n_pos=int(ys.sum())))
        pr, rc, _ = precision_recall_curve(ys, s2)
        curves[(sname, name)] = (rc, pr)
T = pd.DataFrame(table)
T.to_csv("results/k562_benchmark_metrics.csv", index=False)
res["curves"] = curves
print(T.round(3).to_string())

# paired differences on the held-out set, restricted to pairs scored by both (Fig 2c)
diff = []
for sname in ["heldout", "main"]:
    mask = sets[sname]; B = boots[sname]; ys = y[mask]; sp = methods["PACE"][mask]
    for name, s in methods.items():
        if name == "PACE" or name.startswith("PACE,"):
            continue
        so = s[mask]
        both = np.isfinite(sp) & np.isfinite(so)
        a = B.stat(ys, sp, both.astype(float)); b = B.stat(ys, so, both.astype(float))
        d = a - b
        diff.append(dict(set=sname, method=name, n_common=int(both.sum()), pos_common=int(ys[both].sum()),
                         auprc_pace=bench.ap(ys[both], sp[both]), auprc_other=bench.ap(ys[both], so[both]),
                         delta=bench.ap(ys[both], sp[both]) - bench.ap(ys[both], so[both]),
                         lo=np.nanpercentile(d, 2.5), hi=np.nanpercentile(d, 97.5),
                         p_boot=2 * min(np.nanmean(d <= 0), np.nanmean(d >= 0))))
D = pd.DataFrame(diff); D.to_csv("results/k562_benchmark_paired_differences.csv", index=False)
print(D.round(3).to_string())

# distance strata and distance-matched negatives (Supplementary Fig. S2)
strata = [("<50 kb", 0, 50e3), ("50-200 kb", 50e3, 200e3), (">200 kb", 200e3, 1e12)]
rows = []
show = ["PACE", "ABC (ENCODE K562 predictions)", "ABC power-law contact (ENCODE)", "Inverse distance", "GBM (TargetFinder-style)"]
for sname, mask in sets.items():
    B = boots[sname]
    for lab, lo_, hi_ in strata:
        sm = mask & (dist >= lo_) & (dist < hi_)
        ys = y[sm]
        if ys.sum() < 3:
            continue
        for name in show:
            s = methods[name][sm]
            bt = np.array([bench.ap(ys, np.where(np.isfinite(s), s, np.nanmin(s) - 1), w[sm[mask] if False else None]) for w in [None]]) if False else None
            # bootstrap with the mask-restricted replicate weights
            Wm = B.W[:, sm[mask]]
            s2 = np.where(np.isfinite(s), s, np.nanmin(s) - 1)
            bt = np.array([bench.ap(ys, s2, w) if (w * ys).sum() > 0 else np.nan for w in Wm])
            rows.append(dict(set=sname, stratum=lab, method=name, auprc=bench.ap(ys, s2), lo=np.nanpercentile(bt, 2.5),
                             hi=np.nanpercentile(bt, 97.5), n=int(sm.sum()), n_pos=int(ys.sum()), prevalence=ys.mean()))
    # distance-matched negatives: reweight negatives so their log-distance distribution equals that of positives
    ld = np.log10(dist[mask]); ys = y[mask]
    edges = np.quantile(ld[ys == 1], np.linspace(0, 1, 9)); edges[0] -= 1; edges[-1] += 1
    b = np.digitize(ld, edges[1:-1])
    wpos = np.bincount(b[ys == 1], minlength=8) / (ys == 1).sum()
    wneg = np.bincount(b[ys == 0], minlength=8) / (ys == 0).sum()
    ratio = np.where(wneg > 0, wpos / np.maximum(wneg, 1e-12), 0.0)
    w_match = np.where(ys == 1, 1.0, ratio[b] * ((ys == 0).sum() / (ys == 1).sum()) * 0 + ratio[b])
    # scale negatives so the weighted prevalence equals the original prevalence
    w_match[ys == 0] *= (ys == 0).sum() / w_match[ys == 0].sum()
    for name in show:
        s = methods[name][mask]; s2 = np.where(np.isfinite(s), s, np.nanmin(s) - 1)
        bt = np.array([bench.ap(ys, s2, w * w_match) if (w * ys).sum() > 0 else np.nan for w in B.W])
        rows.append(dict(set=sname, stratum="distance-matched negatives", method=name, auprc=bench.ap(ys, s2, w_match),
                         lo=np.nanpercentile(bt, 2.5), hi=np.nanpercentile(bt, 97.5), n=int(mask.sum()), n_pos=int(ys.sum()),
                         prevalence=float((ys * w_match).sum() / w_match.sum())))
S = pd.DataFrame(rows); S.to_csv("results/k562_benchmark_distance_strata.csv", index=False)
print(S.round(3).to_string())
pickle.dump(res, open("results/k562_benchmark.pkl", "wb"))
print("done in %.0fs" % (time.time() - t0))
