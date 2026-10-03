"""Component analysis and sensitivity to the fixed settings (Fig. 2e,f; Supplementary Fig. S5, Table S11)."""
import sys, time, pickle
sys.path.insert(0, ".")
import numpy as np, pandas as pd
from pacefig import k562, bench, catalog, engine

W = k562.W
t0 = time.time()
base = k562.K562()
bench_genes = set(base.pairs.gene_id)
bs = k562.BandSet.load(W + "band/r10000", 10000)
prior = bs.fit()
y = base.pairs.y.to_numpy()
sets = {"main": (base.pairs.set == "main").to_numpy(), "heldout": (base.pairs.set == "heldout").to_numpy()}
boots = {k: bench.Boot(base.pairs.gene_id.to_numpy()[m]) for k, m in sets.items()}

rows, scores = [], {}
def record(tag, param, value, s):
    scores[(param, value)] = s
    for sn, m in sets.items():
        ys = y[m]; ss = s[m]
        mt = bench.metrics(ys, ss)
        bt = boots[sn].stat(ys, ss)
        rows.append(dict(parameter=param, value=value, set=sn, auprc=mt["auprc"], lo=np.nanpercentile(bt, 2.5),
                         hi=np.nanpercentile(bt, 97.5), coverage=mt["coverage"]))

record("default", "default", "default", base.pair_scores(base.run("A", prior, bands=bs)))
for R in [1_000_000, 2_000_000, 10_000_000]:
    record("window", "candidate window (Mb)", R / 1e6, base.pair_scores(base.run("A", prior, bands=bs, radius=R)))
record("window", "candidate window (Mb)", 5.0, scores[("default", "default")])
for d0 in [2000.0, 10000.0]:
    record("d0", "pseudocount distance d0 (kb)", d0 / 1e3, base.pair_scores(base.run("A", prior, bands=bs, d0=d0)))
record("d0", "pseudocount distance d0 (kb)", 5.0, scores[("default", "default")])
print("window/d0 done %.0fs" % (time.time() - t0), flush=True)

# configurations that change the candidate grid or gene set: rebuild catalog and activity from the bigWigs
chrom_sizes = pd.read_csv("/home/user/data/ref/hg38.chrom.sizes", sep="\t", header=None, index_col=0)[1].to_dict()
peaks = pd.read_csv(W + "dnase_peaks_union.bed", sep="\t", header=None, names=["chrom", "start", "end"])
tpm = pd.read_csv(W + "k562_tpm.tsv", sep="\t", index_col=0)["tpm"]
BW = {"dnase_r1": "ENCFF414OGC", "dnase_r2": "ENCFF338LXW", "k27ac_r1": "ENCFF177JXC", "k27ac_r2": "ENCFF745TSK"}
gtf = "/home/user/data/ref/Homo_sapiens.GRCh38.110.gtf.gz"

def workspace(width, thr):
    genes = set(tpm[tpm >= thr].index)
    pr = catalog.promoters_from_gtf(gtf, genes)
    units = catalog.build_units(peaks, pr, chrom_sizes, width)
    obs = []
    for sid, f in BW.items():
        v = catalog.bigwig_mean(f"/home/user/data/k562/{f}.bigWig", units)
        obs.append(pd.DataFrame(dict(element_id=units.element_id, signal=v, sample_id=sid,
                                     measurement_status=np.where(np.isfinite(v), "observed", "unmeasured"))))
    cat = k562.Catalog(units=units, promoters=pr, genes=bench_genes & genes)
    K = k562.K562(cat=cat, obs=pd.concat(obs), cell=width)
    return K

for width in [250, 1000]:
    K = workspace(width, 1.0)
    # element width changes the 10 kb binned contact lookup only through the element anchor
    record("width", "element width (bp)", width, K.pair_scores(K.run("A", prior, bands=bs)))
    print("width", width, "%.0fs" % (time.time() - t0), flush=True)
record("width", "element width (bp)", 500, scores[("default", "default")])
for thr in [0.5, 2.0, 5.0]:
    K = workspace(500, thr)
    record("tpm", "expression threshold (TPM)", thr, K.pair_scores(K.run("A", prior, bands=bs)))
    print("tpm", thr, "%.0fs" % (time.time() - t0), flush=True)
record("tpm", "expression threshold (TPM)", 1.0, scores[("default", "default")])

R = pd.DataFrame(rows)
# Jaccard overlap of the top 10% of benchmark-pair predictions with the default configuration
d_s = scores[("default", "default")]
jac = {}
for (param, value), s in scores.items():
    out = {}
    for sn, m in sets.items():
        a = d_s[m]; b = s[m]
        ka = np.isfinite(a); kb = np.isfinite(b)
        na = max(1, int(round(0.1 * ka.sum()))); nb = max(1, int(round(0.1 * kb.sum())))
        ia = set(np.flatnonzero(m)[np.argsort(-np.where(ka, a, -1))[:na]])
        ib = set(np.flatnonzero(m)[np.argsort(-np.where(kb, b, -1))[:nb]])
        out[sn] = len(ia & ib) / len(ia | ib)
    jac[(param, value)] = out
R["jaccard_top10"] = [jac[(p, v)][s] if (p, v) in jac else (jac[("default", "default")][s] if p == "default" else np.nan)
                      for p, v, s in zip(R.parameter, R.value, R.set)]
R.to_csv("results/k562_sensitivity.csv", index=False)
print(R.round(3).to_string())
pickle.dump({k: v for k, v in scores.items()}, open("results/k562_sensitivity_scores.pkl", "wb"))
print("done %.0fs" % (time.time() - t0))
