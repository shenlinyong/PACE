"""Check that the vectorised engine reproduces `pace run` (Case A on a gene subset)."""
import sys, json, numpy as np, pandas as pd
sys.path.insert(0, ".")
from pacefig import engine, hicband
from pacefig.dataset import Catalog, tss_table
W = "/home/user/work/k562/"
genes = open(W + "sub_genes.txt").read().split()
cat = Catalog(W + "catalog_sub", genes=set(genes))
man = json.load(open(W + "priors/k562_10k/manifest.json"))
prior = engine.Prior(gamma=man["gamma"], a=man["a"], d_ref=man["d_ref"], d_min=man["d_min"], res=man["resolution"])
obs = pd.read_csv(W + "act/all/observed_activity.tsv", sep="\t")
A = engine.element_activity(obs, {"DNase": ["dnase_r1", "dnase_r2"], "H3K27ac": ["k27ac_r1", "k27ac_r2"]}, cat.unit_index)
cands = cat.candidates()
tt = tss_table(cat.promoters)
res_rows = []
for ch, (cu, cg, gl) in cands.items():
    band, w = hicband.load_band(W + "band/r10000", ch)
    genes_ch, t = tt[ch]
    cbar, s_w, s_c = engine.score_chromosome(cat.units, A, cu, cg, pd.DataFrame(index=range(len(gl))), t, "A", prior,
                                             band=band, weight=w, res=10000)
    res_rows.append(pd.DataFrame(dict(element_id=cat.units.element_id.to_numpy()[cu], gene_id=np.array(gl)[cg],
                                      Cbar=cbar, score=s_w, cond=s_c)))
mine = pd.concat(res_rows)
ref = pd.read_csv(W + "run_sub/scores.tsv.gz", sep="\t")
print(ref.columns.tolist())
m = mine.merge(ref[["element_id", "gene_id", "pace_score", "pace_score_conditional"]], on=["element_id", "gene_id"], how="outer", indicator=True)
print(m._merge.value_counts().to_dict())
both = m[m._merge == "both"]
d = np.abs(both.score - both.pace_score)
print("max |Δscore| =", np.nanmax(d), " n_nan engine/CLI:", both.score.isna().sum(), both.pace_score.isna().sum())
print("spearman", both[["score", "pace_score"]].corr("spearman").iloc[0, 1])
m2 = mine.merge(ref[["element_id","gene_id","Cbar","A_used","distance_bp","n_tss"]], on=["element_id","gene_id"], suffixes=("","_cli"))
m2["rel"] = np.abs(m2.Cbar-m2.Cbar_cli)/m2.Cbar_cli
print(m2.rel.describe())
print(m2.sort_values("rel",ascending=False).head(8)[["element_id","gene_id","Cbar","Cbar_cli","distance_bp","rel"]])
a = pd.Series(A, index=cat.unit_index)
m2["A"]=a.reindex(m2.element_id).to_numpy()
print("A diff", np.nanmax(np.abs(m2.A-m2.A_used)))
# ---- Cases B and C
for case, d, pr in [("B", "run_subB", prior), ("C", "run_subC", engine.human_reference_prior())]:
    rows = []
    for ch, (cu, cg, gl) in cands.items():
        genes_ch, t = tt[ch]
        _, s_w, _ = engine.score_chromosome(cat.units, A, cu, cg, pd.DataFrame(index=range(len(gl))), t, case, pr)
        rows.append(pd.DataFrame(dict(element_id=cat.units.element_id.to_numpy()[cu], gene_id=np.array(gl)[cg], score=s_w)))
    mm = pd.concat(rows).merge(pd.read_csv(W + d + "/scores.tsv.gz", sep="\t")[["element_id","gene_id","pace_score"]], on=["element_id","gene_id"])
    print("Case", case, "n", len(mm), "max |Δscore|", np.nanmax(np.abs(mm.score - mm.pace_score)))
