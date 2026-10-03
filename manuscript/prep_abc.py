"""Reduce the ENCODE K562 ABC all-putative predictions (ENCFF681DDZ) to the benchmark genes."""
import sys; sys.path.insert(0, ".")
import pandas as pd
from pacefig import bench
genes = set(bench.load_benchmark().gene_id)
cols = ["#chr", "start", "end", "class", "TargetGene", "TargetGeneEnsembl_ID", "distance", "isSelfPromoter", "ABC.Score", "powerlaw.Score"]
out = []
for ch in pd.read_csv("/home/user/data/k562/e2g/ENCFF681DDZ.bed.gz", sep="\t", usecols=cols, chunksize=2_000_000):
    out.append(ch[ch.TargetGeneEnsembl_ID.isin(genes)])
d = pd.concat(out).rename(columns={"#chr": "chrom", "TargetGeneEnsembl_ID": "gene_id"})
d.to_parquet(bench.WORK + "abc_encode_k562.parquet") if False else d.to_csv(bench.WORK + "abc_encode_k562.tsv.gz", sep="\t", index=False)
print(len(d), d.gene_id.nunique()); print(d.describe().T)
