"""Own catalog builder (same rules as ``pace catalog``) so that grid width and expression threshold can be varied."""
from __future__ import annotations

import gzip
import re

import numpy as np
import pandas as pd
import pyBigWig

CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]


def promoters_from_gtf(gtf_gz: str, genes: set) -> pd.DataFrame:
    """Distinct physical TSS per gene (transcripts sharing a TSS counted once), as promoters.tsv."""
    rows = []
    op = gzip.open if gtf_gz.endswith(".gz") else open
    with op(gtf_gz, "rt") as f:
        for l in f:
            if l[0] == "#":
                continue
            x = l.split("\t")
            chrom = x[0] if x[0].startswith("chr") else "chr" + x[0]
            if x[2] != "transcript" or not re.fullmatch(r"chr([0-9]+|X)", chrom):
                continue
            g = re.search(r'gene_id "([^"]+)"', x[8]).group(1)
            if g not in genes:
                continue
            strand = x[6]
            tss1 = int(x[3]) if strand == "+" else int(x[4])
            rows.append((g, chrom, tss1 - 1, strand))
    d = pd.DataFrame(rows, columns=["gene_id", "chrom", "tss0", "strand"]).drop_duplicates()
    d["promoter_id"] = d.chrom + ":" + d.tss0.astype(str) + ":" + d.strand
    return d.sort_values(["gene_id", "tss0"]).reset_index(drop=True)


def build_units(peaks: pd.DataFrame, promoters: pd.DataFrame, chrom_sizes: dict, width: int) -> pd.DataFrame:
    cells = []
    for ch, p in peaks.groupby("chrom"):
        lo = p.start.to_numpy() // width
        hi = (p.end.to_numpy() - 1) // width
        n = hi - lo + 1
        idx = np.repeat(lo, n) + (np.arange(n.sum()) - np.repeat(np.cumsum(n) - n, n))
        cells.append(pd.DataFrame({"chrom": ch, "cell": idx, "peak": True}))
    pr = promoters.assign(cell=promoters.tss0 // width, peak=False)[["chrom", "cell", "peak"]]
    allc = pd.concat(cells + [pr], ignore_index=True)
    allc["promoter"] = False
    allc.loc[allc.index >= len(allc) - len(pr), "promoter"] = True
    u = allc.groupby(["chrom", "cell"], as_index=False).agg(promoter=("promoter", "max"))
    u["start"] = u.cell * width
    u["end"] = np.minimum(u.start + width, u.chrom.map(chrom_sizes))
    u = u[u.start < u.end]
    u["anchor0"] = (u.start + u.end - 1) // 2
    u["element_id"] = u.chrom + ":" + u.start.astype(str) + "-" + u.end.astype(str)
    return u.sort_values(["chrom", "start"]).reset_index(drop=True)[["element_id", "chrom", "start", "end", "anchor0", "promoter"]]


def bigwig_mean(bw_path: str, units: pd.DataFrame) -> np.ndarray:
    """Mean signal per cell over covered bases (NaN where the bigWig has no data in the cell)."""
    bw = pyBigWig.open(bw_path)
    out = np.full(len(units), np.nan)
    chroms = bw.chroms()
    ch_arr = units.chrom.to_numpy(); s_arr = units.start.to_numpy(); e_arr = units.end.to_numpy()
    for k in range(len(units)):
        ch = ch_arr[k]
        if ch not in chroms:
            continue
        v = bw.stats(ch, int(s_arr[k]), int(min(e_arr[k], chroms[ch])), type="mean", exact=True)[0]
        if v is not None:
            out[k] = v
    bw.close()
    return out
