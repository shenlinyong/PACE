"""Catalog tables -> per-chromosome arrays for the vectorised engine."""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import engine


class Catalog:
    def __init__(self, catalog_dir: str | None = None, genes: set | None = None, radius: int = 5_000_000,
                 units: pd.DataFrame | None = None, promoters: pd.DataFrame | None = None):
        if units is None:
            units = pd.read_csv(f"{catalog_dir}/units.tsv", sep="\t", usecols=["element_id", "chrom", "start", "end", "anchor0"])
            promoters = pd.read_csv(f"{catalog_dir}/promoters.tsv", sep="\t")
        self.units = units.reset_index(drop=True)
        p = promoters
        if genes is not None:
            p = p[p.gene_id.isin(genes)]
        self.promoters = p.reset_index(drop=True)
        self.radius = radius
        self.unit_index = pd.Index(self.units.element_id)
        self._cand = None

    # candidate pairs: every element within `radius` of any TSS of the gene
    def candidates(self, radius=None, tss_mode="all"):
        radius = self.radius if radius is None else radius
        out = {}
        for ch, pg in self.promoters.groupby("chrom"):
            uu = self.units[self.units.chrom == ch]
            order = np.argsort(uu.anchor0.to_numpy(), kind="stable")
            uidx = uu.index.to_numpy()[order]
            anc = uu.anchor0.to_numpy()[order]
            gl, cu, cg = [], [], []
            for gi, (g, t) in enumerate(pg.groupby("gene_id", sort=True)):
                pos = t.tss0.to_numpy()
                lo = np.searchsorted(anc, pos.min() - radius, "left")
                hi = np.searchsorted(anc, pos.max() + radius, "right")
                seg = anc[lo:hi]
                dmin = np.abs(seg[:, None] - pos[None, :]).min(axis=1)
                keep = dmin <= radius
                sel = uidx[lo:hi][keep]
                cu.append(sel)
                cg.append(np.full(len(sel), gi))
                gl.append(g)
            out[ch] = (np.concatenate(cu), np.concatenate(cg), gl)
        return out


def tss_table(promoters: pd.DataFrame, mode="all", canonical: dict | None = None):
    """TSS rows per gene (gene_idx within chromosome) with weights pi (equal by default)."""
    rows = []
    for ch, pg in promoters.groupby("chrom"):
        genes = sorted(pg.gene_id.unique())
        gi = {g: i for i, g in enumerate(genes)}
        t = pg.copy()
        if mode == "canonical":
            t = t[t.promoter_id.isin(set(canonical.values()))]
        t["gene_idx"] = t.gene_id.map(gi)
        t["pi"] = 1.0 / t.groupby("gene_id")["tss0"].transform("size")
        rows.append((ch, genes, t.sort_values(["gene_idx", "tss0"]).reset_index(drop=True)))
    return {ch: (genes, t) for ch, genes, t in rows}
