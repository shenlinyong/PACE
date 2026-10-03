"""K562 benchmark workspace: catalog, activity, Hi-C bands and a single ``run`` entry point."""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from . import bench, engine, hicband
from .dataset import Catalog, tss_table

W = bench.WORK
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
PANELS = {
    "both": {"DNase": ["dnase_r1", "dnase_r2"], "H3K27ac": ["k27ac_r1", "k27ac_r2"]},
    "dnase": {"DNase": ["dnase_r1", "dnase_r2"]},
    "k27ac": {"H3K27ac": ["k27ac_r1", "k27ac_r2"]},
    "both_rep1": {"DNase": ["dnase_r1"], "H3K27ac": ["k27ac_r1"]},
}


class BandSet:
    """All chromosomes of one Hi-C map: banded raw counts and balancing weights."""

    def __init__(self, bands, weights, res):
        self.bands, self.weights, self.res = bands, weights, res

    @classmethod
    def load(cls, d, res, chroms=CHROMS):
        b, w = {}, {}
        for ch in chroms:
            b[ch], w[ch] = hicband.load_band(d, ch)
        return cls(b, w, res)

    def fit(self):
        return engine.fit_prior({c: self.bands[c][:, : 5_000_000 // self.res + 1] for c in self.bands}, self.weights, self.res)[0]

    def thinned(self, frac, seed):
        rng = np.random.default_rng(seed)
        nb = {c: rng.binomial(b, frac).astype(np.int32) for c, b in self.bands.items()}
        return BandSet(nb, {c: engine.ice_balance(b, maxd_bins=5_000_000 // self.res + 2) for c, b in nb.items()}, self.res)

    def rebalanced(self):
        return BandSet(self.bands, {c: engine.ice_balance(b, maxd_bins=5_000_000 // self.res + 2) for c, b in self.bands.items()}, self.res)


class K562:
    def __init__(self, load_bench=True, cat=None, obs=None, cell=500):
        self.cat = cat if cat is not None else Catalog(W + "catalog_bench")
        self.obs = obs if obs is not None else pd.read_csv(W + "act/all/observed_activity.tsv", sep="\t")
        self.canon = pd.read_csv(W + "canonical_promoter.tsv", sep="\t", index_col=0)["promoter_id"].to_dict()
        self._cand = {}
        self._tss = {}
        self._A = {}
        if load_bench:
            self.pairs = bench.load_benchmark()
            self.pm = bench.PairMap(self.pairs, self.cat.unit_index, cell)

    # ---- inputs
    def activity(self, panel="both"):
        if panel not in self._A:
            self._A[panel] = engine.element_activity(self.obs, PANELS[panel], self.cat.unit_index)
        return self._A[panel]

    def cands(self, radius=5_000_000):
        if radius not in self._cand:
            self._cand[radius] = self.cat.candidates(radius)
        return self._cand[radius]

    def tss(self, mode="all"):
        if mode not in self._tss:
            self._tss[mode] = tss_table(self.cat.promoters, mode, self.canon)
        return self._tss[mode]

    def prior_k562(self):
        m = json.load(open(W + "priors/k562_10k/manifest.json"))
        return engine.Prior(gamma=m["gamma"], a=m["a"], d_ref=m["d_ref"], d_min=m["d_min"], res=m["resolution"])

    # ---- scoring
    def run(self, mode, prior, A=None, bands: BandSet | None = None, radius=5_000_000, tss="all",
            kappa=1.0, d0=5000.0, fallback=True, A_override=None):
        A = self.activity("both") if A is None else A
        cands = self.cands(radius)
        tt = self.tss(tss)
        out = []
        for ch, (cu, cg, gl) in cands.items():
            if ch not in tt:
                continue
            genes_ch, t = tt[ch]
            # genes of the candidate table index into the same sorted gene list as tss_table
            assert list(genes_ch) == list(gl)
            band = bands.bands[ch] if bands is not None else None
            w = bands.weights[ch] if bands is not None else None
            res = bands.res if bands is not None else None
            cbar, s_w, s_c = engine.score_chromosome(self.cat.units, A, cu, cg, pd.DataFrame(index=range(len(gl))), t, mode,
                                                     prior, band=band, weight=w, res=res, kappa=kappa, d0=d0, fallback=fallback)
            out.append(pd.DataFrame(dict(unit_idx=cu, gene_id=np.asarray(gl)[cg], cbar=cbar, score=s_w, cond=s_c)))
        return pd.concat(out, ignore_index=True)

    def pair_scores(self, table, col="score"):
        return self.pm.lookup(table, col)
