"""CRISPR benchmark handling: mapping PACE scores to tested pairs, metrics, gene-level bootstrap."""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

from . import engine
from .dataset import Catalog, tss_table

CRISPR = "/home/user/engreitzlab/crispr_comparison/resources/crispr_data/"
WORK = "/home/user/work/k562/"
CELL = 500


def load_benchmark():
    a = pd.read_csv(CRISPR + "EPCrisprBenchmark_combined_data.training_K562.GRCh38.tsv.gz", sep="\t")
    h = pd.read_csv(CRISPR + "EPCrisprBenchmark_combined_data.heldout_5_cell_types.GRCh38.tsv.gz", sep="\t")
    h = h[h.CellType == "K562"]
    a = a.assign(set="main")
    h = h.assign(set="heldout")
    d = pd.concat([a, h], ignore_index=True)
    d["y"] = d.Regulated.astype(int)
    d["gene_id"] = d.measuredGeneEnsemblId
    return d


class PairMap:
    """Maps benchmark element-gene pairs onto the 500 bp grid cells of the catalog."""

    def __init__(self, pairs: pd.DataFrame, unit_index: pd.Index, CELL: int = CELL):
        self.n = len(pairs)
        s = (pairs.chromStart.to_numpy() // CELL) * CELL
        e = ((pairs.chromEnd.to_numpy() - 1) // CELL) * CELL
        ncell = ((e - s) // CELL + 1).astype(int)
        pi = np.repeat(np.arange(self.n), ncell)
        off = np.arange(len(pi)) - np.repeat(np.cumsum(ncell) - ncell, ncell)
        start = s[pi] + off * CELL
        ids = pairs.chrom.to_numpy()[pi].astype(str)
        eid = np.char.add(np.char.add(np.char.add(ids, ":"), start.astype(str)), np.char.add("-", (start + CELL).astype(str)))
        ui = unit_index.get_indexer(eid)
        keep = ui >= 0
        self.pair_idx = pi[keep]
        self.unit_idx = ui[keep]
        self.gene = pairs.gene_id.to_numpy()[self.pair_idx]
        self.gene_codes, self.gene_uniques = pd.factorize(self.gene)

    def lookup(self, table: pd.DataFrame, col: str = "score") -> np.ndarray:
        """Max of ``col`` over the cells of every pair; NaN where no cell has a finite score."""
        key_t = pd.MultiIndex.from_arrays([table.unit_idx.to_numpy(), table.gene_id.to_numpy()])
        key_p = pd.MultiIndex.from_arrays([self.unit_idx, self.gene])
        pos = key_t.get_indexer(key_p)
        val = np.where(pos >= 0, table[col].to_numpy()[np.maximum(pos, 0)], np.nan)
        out = pd.Series(val).groupby(self.pair_idx).max()
        res = np.full(self.n, np.nan)
        res[out.index.to_numpy()] = out.to_numpy()
        return res


def ap(y, s, w=None):
    return average_precision_score(y, s, sample_weight=w)


def metrics(y, s):
    """AUPRC / AUROC with unscored pairs (NaN) assigned the lowest score."""
    s = np.asarray(s, float)
    low = np.nanmin(s) - 1.0 if np.isfinite(s).any() else 0.0
    s2 = np.where(np.isfinite(s), s, low)
    return dict(auprc=ap(y, s2), auroc=roc_auc_score(y, s2), coverage=np.isfinite(s).mean(),
                cov_pos=np.isfinite(s[y == 1]).mean())


class Boot:
    """Gene-level bootstrap with a shared resampling design across methods (paired comparisons)."""

    def __init__(self, genes: np.ndarray, n: int = 1000, seed: int = 11):
        self.codes, uniq = pd.factorize(genes)
        rng = np.random.default_rng(seed)
        G = len(uniq)
        self.W = rng.multinomial(G, np.full(G, 1.0 / G), size=n)[:, self.codes]  # n x pairs replicate weights

    def stat(self, y, s, mask=None):
        s = np.asarray(s, float)
        low = np.nanmin(s) - 1.0
        s2 = np.where(np.isfinite(s), s, low)
        out = np.empty(len(self.W))
        for b, w in enumerate(self.W):
            if mask is not None:
                w = w * mask
            out[b] = ap(y, s2, w) if (w * y).sum() > 0 else np.nan
        return out

    @staticmethod
    def ci(x):
        return np.nanpercentile(x, [2.5, 97.5])
