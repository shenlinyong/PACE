"""Banded Hi-C storage: band[b, k] = raw count between bins b and b+k (k = 0..K) on one chromosome."""
import numpy as np, os, cooler

MAXD = 5_000_000

def extract_band(uri, res, outdir, maxd=MAXD):
    """Stream the cis pixels of a cooler within maxd into dense per-chromosome arrays."""
    c = cooler.Cooler(uri)
    assert c.binsize == res
    os.makedirs(outdir, exist_ok=True)
    K = maxd // res + 2
    bins = c.bins()[:]
    w = bins["weight"].to_numpy(float) if "weight" in bins else np.full(len(bins), np.nan)
    chrom_of = {}
    offs = {ch: c.extent(ch) for ch in c.chromnames}
    bands = {ch: np.zeros((offs[ch][1] - offs[ch][0], K), dtype=np.int32) for ch in c.chromnames
             if ch != "chrY" and ch != "chrM"}
    starts = np.array([offs[ch][0] for ch in c.chromnames]); names = np.array(c.chromnames)
    nnz = int(c.info["nnz"]); step = 2_000_000
    for s in range(0, nnz, step):
        p = c.pixels()[s:s + step]
        i = p["bin1_id"].to_numpy(np.int64); j = p["bin2_id"].to_numpy(np.int64); v = p["count"].to_numpy()
        d = j - i
        keep = d < K
        i, j, v, d = i[keep], j[keep], v[keep], d[keep]
        ci = np.searchsorted(starts, i, side="right") - 1
        cj = np.searchsorted(starts, j, side="right") - 1
        m = ci == cj
        i, d, v, ci = i[m], d[m], v[m], ci[m]
        for code in np.unique(ci):
            ch = names[code]
            if ch not in bands: continue
            sel = ci == code
            bands[ch][i[sel] - offs[ch][0], d[sel]] += v[sel].astype(np.int32)
    for ch, b in bands.items():
        s, e = offs[ch]
        np.savez(os.path.join(outdir, f"{ch}.npz"), band=b, weight=w[s:e])
    return list(bands)

def load_band(outdir, chrom):
    z = np.load(os.path.join(outdir, f"{chrom}.npz"))
    return z["band"], z["weight"]


def coarsen_band(band: np.ndarray, f: int, maxbins: int):
    """Aggregate a banded raw-count array by an integer factor (upper-triangle pairs counted once)."""
    n, K = band.shape
    nc = (n + f - 1) // f
    Kc = min(maxbins, K // f)
    out = np.zeros((nc, Kc), dtype=np.int64)
    pad = np.zeros((nc * f, K + f), dtype=np.int64)
    pad[:n, :K] = band
    for a in range(f):
        for b in range(f):
            for Kk in range(Kc):
                k = f * Kk + b - a
                if k < 0 or (Kk == 0 and k < 0):
                    continue
                if k >= K:
                    continue
                out[:, Kk] += pad[a::f, k][:nc]
    return out.astype(np.int32)


def extract_band_hic(path, res, outdir, maxd=MAXD, norm="KR", chroms=None, stripe=30_000_000):
    """Banded raw counts and normalisation weights from a .hic file (via hictkpy)."""
    import hictkpy
    f = hictkpy.File(path, res)
    os.makedirs(outdir, exist_ok=True)
    K = maxd // res + 2
    for name, size in f.chromosomes().items():
        if name in ("All", "ALL") or name.lstrip("chr") in ("Y", "M", "MT"):
            continue
        ch = name if name.startswith("chr") else "chr" + name
        if chroms and ch not in chroms:
            continue
        n = (size + res - 1) // res
        band = np.zeros((n, K), dtype=np.int32)
        for a in range(0, size, stripe):
            b = min(size, a + stripe)
            m = f.fetch(f"{name}:{a}-{b}", f"{name}:{a}-{min(size, b + maxd + 2 * res)}", normalization="NONE").to_coo().tocoo()
            r = m.row.astype(np.int64) + a // res
            c = m.col.astype(np.int64) + a // res
            lo = np.minimum(r, c); hi = np.maximum(r, c)
            d = hi - lo
            keep = (d < K) & (lo >= a // res) & (lo < b // res)
            np.add.at(band, (lo[keep], d[keep]), m.data[keep].astype(np.int32)) if False else None
            idx = lo[keep] * K + d[keep]
            band.reshape(-1)[:] += np.bincount(idx, weights=m.data[keep], minlength=n * K).astype(np.int32)
        w = np.asarray(f.weights(norm) if False else np.full(n, np.nan))
        try:
            w = f.weights(norm, name)
            w = np.asarray(w, dtype=float)
        except Exception:
            pass
        np.savez_compressed(os.path.join(outdir, f"{ch}.npz"), band=band, weight=w)
