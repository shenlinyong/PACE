import sys,time; sys.path.insert(0,".")
import numpy as np
from pacefig import hicband, engine
t=time.time()
for res,maxd in [(10000,5_000_000),(25000,5_000_000),(5000,5_000_000),(40000,5_000_000)]:
    if res==40000: continue
    out=f"/home/user/work/gm12878/band/r{res}"
    hicband.extract_band_hic("/home/user/data/gm12878/4DNFI1UEG1HD.hic",res,out,maxd=maxd)
    for ch in [f"chr{i}" for i in range(1,23)]+["chrX"]:
        z=np.load(f"{out}/{ch}.npz"); b=z["band"]
        np.savez_compressed(f"{out}/{ch}.npz",band=b,weight=engine.ice_balance(b))
    print(res,time.time()-t,flush=True)
