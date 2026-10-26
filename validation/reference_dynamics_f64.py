"""scVelo's dynamical model on float64 moments (to measure float32 sensitivity)."""
import json, os, sys
import numpy as np
import scvelo as scv

data_dir, out_dir = sys.argv[1], sys.argv[2]
os.makedirs(out_dir, exist_ok=True)
scv.settings.verbosity = 0
adata = scv.datasets.pancreas(os.path.join(data_dir, "endocrinogenesis_day15.h5ad"))
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
scv.tl.velocity(adata, mode="stochastic")
adata.layers["Ms"] = adata.layers["Ms"].astype(np.float64)
adata.layers["Mu"] = adata.layers["Mu"].astype(np.float64)
scv.tl.recover_dynamics(adata, n_jobs=8, t_max=False)
from scvelo.tools.dynamical_model import align_dynamics
align_dynamics(adata)
shapes = {}
for key in ["alpha", "beta", "gamma", "t_", "scaling", "likelihood"]:
    a = np.ascontiguousarray(adata.var[f"fit_{key}"].values, dtype=np.float64)
    a.tofile(os.path.join(out_dir, f"fit_{key}.bin"))
    shapes[f"fit_{key}"] = list(a.shape)
a = np.ascontiguousarray(adata.layers["fit_t"], dtype=np.float64)
a.tofile(os.path.join(out_dir, "fit_t.bin"))
shapes["fit_t"] = list(a.shape)
json.dump(shapes, open(os.path.join(out_dir, "shapes.json"), "w"))
print("done")
