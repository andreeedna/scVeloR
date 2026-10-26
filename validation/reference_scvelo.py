"""Run the reference scVelo pipeline and export inputs/outputs for validation.

Usage: python -I validation/reference_scvelo.py <data_dir> <out_dir> [dynamical]

Arrays are written as raw little-endian float64 (C order) with shapes in
shapes.json, so the Rust/R side can read them without extra dependencies.
"""
import json
import os
import sys
import time

import numpy as np
import scvelo as scv
from scipy.sparse import issparse

data_dir, out_dir = sys.argv[1], sys.argv[2]
run_dynamical = len(sys.argv) > 3 and sys.argv[3] == "dynamical"
os.makedirs(out_dir, exist_ok=True)
shapes = {}
timings = {}


def dump(name, arr):
    arr = arr.toarray() if issparse(arr) else np.asarray(arr)
    arr = np.ascontiguousarray(arr, dtype=np.float64)
    arr.tofile(os.path.join(out_dir, f"{name}.bin"))
    shapes[name] = list(arr.shape)
    with open(os.path.join(out_dir, "shapes.json"), "w") as fh:
        json.dump(shapes, fh, indent=1)


def timed(name, fn, *args, **kwargs):
    t0 = time.perf_counter()
    res = fn(*args, **kwargs)
    timings[name] = time.perf_counter() - t0
    print(f"{name}: {timings[name]:.2f}s", flush=True)
    return res


scv.settings.verbosity = 1

# scvelo 0.2.5 + numpy >= 1.24: compute_cosines asks parallelize() to build a
# ragged numpy array, which now raises. Force list results instead.
_vgm = sys.modules["scvelo.tools.velocity_graph"]  # the package re-exports a function of the same name

_orig_parallelize = _vgm.parallelize


def _parallelize_list(*args, **kwargs):
    kwargs["as_array"] = False
    return _orig_parallelize(*args, **kwargs)


_vgm.parallelize = _parallelize_list
adata = scv.datasets.pancreas(os.path.join(data_dir, "endocrinogenesis_day15.h5ad"))
dump("raw_spliced", adata.layers["spliced"])
dump("raw_unspliced", adata.layers["unspliced"])
with open(os.path.join(out_dir, "raw_genes.txt"), "w") as fh:
    fh.write("\n".join(adata.var_names))
with open(os.path.join(out_dir, "clusters.txt"), "w") as fh:
    fh.write("\n".join(map(str, adata.obs["clusters"])))

timed("filter_and_normalize", scv.pp.filter_and_normalize, adata, min_shared_counts=20, n_top_genes=2000)
with open(os.path.join(out_dir, "genes.txt"), "w") as fh:
    fh.write("\n".join(adata.var_names))
dump("spliced", adata.layers["spliced"])
dump("unspliced", adata.layers["unspliced"])
timed("moments", scv.pp.moments, adata, n_pcs=30, n_neighbors=30)
dump("X_pca", adata.obsm["X_pca"])
dump("X_umap", adata.obsm["X_umap"])
D = adata.obsp["distances"].tocsr()
k = np.diff(D.indptr)
assert (k == k[0]).all(), "non-uniform neighbour counts"
order = np.argsort(D.data.reshape(-1, k[0]), axis=1)
dump("knn_idx", np.take_along_axis(D.indices.reshape(-1, k[0]), order, axis=1))
dump("knn_dist", np.take_along_axis(D.data.reshape(-1, k[0]), order, axis=1))
dump("Ms", adata.layers["Ms"])
dump("Mu", adata.layers["Mu"])

timed("velocity_stochastic", scv.tl.velocity, adata, mode="stochastic")
dump("velocity", adata.layers["velocity"])
dump("velocity_gamma", adata.var["velocity_gamma"].values)
dump("velocity_r2", adata.var["velocity_r2"].values)
dump("velocity_genes", adata.var["velocity_genes"].values.astype(float))
timed("velocity_graph", scv.tl.velocity_graph, adata, n_jobs=8)
G = adata.uns["velocity_graph"].tocoo()
Gn = adata.uns["velocity_graph_neg"].tocoo()
dump("graph_pos", np.vstack([G.row, G.col, G.data]).T)
dump("graph_neg", np.vstack([Gn.row, Gn.col, Gn.data]).T)
dump("self_transition", adata.obs["velocity_self_transition"].values)
timed("velocity_embedding", scv.tl.velocity_embedding, adata, basis="umap")
dump("velocity_umap", adata.obsm["velocity_umap"])

from scvelo.plotting.velocity_embedding_grid import compute_velocity_on_grid

Xg, Vg = compute_velocity_on_grid(adata.obsm["X_umap"], adata.obsm["velocity_umap"], autoscale=False, adjust_for_stream=True)
dump("stream_xs", Xg[0])
dump("stream_ys", Xg[1])
dump("stream_u", Vg[0])
dump("stream_v", Vg[1])
Xq, Vq = compute_velocity_on_grid(adata.obsm["X_umap"], adata.obsm["velocity_umap"])
dump("grid_X", Xq)
dump("grid_V", Vq)

timed("velocity_confidence", scv.tl.velocity_confidence, adata)
dump("velocity_length", adata.obs["velocity_length"].values)
dump("velocity_confidence", adata.obs["velocity_confidence"].values)
dump("velocity_confidence_transition", adata.obs["velocity_confidence_transition"].values)
timed("terminal_states", scv.tl.terminal_states, adata)
dump("root_cells", adata.obs["root_cells"].values)
dump("end_points", adata.obs["end_points"].values)
timed("velocity_pseudotime", scv.tl.velocity_pseudotime, adata)
dump("velocity_pseudotime", adata.obs["velocity_pseudotime"].values)

if run_dynamical:
    # align_dynamics(dm=...) crashes in scvelo 0.2.5 with recent numpy; aligning
    # afterwards without the model object gives the same parameters
    timed("recover_dynamics", scv.tl.recover_dynamics, adata, n_jobs=8, t_max=False)
    from scvelo.tools.dynamical_model import align_dynamics

    align_dynamics(adata)
    fit_keys = ["alpha", "beta", "gamma", "t_", "scaling", "std_u", "std_s", "likelihood", "pval_steady", "steady_u", "steady_s", "variance", "alignment_scaling", "r2"]
    for key in fit_keys:
        if f"fit_{key}" in adata.var:
            dump(f"fit_{key}", adata.var[f"fit_{key}"].values.astype(float))
    dump("fit_t", adata.layers["fit_t"])
    dump("fit_tau", adata.layers["fit_tau"])
    dump("fit_tau_", adata.layers["fit_tau_"])
    timed("velocity_dynamical", scv.tl.velocity, adata, mode="dynamical")
    dump("velocity_dyn", adata.layers["velocity"])
    dump("velocity_dyn_u", adata.layers["velocity_u"])
    dump("velocity_dyn_genes", adata.var["velocity_genes"].values.astype(float))
    timed("velocity_graph_dyn", scv.tl.velocity_graph, adata, n_jobs=8)
    timed("latent_time", scv.tl.latent_time, adata)
    dump("latent_time", adata.obs["latent_time"].values)
    dump("root_cells_dyn", adata.obs["root_cells"].values)
    dump("end_points_dyn", adata.obs["end_points"].values)

with open(os.path.join(out_dir, "shapes.json"), "w") as fh:
    json.dump(shapes, fh, indent=1)
with open(os.path.join(out_dir, "timings.json"), "w") as fh:
    json.dump(timings, fh, indent=1)
print("done")
