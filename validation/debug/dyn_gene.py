"""Trace scVelo's DynamicsRecovery for one gene (debugging aid)."""
import sys, os
import numpy as np
import scvelo as scv
from scvelo.tools.dynamical_model import DynamicsRecovery
from scvelo.preprocessing.moments import get_connectivities

data_dir, gene_idx = sys.argv[1], int(sys.argv[2])
scv.settings.verbosity = 0
adata = scv.datasets.pancreas(os.path.join(data_dir, "endocrinogenesis_day15.h5ad"))
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
gene = adata.var_names[gene_idx]
conn = get_connectivities(adata)
dm = DynamicsRecovery(adata, gene, use_raw=False, max_iter=10, fit_connected_states=conn, fit_scaling=True, fit_steady_states=True)
print("init pars", dm.alpha, dm.beta, dm.gamma, dm.t_, dm.scaling, "pval", dm.pval_steady)
print("init loss", dm.loss)
stages = [("t_and_alpha", lambda: dm.fit_t_and_alpha()), ("scaling", lambda: dm.fit_scaling_()),
          ("rates", lambda: dm.fit_rates()), ("t_", lambda: dm.fit_t_()), ("t_and_rates", lambda: dm.fit_t_and_rates())]
for name, f in stages:
    f()
    print(name, dm.alpha, dm.beta, dm.gamma, dm.t_, dm.scaling, "loss", dm.loss[-1], "nloss", len(dm.loss))
dm.assignment_mode = "projection"
dm.update(adjust_t_=False)
print("proj update", dm.alpha, dm.beta, dm.gamma, dm.t_, dm.scaling, "loss", dm.loss[-1], "nloss", len(dm.loss))
dm.fit_t_and_rates(refit_time=False)
print("t_and_rates fixed", dm.alpha, dm.beta, dm.gamma, dm.t_, dm.scaling, "loss", dm.loss[-1], "nloss", len(dm.loss))
dm.update()
print("final", dm.alpha, dm.beta, dm.gamma, dm.t_, dm.scaling, "loss", dm.loss[-1], "nloss", len(dm.loss))
