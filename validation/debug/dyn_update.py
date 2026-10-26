"""Instrumented scVelo update() for one gene."""
import sys, os
import numpy as np
import scvelo as scv
from scvelo.tools import dynamical_model as dmod
from scvelo.tools.dynamical_model_utils import unspliced
from scvelo.preprocessing.moments import get_connectivities

def update(self, t=None, t_=None, alpha=None, beta=None, gamma=None, scaling=None, u0_=None, s0_=None, adjust_t_=True):
    loss_prev = self.loss[-1] if len(self.loss) > 0 else 1e6
    alpha, beta, gamma, scaling, t_ = self.get_vars(alpha, beta, gamma, scaling, t_, u0_, s0_)
    t, tau, o = self.get_time_assignment(alpha, beta, gamma, scaling, t_, u0_, s0_, t)
    loss = self.get_loss(t, t_, alpha, beta, gamma, scaling)
    perform_update = loss < loss_prev
    on = self.o == 1
    if adjust_t_ and np.any(on):
        if not perform_update:
            alpha, beta, gamma, scaling, t_ = self.get_vars()
            t, tau, o = self.get_time_assignment()
            loss = self.get_loss()
        alt_t_ = t[on].max()
        if 0 < alt_t_ < t_:
            n_eq = np.sum(t == t_)
            alt_t_ += np.max(t) / len(t) * np.sum(t == t_)
            alt_t, alt_tau, alt_o = self.get_time_assignment(alpha, beta, gamma, scaling, alt_t_)
            alt_loss = self.get_loss(alt_t, alt_t_, alpha, beta, gamma, scaling)
            ut_cur = unspliced(t_, 0, alpha, beta)
            ut_alt = unspliced(alt_t_, 0, alpha, beta)
            min_loss = np.min([loss, loss_prev])
            print(f"  alt: t_={t_} alt_t_={alt_t_} tmax={np.max(t)} n_eq={n_eq} loss={loss} prev={loss_prev} alt_loss={alt_loss} ut_cur={ut_cur} ut_alt={ut_alt} dtype={t.dtype}")
            if alt_loss * 0.99 <= min_loss or ut_cur * 0.99 < ut_alt:
                t, tau, o, t_, loss = alt_t, alt_tau, alt_o, alt_t_, alt_loss
                perform_update = True
    if perform_update:
        if scaling is not None:
            self.steady_u *= self.scaling / scaling
            self.u0_ *= self.scaling / scaling
        if u0_ is not None:
            self.u0_ = u0_
        if s0_ is not None:
            self.s0_ = s0_
        self.t, self.tau, self.o = t, tau, o
        self.alpha, self.beta, self.gamma = alpha, beta, gamma
        self.scaling, self.t_ = scaling, t_
        new_pars = np.array([alpha, beta, gamma, t_, scaling])[:, None]
        self.pars = np.c_[self.pars, new_pars]
        self.loss.append(loss)
    return perform_update

dmod.DynamicsRecovery.update = update
data_dir, gene_idx = sys.argv[1], int(sys.argv[2])
scv.settings.verbosity = 0
adata = scv.datasets.pancreas(os.path.join(data_dir, "endocrinogenesis_day15.h5ad"))
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
gene = adata.var_names[gene_idx]
conn = get_connectivities(adata)
dm = dmod.DynamicsRecovery(adata, gene, use_raw=False, max_iter=10, fit_connected_states=conn)
print("after init", dm.alpha, dm.beta, dm.gamma, dm.t_, dm.scaling, dm.loss)
dm.fit_t_and_alpha()
print("t_and_alpha", dm.alpha, dm.beta, dm.gamma, dm.t_, dm.scaling, dm.loss[-1])
