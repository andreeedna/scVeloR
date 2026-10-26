//! R bindings for velocore. All functions take flat vectors plus dimensions;
//! the R wrappers in R/core.R take care of layouts:
//! * a genes x cells R matrix is passed as-is (its column-major memory is the
//!   row-major cells x genes layout velocore expects);
//! * kNN indices are passed 0-based, row-major (n x k);
//! * sparse matrices are passed as dgCMatrix slots (p, i, x) of genes x cells.

#![allow(non_snake_case)]

use extendr_api::prelude::*;
use velocore::dynamics::{self, DynamicsFit, DynamicsOptions, GeneDynamics};
use velocore::embedding::{self, GridOptions};
use velocore::graph::{self, TransitionOptions, VelocityGraph, VelocityGraphOptions};
use velocore::neighbors::{self, Knn};
use velocore::stream::{self, Direction, StreamOptions};
use velocore::velocity::{self, SteadyStateOptions};
use velocore::{Csr, RowMat};

fn named(names: &[&str], values: Vec<Robj>) -> Robj {
    List::from_names_and_values(names.iter().copied(), values).unwrap().into_robj()
}

fn dense(x: &[f64], nrow: i32, ncol: i32) -> RowMat {
    RowMat::from_vec(nrow as usize, ncol as usize, x.to_vec())
}

fn knn_from(idx: &[i32], dist: &[f64], n: i32, k: i32) -> Knn {
    Knn::from_indices(n as usize, k as usize, idx.iter().map(|&v| v as usize).collect(), dist.to_vec())
}

/// dgCMatrix (genes x cells) slots -> CSR cells x genes.
fn csr_from_dgc(p: &[i32], i: &[i32], x: &[f64], ngenes: i32) -> Csr {
    Csr::new(
        p.len() - 1,
        ngenes as usize,
        p.iter().map(|&v| v as usize).collect(),
        i.iter().map(|&v| v as usize).collect(),
        x.to_vec(),
    )
}

/// Triplets (1-based i, j, x) -> CSR.
fn csr_from_triplets(i: &[i32], j: &[i32], x: &[f64], n: i32) -> Csr {
    let r: Vec<usize> = i.iter().map(|&v| v as usize - 1).collect();
    let c: Vec<usize> = j.iter().map(|&v| v as usize - 1).collect();
    Csr::from_triplets(n as usize, n as usize, &r, &c, x)
}

fn triplets_robj(m: &Csr) -> Robj {
    let (r, c, v) = m.to_triplets();
    named(
        &["i", "j", "x"],
        vec![
            r.iter().map(|&a| a as i32 + 1).collect::<Vec<i32>>().into_robj(),
            c.iter().map(|&a| a as i32 + 1).collect::<Vec<i32>>().into_robj(),
            v.into_robj(),
        ],
    )
}

fn set_threads(n_threads: i32) {
    if n_threads > 0 {
        let _ = rayon::ThreadPoolBuilder::new().num_threads(n_threads as usize).build_global();
    }
}

/// Exact kNN on the rows of a cells x dims matrix (passed as t(X)).
/// @noRd
#[extendr]
fn rs_knn(x: &[f64], n: i32, d: i32, k: i32, n_threads: i32) -> Robj {
    set_threads(n_threads);
    let m = RowMat::from_vec(n as usize, d as usize, x.to_vec());
    let res = neighbors::knn(&m, k as usize);
    named(
        &["idx", "dist", "k"],
        vec![
            res.indices.iter().map(|&v| v as i32).collect::<Vec<i32>>().into_robj(),
            res.distances.into_robj(),
            (res.k as i32).into_robj(),
        ],
    )
}

/// Connectivity matrices derived from the kNN graph.
/// mode: "moments" (union graph, scVelo moments) or "distances" (directed kNN).
/// @noRd
#[extendr]
fn rs_connectivities(idx: &[i32], dist: &[f64], n: i32, k: i32, mode: &str, recurse: bool) -> Robj {
    let knn = knn_from(idx, dist, n, k);
    let c = match mode {
        "moments" => neighbors::moment_connectivities(&knn),
        _ => neighbors::knn_connectivities(&knn, recurse),
    };
    triplets_robj(&c)
}

/// First (and optionally second) order moments of spliced/unspliced counts.
/// @noRd
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_moments(
    idx: &[i32],
    dist: &[f64],
    n: i32,
    k: i32,
    s_p: &[i32],
    s_i: &[i32],
    s_x: &[f64],
    u_p: &[i32],
    u_i: &[i32],
    u_x: &[f64],
    ngenes: i32,
    second_order: bool,
) -> Robj {
    let knn = knn_from(idx, dist, n, k);
    let conn = neighbors::moment_connectivities(&knn);
    let s = csr_from_dgc(s_p, s_i, s_x, ngenes);
    let u = csr_from_dgc(u_p, u_i, u_x, ngenes);
    let ms = conn.dot_sparse_dense_out(&s);
    let mu = conn.dot_sparse_dense_out(&u);
    if second_order {
        let (mss, mus) = velocore::moments::second_order(&conn, &s, &u);
        named(&["Ms", "Mu", "Mss", "Mus"], vec![ms.data.into_robj(), mu.data.into_robj(), mss.data.into_robj(), mus.data.into_robj()])
    } else {
        named(&["Ms", "Mu"], vec![ms.data.into_robj(), mu.data.into_robj()])
    }
}

/// Steady-state / stochastic velocity. ms, mu, mss, mus: genes x cells.
/// @noRd
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_velocity_steady(
    ms: &[f64],
    mu: &[f64],
    mss: &[f64],
    mus: &[f64],
    n: i32,
    g: i32,
    stochastic: bool,
    fit_offset: bool,
    fit_offset2: bool,
    perc_lo: f64,
    perc_hi: f64,
    min_r2: f64,
    min_ratio: f64,
    r2_adjusted: bool,
    highly_variable: &[i32],
    n_threads: i32,
) -> Robj {
    set_threads(n_threads);
    let ms = dense(ms, n, g);
    let mu = dense(mu, n, g);
    let opts = SteadyStateOptions {
        fit_offset,
        fit_offset2,
        perc: if perc_hi.is_nan() { None } else { Some((perc_lo, perc_hi)) },
        min_r2,
        min_ratio,
        r2_adjusted,
        constrain_ratio: None,
    };
    let hv: Option<Vec<bool>> = if highly_variable.len() == g as usize { Some(highly_variable.iter().map(|&v| v != 0).collect()) } else { None };
    let res = if stochastic {
        let mss = dense(mss, n, g);
        let mus = dense(mus, n, g);
        velocity::stochastic(&ms, &mu, &mss, &mus, hv.as_deref(), &opts)
    } else {
        velocity::deterministic(&ms, &mu, hv.as_deref(), &opts)
    };
    let var_vel = res.variance_velocity.map(|m| m.data).unwrap_or_default();
    named(
        &["velocity", "variance_velocity", "gamma", "offset", "offset2", "r2", "qreg_ratio", "velocity_genes", "low_signal_r2"],
        vec![
            res.velocity.data.into_robj(),
            var_vel.into_robj(),
            res.gamma.into_robj(),
            res.offset.into_robj(),
            res.offset2.into_robj(),
            res.r2.into_robj(),
            res.qreg_ratio.into_robj(),
            res.velocity_genes.iter().map(|&b| b as i32).collect::<Vec<i32>>().into_robj(),
            res.low_signal_r2_threshold.unwrap_or(f64::NAN).into_robj(),
        ],
    )
}

/// Velocity graph. x, v: genes x cells (velocity genes only).
/// @noRd
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_velocity_graph(
    x: &[f64],
    v: &[f64],
    n: i32,
    g: i32,
    idx: &[i32],
    dist: &[f64],
    k: i32,
    n_recurse: i32,
    sqrt_transform: bool,
    max_neighs: i32,
    seed: i32,
    n_threads: i32,
) -> Robj {
    set_threads(n_threads);
    let knn = knn_from(idx, dist, n, k);
    let opts = VelocityGraphOptions {
        n_recurse_neighbors: n_recurse as usize,
        sqrt_transform,
        max_neighs: if max_neighs > 0 { Some(max_neighs as usize) } else { None },
        seed: seed as u64,
    };
    let vg = graph::velocity_graph(&dense(x, n, g), &dense(v, n, g), &knn, &opts);
    named(&["graph", "self_transition"], vec![triplets_robj(&vg.graph), vg.self_transition.into_robj()])
}

fn vgraph_from(i: &[i32], j: &[i32], x: &[f64], n: i32, self_transition: &[f64]) -> VelocityGraph {
    VelocityGraph { graph: csr_from_triplets(i, j, x, n), self_transition: self_transition.to_vec() }
}

/// Transition matrix from the velocity graph (triplets).
/// @noRd
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_transition_matrix(
    i: &[i32],
    j: &[i32],
    x: &[f64],
    n: i32,
    self_transition: &[f64],
    scale: f64,
    self_transitions: bool,
    use_negative_cosines: bool,
    backward: bool,
) -> Robj {
    let vg = vgraph_from(i, j, x, n, self_transition);
    let t = graph::transition_matrix(
        &vg,
        &TransitionOptions { scale, self_transitions, use_negative_cosines, backward, perc: None, threshold: None },
    );
    triplets_robj(&t)
}

/// Project velocities onto an embedding. x_emb passed as t(X) (dims x cells).
/// @noRd
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_velocity_embedding(
    i: &[i32],
    j: &[i32],
    x: &[f64],
    n: i32,
    self_transition: &[f64],
    scale: f64,
    x_emb: &[f64],
    d: i32,
    autoscale: bool,
) -> Robj {
    let vg = vgraph_from(i, j, x, n, self_transition);
    let t = graph::transition_matrix(
        &vg,
        &TransitionOptions { scale, self_transitions: true, use_negative_cosines: true, backward: false, perc: None, threshold: None },
    );
    let xe = RowMat::from_vec(n as usize, d as usize, x_emb.to_vec());
    let mut ve = embedding::project(&t, &xe);
    if autoscale && d >= 2 {
        let f = 3.0 * embedding::quiver_autoscale(&xe, &ve);
        if f > 0.0 && f.is_finite() {
            ve.data.iter_mut().for_each(|v| *v /= f);
        }
    }
    ve.data.into_robj()
}

/// Gaussian-smoothed velocities on a grid (arrows or stream input).
/// @noRd
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_velocity_grid(
    x_emb: &[f64],
    v_emb: &[f64],
    n: i32,
    density: f64,
    smooth: f64,
    n_neighbors: i32,
    min_mass: f64,
    autoscale: bool,
    stream_mode: bool,
    cutoff_perc: f64,
) -> Robj {
    let xe = RowMat::from_vec(n as usize, 2, x_emb.to_vec());
    let ve = RowMat::from_vec(n as usize, 2, v_emb.to_vec());
    let gv = embedding::velocity_on_grid(
        &xe,
        &ve,
        &GridOptions {
            density,
            smooth,
            n_neighbors: if n_neighbors > 0 { Some(n_neighbors as usize) } else { None },
            min_mass: if min_mass.is_nan() { None } else { Some(min_mass) },
            autoscale,
            adjust_for_stream: stream_mode,
            cutoff_perc,
        },
    );
    named(&["xs", "ys", "points", "velocity"], vec![gv.xs.into_robj(), gv.ys.into_robj(), gv.points.data.into_robj(), gv.velocity.data.into_robj()])
}

/// Streamlines through a gridded velocity field (u, v given as [y, x], x fastest).
/// @noRd
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_streamlines(xs: &[f64], ys: &[f64], u: &[f64], v: &[f64], density_x: f64, density_y: f64, min_length: f64, max_length: f64, direction: &str) -> Robj {
    let dir = match direction {
        "forward" => Direction::Forward,
        "backward" => Direction::Backward,
        _ => Direction::Both,
    };
    let lines = stream::streamlines(
        xs,
        ys,
        u,
        v,
        &StreamOptions { density: (density_x, density_y), min_length, max_length, direction: dir, broken_streamlines: true },
    );
    let mut id = Vec::new();
    let mut lx = Vec::new();
    let mut ly = Vec::new();
    let mut sp = Vec::new();
    let mut a = [Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new()];
    for (k, l) in lines.iter().enumerate() {
        for p in 0..l.x.len() {
            id.push(k as i32 + 1);
            lx.push(l.x[p]);
            ly.push(l.y[p]);
            sp.push(l.speed[p]);
        }
        a[0].push((k + 1) as f64);
        a[1].push(l.arrow_tail.0);
        a[2].push(l.arrow_tail.1);
        a[3].push(l.arrow_head.0);
        a[4].push(l.arrow_head.1);
    }
    let [a0, a1, a2, a3, a4] = a;
    named(
        &["line", "x", "y", "speed", "arrow_line", "arrow_x0", "arrow_y0", "arrow_x1", "arrow_y1"],
        vec![id.into_robj(), lx.into_robj(), ly.into_robj(), sp.into_robj(), a0.into_robj(), a1.into_robj(), a2.into_robj(), a3.into_robj(), a4.into_robj()],
    )
}

/// Velocity length / confidence. v: genes x cells.
/// @noRd
#[extendr]
fn rs_velocity_confidence(v: &[f64], n: i32, g: i32, idx: &[i32], dist: &[f64], k: i32) -> Robj {
    let knn = knn_from(idx, dist, n, k);
    let (len, conf) = velocore::confidence::velocity_confidence(&dense(v, n, g), &knn);
    named(&["length", "confidence"], vec![len.into_robj(), conf.into_robj()])
}

/// Transition confidence. t: triplets of the transition matrix.
/// @noRd
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_confidence_transition(i: &[i32], j: &[i32], tx: &[f64], n: i32, x: &[f64], v: &[f64], g: i32) -> Vec<f64> {
    let t = csr_from_triplets(i, j, tx, n);
    velocore::confidence::confidence_transition(&t, &dense(x, n, g), &dense(v, n, g))
}

/// Dynamical model fit. mu, ms: genes x cells for the genes to fit.
/// @noRd
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_recover_dynamics(
    mu: &[f64],
    ms: &[f64],
    n: i32,
    g: i32,
    idx: &[i32],
    dist: &[f64],
    k: i32,
    max_iter: i32,
    fit_scaling: bool,
    fit_steady_states: bool,
    fit_connected_states: bool,
    projection: bool,
    t_max: f64,
    n_threads: i32,
) -> Robj {
    set_threads(n_threads);
    let knn = knn_from(idx, dist, n, k);
    let conn = neighbors::moment_connectivities(&knn);
    let opts = DynamicsOptions { max_iter: max_iter as usize, fit_scaling, fit_steady_states, fit_connected_states, projection, perc: 99.0 };
    let fit = dynamics::recover_dynamics(&dense(mu, n, g), &dense(ms, n, g), Some(&conn), &opts, if t_max.is_nan() { None } else { Some(t_max) });
    let gs = &fit.genes;
    let col = |f: &dyn Fn(&GeneDynamics) -> f64| -> Robj { gs.iter().map(f).collect::<Vec<f64>>().into_robj() };
    named(
        &[
            "alpha", "beta", "gamma", "t_", "scaling", "std_u", "std_s", "likelihood", "variance", "pval_steady", "steady_u", "steady_s",
            "alignment_scaling", "n_loss", "fit_t", "fit_tau", "fit_tau_",
        ],
        vec![
            col(&|x| x.alpha),
            col(&|x| x.beta),
            col(&|x| x.gamma),
            col(&|x| x.t_),
            col(&|x| x.scaling),
            col(&|x| x.std_u),
            col(&|x| x.std_s),
            col(&|x| x.likelihood),
            col(&|x| x.variance),
            col(&|x| x.pval_steady),
            col(&|x| x.steady_u),
            col(&|x| x.steady_s),
            fit.alignment_scaling.clone().into_robj(),
            col(&|x| x.loss.len() as f64),
            // cells x genes row-major == genes x cells in R
            fit.fit_t.data.into_robj(),
            fit.fit_tau.data.into_robj(),
            fit.fit_tau_.data.into_robj(),
        ],
    )
}

/// Velocities from fitted dynamics. All matrices genes x cells.
/// @noRd
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_dynamical_velocity(
    mu: &[f64],
    ms: &[f64],
    n: i32,
    g: i32,
    alpha: &[f64],
    beta: &[f64],
    gamma: &[f64],
    t_: &[f64],
    scaling: &[f64],
    std_u: &[f64],
    std_s: &[f64],
    fit_tau: &[f64],
    fit_tau_: &[f64],
    idx: &[i32],
    dist: &[f64],
    k: i32,
    use_connectivities: bool,
) -> Robj {
    let (nu, gu) = (n as usize, g as usize);
    let genes: Vec<GeneDynamics> = (0..gu)
        .map(|j| GeneDynamics {
            recoverable: !alpha[j].is_nan(),
            alpha: alpha[j],
            beta: beta[j],
            gamma: gamma[j],
            t_: t_[j],
            scaling: scaling[j],
            std_u: std_u[j],
            std_s: std_s[j],
            likelihood: f64::NAN,
            variance: f64::NAN,
            pval_steady: f64::NAN,
            steady_u: f64::NAN,
            steady_s: f64::NAN,
            t: vec![],
            tau: vec![],
            tau_: vec![],
            loss: vec![],
        })
        .collect();
    let fit = DynamicsFit {
        genes,
        fit_t: RowMat::zeros(nu, gu),
        fit_tau: dense(fit_tau, n, g),
        fit_tau_: dense(fit_tau_, n, g),
        alignment_scaling: vec![f64::NAN; gu],
    };
    let knn = knn_from(idx, dist, n, k);
    let conn = neighbors::moment_connectivities(&knn);
    let (vs, vu) = dynamics::dynamical_velocity(&dense(mu, n, g), &dense(ms, n, g), &fit, if use_connectivities { Some(&conn) } else { None });
    named(&["velocity", "velocity_u"], vec![vs.data.into_robj(), vu.data.into_robj()])
}

/// Latent time. fit_t: genes x cells (selected genes); roots/fates 1-based.
/// @noRd
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_latent_time(fit_t: &[f64], n: i32, g: i32, roots: &[i32], fates: &[i32], idx: &[i32], dist: &[f64], k: i32, min_confidence: f64) -> Vec<f64> {
    let knn = knn_from(idx, dist, n, k);
    let conn = neighbors::moment_connectivities(&knn);
    let r: Vec<usize> = roots.iter().map(|&v| v as usize - 1).collect();
    let f: Vec<usize> = fates.iter().map(|&v| v as usize - 1).collect();
    velocore::latent::latent_time(&dense(fit_t, n, g), &r, &f, &conn, min_confidence)
}

/// Symmetric velocity kernel for the velocity pseudotime diffusion map.
/// @noRd
#[extendr]
fn rs_pseudotime_kernel(i: &[i32], j: &[i32], x: &[f64], n: i32) -> Robj {
    let vg = VelocityGraph { graph: csr_from_triplets(i, j, x, n), self_transition: vec![0.0; n as usize] };
    triplets_robj(&graph::pseudotime_kernel(&vg))
}

/// DPT distance from a root using eigenpairs (evecs: k x cells, i.e. t(V)).
/// @noRd
#[extendr]
fn rs_dpt(evals: &[f64], evecs: &[f64], n: i32, root: i32) -> Vec<f64> {
    let k = evals.len();
    let ev = RowMat::from_vec(n as usize, k, evecs.to_vec());
    velocore::latent::dpt_from_root(evals, &ev, root as usize - 1)
}

/// Model trajectory of one gene for phase portraits (unscaled units).
/// @noRd
#[extendr]
fn rs_gene_trajectory(alpha: f64, beta: f64, gamma: f64, t_: f64, scaling: f64, tmax: f64, num: i32) -> Robj {
    let gd = GeneDynamics {
        recoverable: true,
        alpha,
        beta,
        gamma,
        t_,
        scaling,
        std_u: 1.0,
        std_s: 1.0,
        likelihood: f64::NAN,
        variance: f64::NAN,
        pval_steady: f64::NAN,
        steady_u: f64::NAN,
        steady_s: f64::NAN,
        t: vec![tmax],
        tau: vec![],
        tau_: vec![],
        loss: vec![],
    };
    let (t, u, s) = dynamics::gene_trajectory(&gd, num as usize);
    named(&["t", "u", "s"], vec![t.into_robj(), u.into_robj(), s.into_robj()])
}

extendr_module! {
    mod scVelo2;
    fn rs_knn;
    fn rs_connectivities;
    fn rs_moments;
    fn rs_velocity_steady;
    fn rs_velocity_graph;
    fn rs_transition_matrix;
    fn rs_velocity_embedding;
    fn rs_velocity_grid;
    fn rs_streamlines;
    fn rs_velocity_confidence;
    fn rs_confidence_transition;
    fn rs_recover_dynamics;
    fn rs_dynamical_velocity;
    fn rs_latent_time;
    fn rs_pseudotime_kernel;
    fn rs_dpt;
    fn rs_gene_trajectory;
}
