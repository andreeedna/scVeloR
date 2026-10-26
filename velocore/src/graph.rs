//! Velocity graph (cosine correlations between velocities and cell state
//! displacements) and transition matrices (scvelo.tl.velocity_graph,
//! scvelo.tl.transition_matrix).

use crate::matrix::{Csr, RowMat};
use crate::neighbors::Knn;
use crate::stats::percentile;
use rayon::prelude::*;

#[derive(Clone, Debug)]
pub struct VelocityGraphOptions {
    /// Recursion depth of the neighbour search (scVelo default: 2).
    pub n_recurse_neighbors: usize,
    /// Signed square-root transform of velocities and displacements
    /// (scVelo applies it by default for the stochastic model).
    pub sqrt_transform: bool,
    /// Randomly subsample neighbourhoods larger than this (None: keep all).
    pub max_neighs: Option<usize>,
    pub seed: u64,
}

impl Default for VelocityGraphOptions {
    fn default() -> Self {
        VelocityGraphOptions { n_recurse_neighbors: 2, sqrt_transform: false, max_neighs: None, seed: 0 }
    }
}

#[derive(Clone, Debug)]
pub struct VelocityGraph {
    /// Signed cosine correlations (cells x cells); positive and negative parts
    /// of scVelo's `velocity_graph` / `velocity_graph_neg`.
    pub graph: Csr,
    /// Self-transition probability per cell.
    pub self_transition: Vec<f64>,
}

impl VelocityGraph {
    pub fn positive(&self) -> Csr {
        let mut g = self.graph.map(|v| v.clamp(0.0, 1.0));
        g.eliminate_zeros();
        g
    }
    pub fn negative(&self) -> Csr {
        let mut g = self.graph.map(|v| v.clamp(-1.0, 0.0));
        g.eliminate_zeros();
        g
    }
}

#[inline]
fn signed_sqrt(v: f64) -> f64 {
    v.abs().sqrt() * v.signum_or_zero()
}

trait SignumOrZero {
    fn signum_or_zero(self) -> f64;
}
impl SignumOrZero for f64 {
    #[inline]
    fn signum_or_zero(self) -> f64 {
        if self > 0.0 {
            1.0
        } else if self < 0.0 {
            -1.0
        } else {
            0.0
        }
    }
}

/// Neighbours up to the given recursion depth, including the cell itself,
/// sorted and unique (scvelo get_iterative_indices).
pub fn iterative_neighbors(knn: &Knn, i: usize, depth: usize) -> Vec<usize> {
    let mut cur: Vec<usize> = vec![i];
    for _ in 0..depth {
        let mut next = cur.clone();
        for &j in &cur {
            next.extend_from_slice(knn.neighbors(j));
        }
        next.sort_unstable();
        next.dedup();
        cur = next;
    }
    cur
}

/// Small deterministic xorshift RNG for neighbour subsampling.
struct XorShift(u64);
impl XorShift {
    fn next(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }
}

/// Compute the velocity graph.
/// `x`: cells x genes expression (Ms) restricted to velocity genes;
/// `v`: cells x genes velocities on the same genes.
/// Genes with any NaN velocity are dropped beforehand by the caller or here.
pub fn velocity_graph(x: &RowMat, v: &RowMat, knn: &Knn, opts: &VelocityGraphOptions) -> VelocityGraph {
    let (n, g0) = (x.nrow, x.ncol);
    assert_eq!((v.nrow, v.ncol), (n, g0));
    // drop genes with NaN velocities
    let keep: Vec<usize> = (0..g0).filter(|&j| (0..n).all(|i| !v.get(i, j).is_nan())).collect();
    let (x, v) = if keep.len() < g0 { (x.select_cols(&keep), v.select_cols(&keep)) } else { (x.clone(), v.clone()) };
    let g = x.ncol;

    // transform + center velocities per cell
    let mut vt = v.clone();
    vt.data.par_chunks_mut(g.max(1)).for_each(|row| {
        if opts.sqrt_transform {
            for e in row.iter_mut() {
                *e = signed_sqrt(*e);
            }
        }
        let m = row.iter().sum::<f64>() / g as f64;
        for e in row.iter_mut() {
            *e -= m;
        }
    });

    let rows: Vec<Vec<(usize, f64)>> = (0..n)
        .into_par_iter()
        .map(|i| {
            let vi = vt.row(i);
            let vmax = vi.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            let vmin = vi.iter().copied().fold(f64::INFINITY, f64::min);
            if !(vmax != 0.0 || vmin != 0.0) {
                return Vec::new();
            }
            let vnorm = vi.iter().map(|a| a * a).sum::<f64>().sqrt();
            if vnorm == 0.0 || !vnorm.is_finite() {
                return Vec::new();
            }
            let mut neighs = iterative_neighbors(knn, i, opts.n_recurse_neighbors);
            if let Some(m) = opts.max_neighs {
                if neighs.len() > m {
                    let mut rng = XorShift(opts.seed ^ (i as u64).wrapping_mul(0x9E3779B97F4A7C15) | 1);
                    // partial Fisher-Yates
                    for k in 0..m {
                        let r = k + (rng.next() as usize) % (neighs.len() - k);
                        neighs.swap(k, r);
                    }
                    neighs.truncate(m);
                    neighs.sort_unstable();
                }
            }
            let xi = x.row(i);
            let mut out = Vec::with_capacity(neighs.len());
            let mut dx = vec![0.0; g];
            for &j in &neighs {
                let xj = x.row(j);
                let mut s = 0.0;
                for k in 0..g {
                    let mut d = xj[k] - xi[k];
                    if opts.sqrt_transform {
                        d = signed_sqrt(d);
                    }
                    dx[k] = d;
                    s += d;
                }
                let m = s / g as f64;
                let mut dot = 0.0;
                let mut nn = 0.0;
                for k in 0..g {
                    let d = dx[k] - m;
                    dot += d * vi[k];
                    nn += d * d;
                }
                let c = dot / (nn.sqrt() * vnorm);
                if c.is_finite() && c != 0.0 {
                    out.push((j, c));
                }
            }
            out
        })
        .collect();
    let graph = Csr::from_rows(n, rows);
    let pos_max: Vec<f64> = (0..n)
        .map(|i| {
            let (_, val) = graph.row(i);
            val.iter().copied().fold(0.0f64, |a, b| a.max(b.min(1.0)))
        })
        .collect();
    let ub = percentile(&pos_max, 98.0);
    let self_transition = pos_max.iter().map(|c| (ub - c).clamp(0.0, 1.0)).collect();
    VelocityGraph { graph, self_transition }
}

#[derive(Clone, Debug)]
pub struct TransitionOptions {
    pub scale: f64,
    pub self_transitions: bool,
    pub use_negative_cosines: bool,
    pub backward: bool,
    /// Drop transitions below this percentile of all transition values.
    pub perc: Option<f64>,
    pub threshold: Option<f64>,
}

impl Default for TransitionOptions {
    fn default() -> Self {
        TransitionOptions {
            scale: 10.0,
            self_transitions: true,
            use_negative_cosines: false,
            backward: false,
            perc: None,
            threshold: None,
        }
    }
}

/// Transition matrix T_ij ∝ exp(scale · cos_ij), row normalised
/// (scvelo.tl.transition_matrix without an embedding-based diffusion kernel).
pub fn transition_matrix(vg: &VelocityGraph, opts: &TransitionOptions) -> Csr {
    let n = vg.graph.nrow;
    let scale = opts.scale;
    let rows: Vec<Vec<(usize, f64)>> = (0..n)
        .into_par_iter()
        .map(|i| {
            let (idx, val) = vg.graph.row(i);
            let mut out: Vec<(usize, f64)> = Vec::with_capacity(idx.len() + 1);
            let mut has_diag = false;
            for (&j, &c) in idx.iter().zip(val) {
                let c = if j == i && opts.self_transitions { vg.self_transition[i] } else { c };
                if j == i {
                    has_diag = true;
                }
                let t = if c > 0.0 {
                    let c = c.min(1.0);
                    if opts.use_negative_cosines { (c * scale).exp_m1() } else { (c * scale).exp() }
                } else if c < 0.0 {
                    let c = c.max(-1.0);
                    if opts.use_negative_cosines { -(-c * scale).exp_m1() } else { (c * scale).exp() }
                } else {
                    0.0
                };
                if t != 0.0 {
                    out.push((j, t));
                }
            }
            if opts.self_transitions && !has_diag {
                let c = vg.self_transition[i];
                if c > 0.0 {
                    let t = if opts.use_negative_cosines { (c * scale).exp_m1() } else { (c * scale).exp() };
                    out.push((i, t));
                }
            }
            out
        })
        .collect();
    let mut t = Csr::from_rows(n, rows);
    if opts.perc.is_some() || opts.threshold.is_some() {
        let thr = match opts.threshold {
            Some(t) => t,
            None => percentile(&t.data, opts.perc.unwrap()),
        };
        for v in t.data.iter_mut() {
            if *v < thr {
                *v = 0.0;
            }
        }
        t.eliminate_zeros();
    }
    if opts.backward {
        t = t.transpose();
    }
    t.normalize_rows_abs();
    t
}

/// Symmetric velocity-graph kernel |cos| + |cos|ᵀ used for the velocity
/// pseudotime diffusion map, followed by scanpy's density normalisation:
/// returns the symmetric matrix Z K Z ready for an eigen decomposition.
pub fn pseudotime_kernel(vg: &VelocityGraph) -> Csr {
    let a = vg.graph.map(|v| v.abs());
    let at = a.transpose();
    let n = a.nrow;
    // T = A + Aᵀ
    let rows: Vec<Vec<(usize, f64)>> = (0..n)
        .into_par_iter()
        .map(|i| {
            let mut m: std::collections::BTreeMap<usize, f64> = std::collections::BTreeMap::new();
            let (ai, av) = a.row(i);
            for (&j, &v) in ai.iter().zip(av) {
                *m.entry(j).or_insert(0.0) += v;
            }
            let (bi, bv) = at.row(i);
            for (&j, &v) in bi.iter().zip(bv) {
                *m.entry(j).or_insert(0.0) += v;
            }
            m.into_iter().filter(|e| e.1 != 0.0).collect()
        })
        .collect();
    let t = Csr::from_rows(n, rows);
    // density normalisation K = Q T Q with q = column sums
    let mut q = vec![0.0; n];
    for i in 0..n {
        let (idx, val) = t.row(i);
        for (&j, &v) in idx.iter().zip(val) {
            q[j] += v;
        }
    }
    for v in q.iter_mut() {
        if *v == 0.0 {
            *v = 1.0;
        }
    }
    let mut k = t.clone();
    for i in 0..n {
        let (a0, b0) = (k.indptr[i], k.indptr[i + 1]);
        for p in a0..b0 {
            let j = k.indices[p];
            k.data[p] /= q[i] * q[j];
        }
    }
    let mut z = vec![0.0; n];
    for i in 0..n {
        let (idx, val) = k.row(i);
        for (&j, &v) in idx.iter().zip(val) {
            z[j] += v;
        }
    }
    let z: Vec<f64> = z.iter().map(|v| v.sqrt()).collect();
    for i in 0..n {
        let (a0, b0) = (k.indptr[i], k.indptr[i + 1]);
        for p in a0..b0 {
            let j = k.indices[p];
            k.data[p] /= z[i] * z[j];
        }
    }
    k
}
