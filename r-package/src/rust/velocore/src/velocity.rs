//! Steady-state ("deterministic") and second-order ("stochastic") velocity
//! models, ported from scvelo.tools.velocity / scvelo.core.LinearRegression /
//! scvelo.tools.optimization.leastsq_generalized.

use crate::matrix::RowMat;
use crate::stats::{mean, percentile, std};
use rayon::prelude::*;

#[derive(Clone, Debug)]
pub struct SteadyStateOptions {
    /// Fit an intercept (offset) for the first order moments.
    pub fit_offset: bool,
    /// Fit an intercept for the second order moments (stochastic only).
    pub fit_offset2: bool,
    /// Extreme-quantile percentiles (lower, upper). Only the upper is used
    /// without an offset. `None` fits on all cells.
    pub perc: Option<(f64, f64)>,
    pub min_r2: f64,
    pub min_ratio: f64,
    /// Compute R² on the full-data OLS fit instead of the quantile fit.
    pub r2_adjusted: bool,
    /// Optional clipping bounds for gamma.
    pub constrain_ratio: Option<(f64, f64)>,
}

impl Default for SteadyStateOptions {
    fn default() -> Self {
        SteadyStateOptions {
            fit_offset: false,
            fit_offset2: false,
            perc: Some((5.0, 95.0)),
            min_r2: 0.01,
            min_ratio: 0.01,
            r2_adjusted: false,
            constrain_ratio: None,
        }
    }
}

#[derive(Clone, Debug)]
pub struct VelocityResult {
    /// cells x genes velocities
    pub velocity: RowMat,
    /// cells x genes variance velocities (stochastic only, zero elsewhere)
    pub variance_velocity: Option<RowMat>,
    pub gamma: Vec<f64>,
    pub offset: Vec<f64>,
    pub offset2: Vec<f64>,
    pub r2: Vec<f64>,
    pub qreg_ratio: Vec<f64>,
    pub velocity_genes: Vec<bool>,
    /// Set when too few genes passed and the R² threshold was lowered.
    pub low_signal_r2_threshold: Option<f64>,
}

/// Per-gene extreme quantile linear regression y ~ gamma * x (+ offset).
/// `x`, `y` are the gene's values across cells.
fn quantile_regression(x: &[f64], y: &[f64], fit_offset: bool, perc: Option<(f64, f64)>) -> (f64, f64) {
    let n = x.len();
    let mut keep = vec![true; n];
    let mut n_obs = n as f64;
    if let Some((lo, hi)) = perc {
        let xmax = x.iter().copied().fold(f64::NEG_INFINITY, f64::max).max(1e-3);
        let ymax = y.iter().copied().fold(f64::NEG_INFINITY, f64::max).max(1e-3);
        let normed: Vec<f64> = x.iter().zip(y).map(|(a, b)| a / xmax + b / ymax).collect();
        if fit_offset {
            let lb = percentile(&normed, lo);
            let ub = percentile(&normed, hi);
            for i in 0..n {
                keep[i] = normed[i] <= lb || normed[i] >= ub;
            }
        } else {
            let ub = percentile(&normed, hi);
            for i in 0..n {
                keep[i] = normed[i] >= ub;
            }
        }
        // getnnz of the trimmer: counts True entries
        n_obs = keep.iter().filter(|&&k| k).count() as f64;
    }
    let (mut sxx, mut sxy, mut sx, mut sy) = (0.0, 0.0, 0.0, 0.0);
    for i in 0..n {
        if keep[i] {
            sxx += x[i] * x[i];
            sxy += x[i] * y[i];
            sx += x[i];
            sy += y[i];
        }
    }
    let (mut gamma, mut offset);
    if fit_offset {
        let mx = sx / n_obs;
        let my = sy / n_obs;
        gamma = (sxy / n_obs - mx * my) / (sxx / n_obs - mx * mx);
        offset = my - gamma * mx;
        if offset < 0.0 {
            gamma = sxy / sxx;
            offset = 0.0;
        }
    } else {
        gamma = sxy / sxx;
        offset = 0.0;
    }
    if gamma.is_nan() {
        gamma = 0.0;
    }
    if offset.is_nan() {
        offset = 0.0;
    }
    (gamma, offset)
}

fn r_squared(residual: &[f64], y: &[f64]) -> f64 {
    let my = mean(y);
    let ss_res: f64 = residual.iter().map(|r| r * r).sum();
    let ss_tot: f64 = y.iter().map(|v| (v - my) * (v - my)).sum();
    let r2 = 1.0 - ss_res / ss_tot;
    if r2.is_nan() {
        0.0
    } else {
        r2
    }
}

struct GeneFit {
    gamma: f64,
    offset: f64,
    r2: f64,
}

fn fit_deterministic_gene(x: &[f64], y: &[f64], opts: &SteadyStateOptions) -> GeneFit {
    let (mut gamma, offset) = quantile_regression(x, y, opts.fit_offset, opts.perc);
    if let Some((lo, hi)) = opts.constrain_ratio {
        gamma = gamma.clamp(lo, hi);
    }
    let r2 = if opts.r2_adjusted {
        let (g, o) = quantile_regression(x, y, opts.fit_offset, None);
        let res: Vec<f64> = x.iter().zip(y).map(|(a, b)| b - g * a - if opts.fit_offset { o } else { 0.0 }).collect();
        r_squared(&res, y)
    } else {
        let res: Vec<f64> = x
            .iter()
            .zip(y)
            .map(|(a, b)| b - gamma * a - if opts.fit_offset { offset } else { 0.0 })
            .collect();
        r_squared(&res, y)
    };
    GeneFit { gamma, offset, r2 }
}

/// Steady-state model. `ms`, `mu` are cells x genes first order moments.
/// `highly_variable` optionally restricts velocity genes.
pub fn deterministic(ms: &RowMat, mu: &RowMat, highly_variable: Option<&[bool]>, opts: &SteadyStateOptions) -> VelocityResult {
    let (n, g) = (ms.nrow, ms.ncol);
    let mst = ms.transpose();
    let mut_ = mu.transpose();
    let fits: Vec<GeneFit> = (0..g)
        .into_par_iter()
        .map(|j| fit_deterministic_gene(mst.row(j), mut_.row(j), opts))
        .collect();
    let gamma: Vec<f64> = fits.iter().map(|f| f.gamma).collect();
    let offset: Vec<f64> = fits.iter().map(|f| if opts.fit_offset { f.offset } else { 0.0 }).collect();
    let r2: Vec<f64> = fits.iter().map(|f| f.r2).collect();
    let mut velocity_genes: Vec<bool> = (0..g)
        .map(|j| {
            let any_s = mst.row(j).iter().any(|&v| v > 0.0);
            let any_u = mut_.row(j).iter().any(|&v| v > 0.0);
            r2[j] > opts.min_r2 && gamma[j] > opts.min_ratio && any_s && any_u
        })
        .collect();
    if let Some(hv) = highly_variable {
        for j in 0..g {
            velocity_genes[j] &= hv[j];
        }
    }
    let mut low = None;
    if velocity_genes.iter().filter(|&&b| b).count() < 2 {
        let thr = percentile(&r2, 80.0);
        velocity_genes = r2.iter().map(|&r| r > thr).collect();
        low = Some(thr);
    }
    let mut velocity = RowMat::zeros(n, g);
    velocity.data.par_chunks_mut(g).enumerate().for_each(|(i, row)| {
        let (s, u) = (ms.row(i), mu.row(i));
        for j in 0..g {
            row[j] = u[j] - gamma[j] * s[j] - offset[j];
        }
    });
    VelocityResult {
        velocity,
        variance_velocity: None,
        qreg_ratio: gamma.clone(),
        gamma,
        offset,
        offset2: vec![0.0; g],
        r2,
        velocity_genes,
        low_signal_r2_threshold: low,
    }
}

/// get_weight from scvelo.tools.optimization for a single gene. With `lower`,
/// both tails (<= lower percentile or >= upper percentile) are kept.
fn extreme_weights(x: &[f64], y: Option<&[f64]>, lower: Option<f64>, perc: f64) -> Vec<bool> {
    let normed: Vec<f64> = match y {
        Some(y) => {
            let xmax = x.iter().copied().fold(f64::NEG_INFINITY, f64::max).max(1e-3);
            let ymax = y.iter().copied().fold(f64::NEG_INFINITY, f64::max).max(1e-3);
            x.iter().zip(y).map(|(a, b)| a / xmax + b / ymax).collect()
        }
        None => x.to_vec(),
    };
    let ub = percentile(&normed, perc);
    match lower {
        Some(lo) => {
            let lb = percentile(&normed, lo);
            normed.iter().map(|&v| v <= lb || v >= ub).collect()
        }
        None => normed.iter().map(|&v| v >= ub).collect(),
    }
}

/// Stochastic (second order) model. `mss`, `mus` are second order moments of
/// the same genes (cells x genes); they are only read for velocity genes.
pub fn stochastic(
    ms: &RowMat,
    mu: &RowMat,
    mss: &RowMat,
    mus: &RowMat,
    highly_variable: Option<&[bool]>,
    opts: &SteadyStateOptions,
) -> VelocityResult {
    let mut det = deterministic(ms, mu, highly_variable, opts);
    let (n, g) = (ms.nrow, ms.ncol);
    let idx: Vec<usize> = if det.velocity_genes.iter().any(|&b| b) {
        (0..g).filter(|&j| det.velocity_genes[j]).collect()
    } else {
        (0..g).collect()
    };
    let mst = ms.transpose();
    let mut_ = mu.transpose();
    let msst = mss.transpose();
    let must = mus.transpose();
    let vt = det.velocity.transpose();
    // without an offset only the upper percentile is used (scvelo semantics)
    let perc = opts.perc.map(|(lo, hi)| (if opts.fit_offset { Some(lo) } else { None }, hi));

    let fitted: Vec<(usize, f64, f64, f64)> = idx
        .par_iter()
        .map(|&j| {
            let (x, y) = (mst.row(j), mut_.row(j));
            let var_ss: Vec<f64> = (0..n).map(|i| 2.0 * msst.get(j, i) - x[i]).collect();
            let cov_us: Vec<f64> = (0..n).map(|i| 2.0 * must.get(j, i) + y[i]).collect();
            // OLS on second order moments to initialise the covariance
            let (gamma2, offset2_init) = quantile_regression(&var_ss, &cov_us, opts.fit_offset2, None);
            let res_std = std(vt.row(j));
            let res2: Vec<f64> = (0..n).map(|i| cov_us[i] - gamma2 * var_ss[i] - offset2_init).collect();
            let res2_std = std(&res2);

            let w: Vec<bool> = match perc {
                Some((lo, hi)) => {
                    let a = extreme_weights(x, Some(y), lo, hi);
                    let b = extreme_weights(x, None, lo, hi);
                    a.iter().zip(&b).map(|(a, b)| *a || *b).collect()
                }
                None => vec![true; n],
            };
            // stacked weighted system, solved by pinv(AᵀA)Aᵀy
            let mut a_rows: Vec<[f64; 3]> = Vec::with_capacity(2 * n);
            let mut b_vals: Vec<f64> = Vec::with_capacity(2 * n);
            for i in 0..n {
                let xi = if w[i] { x[i] } else { 0.0 };
                let yi = if w[i] { y[i] } else { 0.0 };
                a_rows.push([1.0 / res_std, 0.0, xi / res_std]);
                b_vals.push(yi / res_std);
            }
            for i in 0..n {
                a_rows.push([0.0, 1.0 / res2_std, var_ss[i] / res2_std]);
                b_vals.push(cov_us[i] / res2_std);
            }
            let cols: Vec<usize> = match (opts.fit_offset, opts.fit_offset2) {
                (true, true) => vec![0, 1, 2],
                (true, false) => vec![0, 2],
                (false, true) => vec![1, 2],
                (false, false) => vec![2],
            };
            let sol = solve_normal_pinv(&a_rows, &b_vals, &cols);
            let mut off = 0.0;
            let mut off2 = 0.0;
            let mut gam = *sol.last().unwrap();
            for (k, &c) in cols.iter().enumerate() {
                match c {
                    0 => off = sol[k],
                    1 => off2 = sol[k],
                    _ => {}
                }
            }
            if off.is_nan() {
                off = 0.0;
            }
            if off2.is_nan() {
                off2 = 0.0;
            }
            if gam.is_nan() {
                gam = 0.0;
            }
            (j, off, off2, gam)
        })
        .collect();

    for &(j, off, off2, gam) in &fitted {
        det.offset[j] = off;
        det.offset2[j] = off2;
        det.gamma[j] = gam;
    }
    let gamma = det.gamma.clone();
    let offset = det.offset.clone();
    let fit_offset = opts.fit_offset;
    det.velocity.data.par_chunks_mut(g).enumerate().for_each(|(i, row)| {
        let (s, u) = (ms.row(i), mu.row(i));
        for j in 0..g {
            row[j] = u[j] - gamma[j] * s[j] - if fit_offset { offset[j] } else { 0.0 };
        }
    });
    let mut var_vel = RowMat::zeros(n, g);
    let offset2 = det.offset2.clone();
    let fit_offset2 = opts.fit_offset2;
    let in_idx: Vec<bool> = {
        let mut v = vec![false; g];
        for &j in &idx {
            v[j] = true;
        }
        v
    };
    var_vel.data.par_chunks_mut(g).enumerate().for_each(|(i, row)| {
        let (s, u, ss, us) = (ms.row(i), mu.row(i), mss.row(i), mus.row(i));
        for j in 0..g {
            if !in_idx[j] {
                continue;
            }
            let var_ss = 2.0 * ss[j] - s[j];
            let cov_us = 2.0 * us[j] + u[j];
            let mut r = (cov_us - 2.0 * s[j] * u[j]) - gamma[j] * (var_ss - 2.0 * s[j] * s[j]);
            if fit_offset {
                r += 2.0 * offset[j] * s[j];
            }
            if fit_offset2 {
                r -= offset2[j];
            }
            row[j] = r;
        }
    });
    det.variance_velocity = Some(var_vel);
    det
}

/// Solve the least squares problem restricted to columns `cols` of the
/// 3-column design matrix via the pseudo-inverse of the normal equations.
fn solve_normal_pinv(a: &[[f64; 3]], b: &[f64], cols: &[usize]) -> Vec<f64> {
    let k = cols.len();
    let mut ata = vec![vec![0.0; k]; k];
    let mut atb = vec![0.0; k];
    for (row, &bv) in a.iter().zip(b) {
        for p in 0..k {
            let ap = row[cols[p]];
            atb[p] += ap * bv;
            for q in 0..k {
                ata[p][q] += ap * row[cols[q]];
            }
        }
    }
    let inv = pinv_small(&ata);
    (0..k).map(|p| (0..k).map(|q| inv[p][q] * atb[q]).sum()).collect()
}

/// Moore-Penrose pseudo-inverse of a small symmetric PSD matrix through a
/// Jacobi eigen decomposition (numpy's pinv cutoff rcond = 1e-15 · max σ).
pub fn pinv_small(m: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let k = m.len();
    if m.iter().flatten().any(|v| !v.is_finite()) {
        return vec![vec![f64::NAN; k]; k];
    }
    let (evals, evecs) = crate::linalg::jacobi_eigen(m);
    let smax = evals.iter().fold(0.0f64, |a, &b| a.max(b.abs()));
    let cutoff = 1e-15 * smax;
    let mut out = vec![vec![0.0; k]; k];
    for l in 0..k {
        if evals[l].abs() > cutoff && evals[l] != 0.0 {
            let inv = 1.0 / evals[l];
            for p in 0..k {
                for q in 0..k {
                    out[p][q] += evecs[p][l] * inv * evecs[q][l];
                }
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deterministic_recovers_ratio() {
        // u = 0.5 s exactly at the upper extreme
        let n = 100;
        let s: Vec<f64> = (0..n).map(|i| i as f64 / 10.0).collect();
        let u: Vec<f64> = s.iter().map(|v| 0.5 * v).collect();
        let ms = RowMat::from_vec(n, 1, s);
        let mu = RowMat::from_vec(n, 1, u);
        let r = deterministic(&ms, &mu, None, &SteadyStateOptions::default());
        assert!((r.gamma[0] - 0.5).abs() < 1e-12);
        assert!(r.velocity.data.iter().all(|v| v.abs() < 1e-12));
    }
}
