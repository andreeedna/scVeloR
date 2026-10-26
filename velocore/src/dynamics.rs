//! Dynamical model of splicing kinetics (scvelo.tl.recover_dynamics,
//! align_dynamics and the dynamical branch of scvelo.tl.velocity).
//!
//! The per-gene expectation-maximisation follows scVelo's `DynamicsRecovery`
//! step by step (same initialisation, same Nelder-Mead schedule, same time
//! assignment rules) so that results are comparable; genes are fitted in
//! parallel.

use crate::matrix::{Csr, RowMat};
use crate::optim::nelder_mead;
use crate::stats::{gaussian_kde_1d, linspace, mean, norm_sf, percentile, std};
use rayon::prelude::*;

#[derive(Clone, Debug)]
pub struct DynamicsOptions {
    /// Total EM iterations, split over 5 optimisation stages (scVelo: 10).
    pub max_iter: usize,
    pub fit_scaling: bool,
    pub fit_steady_states: bool,
    /// Smooth state assignments over the neighbour graph (scVelo default when
    /// fitting on moments).
    pub fit_connected_states: bool,
    /// Orthogonal projection onto the trajectory in the final EM phase.
    pub projection: bool,
    /// Upper percentile for outlier removal in the fit weights.
    pub perc: f64,
}

impl Default for DynamicsOptions {
    fn default() -> Self {
        DynamicsOptions { max_iter: 10, fit_scaling: true, fit_steady_states: true, fit_connected_states: true, projection: true, perc: 99.0 }
    }
}

#[derive(Clone, Debug)]
pub struct GeneDynamics {
    pub recoverable: bool,
    pub alpha: f64,
    /// Splicing rate, *unscaled* (divided by `scaling`, as stored by scVelo).
    pub beta: f64,
    pub gamma: f64,
    pub t_: f64,
    pub scaling: f64,
    pub std_u: f64,
    pub std_s: f64,
    pub likelihood: f64,
    pub variance: f64,
    pub pval_steady: f64,
    pub steady_u: f64,
    pub steady_s: f64,
    pub t: Vec<f64>,
    pub tau: Vec<f64>,
    pub tau_: Vec<f64>,
    pub loss: Vec<f64>,
}

impl GeneDynamics {
    fn unrecoverable(n: usize) -> Self {
        GeneDynamics {
            recoverable: false,
            alpha: f64::NAN,
            beta: f64::NAN,
            gamma: f64::NAN,
            t_: f64::NAN,
            scaling: f64::NAN,
            std_u: f64::NAN,
            std_s: f64::NAN,
            likelihood: f64::NAN,
            variance: f64::NAN,
            pval_steady: f64::NAN,
            steady_u: f64::NAN,
            steady_s: f64::NAN,
            t: vec![f64::NAN; n],
            tau: vec![f64::NAN; n],
            tau_: vec![f64::NAN; n],
            loss: vec![],
        }
    }
}

// ---------------------------------------------------------------------------
// analytic solution of du/dt = alpha - beta u, ds/dt = beta u - gamma s
// ---------------------------------------------------------------------------

#[inline]
fn invert(x: f64) -> f64 {
    if x != 0.0 {
        1.0 / x
    } else {
        f64::NAN
    }
}

#[inline]
pub fn unspliced(tau: f64, u0: f64, alpha: f64, beta: f64) -> f64 {
    let expu = (-beta * tau).exp();
    u0 * expu + alpha / beta * (1.0 - expu)
}

#[inline]
pub fn solution(tau: f64, u0: f64, s0: f64, alpha: f64, beta: f64, gamma: f64) -> (f64, f64) {
    let expu = (-beta * tau).exp();
    let exps = (-gamma * tau).exp();
    let u = u0 * expu + alpha / beta * (1.0 - expu);
    let c = (alpha - u0 * beta) * invert(gamma - beta);
    let s = s0 * exps + alpha / gamma * (1.0 - exps) + c * (exps - expu);
    (u, s)
}

#[inline]
fn np_max(a: f64, b: f64) -> f64 {
    if a.is_nan() || b.is_nan() {
        f64::NAN
    } else {
        a.max(b)
    }
}

#[inline]
fn np_min(a: f64, b: f64) -> f64 {
    if a.is_nan() || b.is_nan() {
        f64::NAN
    } else {
        a.min(b)
    }
}

#[inline]
fn np_clip(x: f64, lo: f64, hi: f64) -> f64 {
    np_min(np_max(x, lo), hi)
}

#[inline]
fn clipped_log(x: f64) -> f64 {
    let eps = 1e-6;
    np_clip(x, eps, 1.0 - eps).ln()
}

/// Inverse of the solution: time needed to reach (u, s) from (u0, s0).
#[inline]
fn tau_inv(u: f64, s: f64, u0: f64, s0: f64, alpha: f64, beta: f64, gamma: f64) -> f64 {
    if gamma >= beta || gamma.is_nan() {
        let uinf = alpha / beta;
        -1.0 / beta * clipped_log((u - uinf) / (u0 - uinf))
    } else {
        let beta_ = beta * invert(gamma - beta);
        let xinf = alpha / gamma - beta_ * (alpha / beta);
        -1.0 / gamma * clipped_log((s - beta_ * u - xinf) / (s0 - beta_ * u0 - xinf))
    }
}

/// tau_inv using only the unspliced equation.
#[inline]
fn tau_inv_u(u: f64, u0: f64, alpha: f64, beta: f64) -> f64 {
    let uinf = alpha / beta;
    -1.0 / beta * clipped_log((u - uinf) / (u0 - uinf))
}

/// scvelo adjust_increments for one (on) or two (on, off) sets of times.
fn adjust_increments(tau: &mut [f64], tau_: Option<&mut [f64]>) {
    fn sorted_diffs(v: &[f64]) -> (Vec<f64>, Vec<f64>) {
        let mut ord = v.to_vec();
        ord.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let mut d = Vec::with_capacity(ord.len());
        let mut prev = 0.0;
        for &x in &ord {
            d.push(x - prev);
            prev = x;
        }
        (ord, d)
    }
    fn apply(orig: &[f64], new: &mut [f64], ord: &[f64], d: &[f64], ub: f64) {
        for i in 0..d.len() {
            if d[i] > ub {
                let (ti, dti) = (ord[i], d[i]);
                for k in 0..orig.len() {
                    if orig[k] >= ti {
                        new[k] -= dti;
                    }
                }
            }
        }
    }
    let orig = tau.to_vec();
    let (ord, d) = sorted_diffs(&orig);
    match tau_ {
        None => {
            let ub = 3.0 * percentile(&d, 99.5);
            apply(&orig, tau, &ord, &d, ub);
        }
        Some(tau_) => {
            let orig_ = tau_.to_vec();
            let (ord_, d_) = sorted_diffs(&orig_);
            let mut both = d.clone();
            both.extend_from_slice(&d_);
            let ub = 3.0 * percentile(&both, 99.5);
            apply(&orig_, tau_, &ord_, &d_, ub);
            apply(&orig, tau, &ord, &d, ub);
        }
    }
}

/// scvelo test_bimodality (kde=True, bins=30). Returns (pval, means).
fn test_bimodality(x: &[f64]) -> (f64, [f64; 2]) {
    let bins = 30usize;
    let lb = x.iter().copied().fold(f64::INFINITY, f64::min);
    let ub = percentile(x, 99.9);
    let hi = if ub <= lb { ub } else { x.iter().copied().fold(f64::NEG_INFINITY, f64::max) };
    let grid = linspace(lb, hi, bins);
    let kde = match gaussian_kde_1d(x, &grid) {
        Some(k) => k,
        None => return (1.0, [0.0, 0.0]),
    };
    let mut idx = bins / 2 - 2;
    let end = idx + 4;
    idx += crate::stats::argmin(&kde[idx..end]);
    let peak_0 = crate::stats::argmax(&kde[..idx]);
    let peak_1 = crate::stats::argmax(&kde[idx..]);
    let kde_peak = kde[idx + peak_1];
    let kde_mid = mean(&kde[idx..]);
    let t_stat = (kde_peak - kde_mid) / (std(&kde) / (bins as f64).sqrt()).max(1.0);
    let p_val = norm_sf(t_stat);
    let g0 = &grid[..idx];
    let g1 = &grid[idx..];
    let means = [
        (g0[peak_0] + g0[(peak_0 + 1).min(g0.len() - 1)]) / 2.0,
        (g1[peak_1] + g1[(peak_1 + 1).min(g1.len() - 1)]) / 2.0,
    ];
    (p_val, means)
}

#[derive(Clone, Copy, Debug)]
struct Pars {
    alpha: f64,
    beta: f64,
    gamma: f64,
    scaling: f64,
    t_: f64,
}

#[derive(Clone, Copy, PartialEq)]
enum Weighting {
    Fit,
    Upper,
}

struct Timing {
    t: Vec<f64>,
    tau: Vec<f64>,
    o: Vec<u8>,
}

struct Model<'a> {
    u: Vec<f64>,
    s: Vec<f64>,
    n: usize,
    weights: Vec<bool>,
    weights_upper: Vec<bool>,
    std_u: f64,
    std_s: f64,
    conn: Option<&'a Csr>,
    fit_steady_states: bool,
    projection: bool,
    p: Pars,
    t: Vec<f64>,
    tau: Vec<f64>,
    o: Vec<u8>,
    loss: Vec<f64>,
    steady_state_ratio: Option<f64>,
    steady_u: f64,
    steady_s: f64,
    pval_steady: f64,
    simplex_maxiter: usize,
}

impl<'a> Model<'a> {
    fn smooth(&self, v: Vec<f64>) -> Vec<f64> {
        match self.conn {
            Some(c) => c.dot_vec(&v),
            None => v,
        }
    }

    fn resolve(&self, alpha: Option<f64>, beta: Option<f64>, gamma: Option<f64>, scaling: Option<f64>, t_: Option<f64>) -> Pars {
        let t_ = match t_ {
            Some(v) if v != 0.0 => v,
            _ => self.p.t_,
        };
        Pars {
            alpha: alpha.unwrap_or(self.p.alpha),
            beta: beta.unwrap_or(self.p.beta),
            gamma: gamma.unwrap_or(self.p.gamma),
            scaling: scaling.unwrap_or(self.p.scaling),
            t_,
        }
    }

    /// assign_tau: (tau, tau_) for all cells. `u` is already divided by scaling.
    fn assign_tau(&self, u: &[f64], p: &Pars, u0_: f64, s0_: f64, projection: bool) -> (Vec<f64>, Vec<f64>) {
        let n = self.n;
        let s = &self.s;
        if projection && p.beta < p.gamma {
            let umin = (0..n).filter(|&i| s[i] > 0.0).map(|i| u[i]).fold(f64::INFINITY, f64::min);
            let t0 = tau_inv_u(umin, u0_, 0.0, p.beta);
            let num = ((n / 5) as usize).clamp(200, 500);
            let tpoints = linspace(0.0, p.t_, num);
            let tpoints_: Vec<f64> = linspace(0.0, t0, num)[1..].to_vec();
            let xt: Vec<(f64, f64)> = tpoints.iter().map(|&t| solution(t, 0.0, 0.0, p.alpha, p.beta, p.gamma)).collect();
            let xt_: Vec<(f64, f64)> = tpoints_.iter().map(|&t| solution(t, u0_, s0_, 0.0, p.beta, p.gamma)).collect();
            let nearest = |pts: &[(f64, f64)], ui: f64, si: f64| -> usize {
                let mut best = 0;
                let mut bd = f64::INFINITY;
                for (k, &(a, b)) in pts.iter().enumerate() {
                    let d = (a - ui) * (a - ui) + (b - si) * (b - si);
                    if d.is_nan() {
                        return k;
                    }
                    if d < bd {
                        bd = d;
                        best = k;
                    }
                }
                best
            };
            let tau: Vec<f64> = (0..n).map(|i| tpoints[nearest(&xt, u[i], s[i])]).collect();
            let tau_: Vec<f64> = if tpoints_.is_empty() {
                vec![0.0; n]
            } else {
                (0..n).map(|i| tpoints_[nearest(&xt_, u[i], s[i])]).collect()
            };
            (tau, tau_)
        } else {
            let tau: Vec<f64> = (0..n)
                .map(|i| np_clip(tau_inv(u[i], s[i], 0.0, 0.0, p.alpha, p.beta, p.gamma), 0.0, p.t_))
                .collect();
            let tau_raw: Vec<f64> = (0..n).map(|i| tau_inv(u[i], s[i], u0_, s0_, 0.0, p.beta, p.gamma)).collect();
            let mut mx = f64::NEG_INFINITY;
            for i in 0..n {
                if s[i] > 0.0 {
                    mx = np_max(mx, tau_raw[i]);
                }
            }
            let tau_ = tau_raw.iter().map(|&v| np_clip(v, 0.0, mx)).collect();
            (tau, tau_)
        }
    }

    /// compute_divergence(mode="assign_timepoints") (and "tau" via `taus_only`).
    fn divergence_time(&self, p: &Pars, projection: bool) -> (Timing, Vec<f64>, Vec<f64>) {
        let n = self.n;
        let u: Vec<f64> = self.u.iter().map(|v| v / p.scaling).collect();
        let s = &self.s;
        let (u0_, s0_) = solution(p.t_, 0.0, 0.0, p.alpha, p.beta, p.gamma);
        let (mut tau, mut tau_) = self.assign_tau(&u, p, u0_, s0_, projection);
        let std_u = self.std_u / p.scaling;
        let std_s = self.std_s;

        // constrain time increments
        {
            let mut r_off = vec![0.0; n];
            let mut r_on = vec![0.0; n];
            for i in 0..n {
                let (ut, st) = solution(tau[i], 0.0, 0.0, p.alpha, p.beta, p.gamma);
                let (ut_, st_) = solution(tau_[i], u0_, s0_, 0.0, p.beta, p.gamma);
                let (du, du_) = ((u[i] - ut) / std_u, (u[i] - ut_) / std_u);
                let (ds, ds_) = ((s[i] - st) / std_s, (s[i] - st_) / std_s);
                r_off[i] = du_ * du_ + ds_ * ds_;
                r_on[i] = du * du + ds * ds;
            }
            let r_off = self.smooth(r_off);
            let r_on = self.smooth(r_on);
            let on: Vec<usize> = (0..n).filter(|&i| state_argmin(&[r_off[i], r_on[i]]) == 1).collect();
            let off: Vec<usize> = (0..n).filter(|&i| state_argmin(&[r_off[i], r_on[i]]) == 0).collect();
            let mut t_on: Vec<f64> = on.iter().map(|&i| tau[i]).collect();
            let mut t_off: Vec<f64> = off.iter().map(|&i| tau_[i]).collect();
            if !on.is_empty() && !off.is_empty() {
                adjust_increments(&mut t_on, Some(&mut t_off));
            } else if !on.is_empty() {
                adjust_increments(&mut t_on, None);
            } else if !off.is_empty() {
                adjust_increments(&mut t_off, None);
            }
            for (k, &i) in on.iter().enumerate() {
                tau[i] = t_on[k];
            }
            for (k, &i) in off.iter().enumerate() {
                tau_[i] = t_off[k];
            }
        }

        let mut res: Vec<Vec<f64>> = vec![vec![0.0; n]; if self.fit_steady_states { 4 } else { 2 }];
        let (u_inf, s_inf) = (p.alpha / p.beta, p.alpha / p.gamma);
        for i in 0..n {
            let (ut, st) = solution(tau[i], 0.0, 0.0, p.alpha, p.beta, p.gamma);
            let (ut_, st_) = solution(tau_[i], u0_, s0_, 0.0, p.beta, p.gamma);
            let (du, du_) = ((u[i] - ut) / std_u, (u[i] - ut_) / std_u);
            let (ds, ds_) = ((s[i] - st) / std_s, (s[i] - st_) / std_s);
            res[0][i] = du_ * du_ + ds_ * ds_;
            res[1][i] = du * du + ds * ds;
            if self.fit_steady_states {
                res[2][i] = (u[i] / std_u).powi(2) + (s[i] / std_s).powi(2);
                res[3][i] = ((u[i] - u_inf) / std_u).powi(2) + ((s[i] - s_inf) / std_s).powi(2);
            }
        }
        let res: Vec<Vec<f64>> = res.into_iter().map(|r| self.smooth(r)).collect();
        let taus = (tau.clone(), tau_.clone());
        let mut t = vec![0.0; n];
        let mut o = vec![0u8; n];
        for i in 0..n {
            let vals: Vec<f64> = res.iter().map(|r| r[i]).collect();
            let st = state_argmin(&vals);
            if st != 0 {
                tau_[i] = 0.0;
            }
            if st != 1 {
                tau[i] = 0.0;
            }
            let oo = match st {
                2 => 1u8,
                3 => 0u8,
                x => x as u8,
            };
            o[i] = oo;
            t[i] = if oo == 1 { tau[i] } else { tau_[i] + p.t_ };
        }
        (Timing { t, tau, o }, taus.0, taus.1)
    }

    fn time_assignment(&self, p: &Pars, refit: bool, given_t: Option<&[f64]>) -> Timing {
        if let Some(t) = given_t {
            let o: Vec<u8> = t.iter().map(|&v| (v < p.t_) as u8).collect();
            let tau = t.iter().zip(&o).map(|(&v, &oo)| if oo == 1 { v } else { v - p.t_ }).collect();
            Timing { t: t.to_vec(), tau, o }
        } else if refit {
            self.divergence_time(p, self.projection).0
        } else {
            Timing { t: self.t.clone(), tau: self.tau.clone(), o: self.o.clone() }
        }
    }

    fn weights(&self, w: Weighting) -> Vec<bool> {
        match w {
            Weighting::Fit => self.weights.clone(),
            Weighting::Upper => self.weights_upper.clone(),
        }
    }

    /// get_dists: (udiff, sdiff, reg) on the weighted cells.
    fn dists(&self, p: &Pars, t: &[f64], w: Weighting) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let wt = self.weights(w);
        let (u0_, s0_) = solution(p.t_, 0.0, 0.0, p.alpha, p.beta, p.gamma);
        let mut ud = Vec::new();
        let mut sd = Vec::new();
        let mut rg = Vec::new();
        for i in 0..self.n {
            if !wt[i] {
                continue;
            }
            let on = t[i] < p.t_;
            let (tau, a, u0, s0) = if on { (t[i], p.alpha, 0.0, 0.0) } else { (t[i] - p.t_, 0.0, u0_, s0_) };
            let (ut, st) = solution(tau, u0, s0, a, p.beta, p.gamma);
            let ui = self.u[i] / p.scaling;
            let si = self.s[i];
            ud.push((ut - ui) / self.std_u * p.scaling);
            sd.push((st - si) / self.std_s);
            rg.push(match self.steady_state_ratio {
                Some(r) => (p.gamma / p.beta - r) * si / self.std_s,
                None => 0.0,
            });
        }
        (ud, sd, rg)
    }

    fn se(&self, p: &Pars, t: &[f64]) -> f64 {
        let (ud, sd, rg) = self.dists(p, t, Weighting::Fit);
        ud.iter().zip(&sd).zip(&rg).map(|((a, b), c)| a * a + b * b + c * c).sum()
    }

    /// get_mse with time re-assignment (or fixed times if `refit` is false).
    fn mse(&self, p: &Pars, refit: bool) -> f64 {
        let tm = self.time_assignment(p, refit, None);
        let (ud, sd, rg) = self.dists(p, &tm.t, Weighting::Fit);
        let n = ud.len() as f64;
        ud.iter().zip(&sd).zip(&rg).map(|((a, b), c)| a * a + b * b + c * c).sum::<f64>() / n
    }

    fn update(&mut self, alpha: Option<f64>, beta: Option<f64>, gamma: Option<f64>, scaling: Option<f64>, t_: Option<f64>, adjust_t_: bool) -> bool {
        let loss_prev = *self.loss.last().unwrap_or(&1e6);
        let mut p = self.resolve(alpha, beta, gamma, scaling, t_);
        let mut tm = self.time_assignment(&p, true, None);
        let mut loss = self.se(&p, &tm.t);
        let mut perform = loss < loss_prev;
        let on: Vec<bool> = self.o.iter().map(|&v| v == 1).collect();
        if adjust_t_ && on.iter().any(|&b| b) {
            if !perform {
                p = self.p;
                tm = self.time_assignment(&p, true, None);
                loss = self.se(&p, &tm.t);
            }
            let mut alt_t_ = f64::NEG_INFINITY;
            for i in 0..self.n {
                if on[i] {
                    alt_t_ = np_max(alt_t_, tm.t[i]);
                }
            }
            if 0.0 < alt_t_ && alt_t_ < p.t_ {
                let tmax = tm.t.iter().copied().fold(f64::NEG_INFINITY, np_max);
                let n_eq = tm.t.iter().filter(|&&v| v == p.t_).count() as f64;
                alt_t_ += tmax / self.n as f64 * n_eq;
                let mut pa = p;
                pa.t_ = alt_t_;
                let tma = self.time_assignment(&pa, true, None);
                let alt_loss = self.se(&pa, &tma.t);
                let ut_cur = unspliced(p.t_, 0.0, p.alpha, p.beta);
                let ut_alt = unspliced(alt_t_, 0.0, p.alpha, p.beta);
                let min_loss = loss.min(loss_prev);
                if std::env::var("VELOCORE_TRACE").is_ok() {
                    eprintln!("  alt: t_={} alt_t_={} tmax={} n_eq={} loss={} prev={} alt_loss={} ut_cur={} ut_alt={}", p.t_, alt_t_, tmax, n_eq, loss, loss_prev, alt_loss, ut_cur, ut_alt);
                }
                if alt_loss * 0.99 <= min_loss || ut_cur * 0.99 < ut_alt {
                    p = pa;
                    tm = tma;
                    loss = alt_loss;
                    perform = true;
                }
            }
        }
        if perform {
            let ratio = self.p.scaling / p.scaling;
            self.steady_u *= ratio;
            self.t = tm.t;
            self.tau = tm.tau;
            self.o = tm.o;
            self.p = p;
            self.loss.push(loss);
        }
        perform
    }

    fn simplex(&self, x0: &[f64], tol: Option<f64>, f: impl Fn(&[f64]) -> f64) -> Vec<f64> {
        let tol = tol.unwrap_or(1e-4);
        nelder_mead(f, x0, self.simplex_maxiter, tol, tol)
    }

    fn fit_t_and_alpha(&mut self) {
        let a0 = self.p.alpha;
        for k in linspace(-1.0, 1.0, 5) {
            self.update(Some(a0 + k * a0 / 10.0), None, None, None, None, true);
        }
        let x0 = [self.p.t_, self.p.alpha];
        let x = self.simplex(&x0, None, |x| {
            let p = self.resolve(Some(x[1]), None, None, None, Some(x[0]));
            self.mse(&p, true)
        });
        self.update(Some(x[1]), None, None, None, Some(x[0]), true);
    }

    fn fit_scaling_(&mut self) {
        let x0 = [self.p.t_, self.p.beta, self.p.scaling];
        let x = self.simplex(&x0, None, |x| {
            let p = self.resolve(None, Some(x[1]), None, Some(x[2]), Some(x[0]));
            self.mse(&p, true)
        });
        self.update(None, Some(x[1]), None, Some(x[2]), Some(x[0]), true);
    }

    fn fit_rates(&mut self) {
        let x0 = [self.p.alpha, self.p.gamma];
        let x = self.simplex(&x0, Some(1e-2), |x| {
            let p = self.resolve(Some(x[0]), None, Some(x[1]), None, None);
            self.mse(&p, true)
        });
        self.update(Some(x[0]), None, Some(x[1]), None, None, true);
    }

    fn fit_t_(&mut self) {
        let x0 = [self.p.t_];
        let x = self.simplex(&x0, None, |x| {
            let p = self.resolve(None, None, None, None, Some(x[0]));
            self.mse(&p, true)
        });
        self.update(None, None, None, None, Some(x[0]), true);
    }

    fn fit_t_and_rates(&mut self, refit: bool) {
        let x0 = [self.p.t_, self.p.alpha, self.p.beta, self.p.gamma];
        let x = self.simplex(&x0, Some(1e-2), |x| {
            let p = self.resolve(Some(x[1]), Some(x[2]), Some(x[3]), None, Some(x[0]));
            self.mse(&p, refit)
        });
        self.update(Some(x[1]), Some(x[2]), Some(x[3]), None, Some(x[0]), true);
    }

    fn loglikelihood_upper(&self, refit: bool) -> f64 {
        let tm = self.time_assignment(&self.p, refit, None);
        let (ud, sd, rg) = self.dists(&self.p, &tm.t, Weighting::Upper);
        let distx: Vec<f64> = ud.iter().zip(&sd).zip(&rg).map(|((a, b), c)| a * a + b * b + c * c).collect();
        let n = (distx.len() as f64 - self.n as f64 * 0.01).max(2.0);
        let m1 = mean(&distx);
        let m2 = mean(&distx.iter().zip(&sd).map(|(d, s)| sign(*s) * d.sqrt()).collect::<Vec<_>>());
        let mut varx = m1 - m2 * m2;
        if varx == 0.0 {
            varx += 1.0;
        }
        -0.5 / n * distx.iter().sum::<f64>() / varx - 0.5 * (2.0 * std::f64::consts::PI * varx).ln()
    }

    fn variance_upper(&self) -> f64 {
        let tm = self.time_assignment(&self.p, true, None);
        let (ud, sd, _) = self.dists(&self.p, &tm.t, Weighting::Upper);
        let distx: Vec<f64> = ud.iter().zip(&sd).map(|(a, b)| a * a + b * b).collect();
        let m2 = mean(&distx.iter().zip(&sd).map(|(d, s)| sign(*s) * d.sqrt()).collect::<Vec<_>>());
        mean(&distx) - m2 * m2
    }
}

#[inline]
fn sign(v: f64) -> f64 {
    if v > 0.0 {
        1.0
    } else if v < 0.0 {
        -1.0
    } else {
        0.0
    }
}

/// numpy argmin over a short list of states (first minimum wins, NaN wins).
#[inline]
fn state_argmin(v: &[f64]) -> usize {
    crate::stats::argmin(v)
}

/// Fit the dynamical model for a single gene.
/// `u`, `s`: unspliced / spliced moments of the gene across cells.
pub fn fit_gene(u: &[f64], s: &[f64], conn: Option<&Csr>, opts: &DynamicsOptions) -> GeneDynamics {
    let n = u.len();
    // --- BaseDynamics.__init__ / initialize_weights
    let nonzero: Vec<bool> = (0..n).map(|i| s[i] > 0.0 && u[i] > 0.0).collect();
    if nonzero.iter().filter(|&&b| b).count() <= 2 {
        return GeneDynamics::unrecoverable(n);
    }
    let mut weights = nonzero.clone();
    let sw: Vec<f64> = (0..n).filter(|&i| weights[i]).map(|i| s[i]).collect();
    let uw: Vec<f64> = (0..n).filter(|&i| weights[i]).map(|i| u[i]).collect();
    let ub_s = percentile(&sw, opts.perc);
    let ub_u = percentile(&uw, opts.perc);
    for i in 0..n {
        if ub_s > 0.0 {
            weights[i] &= s[i] <= ub_s;
        }
        if ub_u > 0.0 {
            weights[i] &= u[i] <= ub_u;
        }
    }
    let u_w: Vec<f64> = (0..n).filter(|&i| weights[i]).map(|i| u[i]).collect();
    let s_w: Vec<f64> = (0..n).filter(|&i| weights[i]).map(|i| s[i]).collect();
    let mut weights_upper = weights.clone();
    if !u_w.is_empty() {
        let umax = u_w.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let smax = s_w.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        for i in 0..n {
            weights_upper[i] &= u[i] > umax / 3.0 && s[i] > smax / 3.0;
        }
    }

    // --- DynamicsRecovery.initialize
    let (mut std_u, mut std_s) = (std(&u_w), std(&s_w));
    if std_u == 0.0 || std_s == 0.0 {
        std_u = 1.0;
        std_s = 1.0;
    }
    let scaling = if opts.fit_scaling { std_u / std_s } else { 1.0 };
    let u_ws: Vec<f64> = u_w.iter().map(|v| v / scaling).collect();
    let perc = 98.0;
    let ps = percentile(&s_w, perc);
    let pu = percentile(&u_ws, perc);
    let weights_s: Vec<bool> = s_w.iter().map(|&v| v >= ps).collect();
    let weights_u: Vec<bool> = u_ws.iter().map(|&v| v >= pu).collect();
    let beta = 1.0;
    let (mut sus, mut sss) = (0.0, 0.0);
    for k in 0..s_w.len() {
        if weights_s[k] {
            sus += s_w[k] * u_ws[k];
            sss += s_w[k] * s_w[k];
        }
    }
    let mut gamma = sus / sss + 1e-6;
    if gamma < 0.05 / scaling {
        gamma *= 1.2;
    } else if gamma > 1.5 / scaling {
        gamma /= 1.2;
    }
    let mut u_inf = mean(&(0..u_ws.len()).filter(|&k| weights_u[k] || weights_s[k]).map(|k| u_ws[k]).collect::<Vec<_>>());
    let s_inf = mean(&(0..s_w.len()).filter(|&k| weights_s[k]).map(|k| s_w[k]).collect::<Vec<_>>());
    let (mut u0_, mut s0_) = (u_inf, s_inf);
    let mut alpha = u_inf * beta;
    let mut beta = beta;
    let (pval_u, means_u) = test_bimodality(&u_ws);
    let (pval_s, means_s) = test_bimodality(&s_w);
    let pval_steady = pval_u.max(pval_s);
    let steady_u = means_u[1];
    let steady_s = means_s[1];
    if pval_steady < 1e-3 {
        u_inf = (u_inf + steady_u) / 2.0;
        alpha = gamma * s_inf;
        beta = alpha / u_inf;
        u0_ = u_inf;
        s0_ = s_inf;
    }
    let t_ = tau_inv(u0_, s0_, 0.0, 0.0, alpha, beta, gamma);

    let mut m = Model {
        u: u.to_vec(),
        s: s.to_vec(),
        n,
        weights,
        weights_upper,
        std_u,
        std_s,
        conn: if opts.fit_connected_states { conn } else { None },
        fit_steady_states: opts.fit_steady_states,
        projection: false,
        p: Pars { alpha, beta, gamma, scaling, t_ },
        t: vec![],
        tau: vec![],
        o: vec![],
        loss: vec![],
        steady_state_ratio: None,
        steady_u,
        steady_s,
        pval_steady,
        simplex_maxiter: opts.max_iter / 5,
    };
    let tm = m.time_assignment(&m.p, true, None);
    m.t = tm.t;
    m.tau = tm.tau;
    m.o = tm.o;
    let l0 = m.se(&m.p, &m.t.clone());
    m.loss.push(l0);
    if opts.fit_scaling {
        for sight in [0.5, 0.1] {
            let sc = m.p.scaling;
            for k in linspace(-1.0, 1.0, 4) {
                let z = sc + k * sc * sight;
                let b = m.p.beta / m.p.scaling * z;
                m.update(None, Some(b), None, Some(z), None, true);
            }
        }
    }
    m.steady_state_ratio = Some(m.p.gamma / m.p.beta);

    // set VELOCORE_TRACE=1 to print the EM path (debugging aid)
    let trace = std::env::var("VELOCORE_TRACE").is_ok();
    let show = |m: &Model, stage: &str| {
        if trace {
            eprintln!(
                "{stage} {} {} {} {} {} loss {} nloss {}",
                m.p.alpha, m.p.beta, m.p.gamma, m.p.t_, m.p.scaling, m.loss.last().unwrap(), m.loss.len()
            );
        }
    };
    if trace {
        eprintln!("init pars {} {} {} {} {} pval {}", alpha, beta, gamma, t_, scaling, pval_steady);
        eprintln!("init loss {:?}", m.loss);
    }

    // --- DynamicsRecovery.fit
    if opts.max_iter > 0 {
        m.fit_t_and_alpha();
        show(&m, "t_and_alpha");
        if opts.fit_scaling {
            m.fit_scaling_();
            show(&m, "scaling");
        }
        m.fit_rates();
        show(&m, "rates");
        m.fit_t_();
        show(&m, "t_");
        m.fit_t_and_rates(true);
        show(&m, "t_and_rates");
        m.projection = opts.projection;
        m.update(None, None, None, None, None, false);
        show(&m, "proj update");
        m.fit_t_and_rates(false);
        show(&m, "t_and_rates fixed");
    }
    m.update(None, None, None, None, None, true);
    show(&m, "final");
    let (_, tau, tau_) = m.divergence_time(&m.p, m.projection);
    let likelihood = m.loglikelihood_upper(false).exp();
    let variance = m.variance_upper();

    GeneDynamics {
        recoverable: true,
        alpha: m.p.alpha,
        beta: m.p.beta / m.p.scaling,
        gamma: m.p.gamma,
        t_: m.p.t_,
        scaling: m.p.scaling,
        std_u: m.std_u,
        std_s: m.std_s,
        likelihood,
        variance,
        pval_steady: m.pval_steady,
        steady_u: m.steady_u * m.p.scaling,
        steady_s: m.steady_s,
        t: m.t,
        tau,
        tau_,
        loss: m.loss,
    }
}

/// Fitted parameters for many genes plus the cell x gene time matrices.
#[derive(Clone, Debug)]
pub struct DynamicsFit {
    pub genes: Vec<GeneDynamics>,
    /// cells x genes latent time (smoothed over connectivities), tau, tau_
    pub fit_t: RowMat,
    pub fit_tau: RowMat,
    pub fit_tau_: RowMat,
    /// per gene time alignment factor (NaN where not applicable)
    pub alignment_scaling: Vec<f64>,
}

/// recover_dynamics + align_dynamics for the selected genes.
/// `mu`, `ms`: cells x genes moments (only the genes to fit).
/// `t_max`: total time for alignment (scVelo default 20); None disables alignment.
pub fn recover_dynamics(mu: &RowMat, ms: &RowMat, conn: Option<&Csr>, opts: &DynamicsOptions, t_max: Option<f64>) -> DynamicsFit {
    let (n, g) = (ms.nrow, ms.ncol);
    let mut_ = mu.transpose();
    let mst = ms.transpose();
    let mut genes: Vec<GeneDynamics> = (0..g)
        .into_par_iter()
        .map(|j| fit_gene(mut_.row(j), mst.row(j), conn, opts))
        .collect();

    // cells x genes time matrices
    let mut tt = RowMat::zeros(g, n);
    let mut tau = RowMat::zeros(g, n);
    let mut tau_ = RowMat::zeros(g, n);
    for (j, gd) in genes.iter().enumerate() {
        tt.row_mut(j).copy_from_slice(&gd.t);
        tau.row_mut(j).copy_from_slice(&gd.tau);
        tau_.row_mut(j).copy_from_slice(&gd.tau_);
    }
    // fit_t = conn · T
    if let Some(c) = conn {
        let rows: Vec<Vec<f64>> = (0..g).into_par_iter().map(|j| c.dot_vec(tt.row(j))).collect();
        for (j, r) in rows.into_iter().enumerate() {
            tt.row_mut(j).copy_from_slice(&r);
        }
    }

    let mut alignment = vec![f64::NAN; g];
    if let Some(tmax) = t_max {
        for j in 0..g {
            let gd = &mut genes[j];
            if !gd.recoverable {
                continue;
            }
            let t = tt.row(j);
            let t_ = gd.t_;
            let mut on_max = f64::NEG_INFINITY;
            let mut off_max = f64::NEG_INFINITY;
            let mut n_eq = 0usize;
            for &v in t {
                on_max = np_max(on_max, v * if v < t_ { 1.0 } else { 0.0 });
                off_max = np_max(off_max, (v - t_) * if v > t_ { 1.0 } else { 0.0 });
                if v == t_ || v == 0.0 {
                    n_eq += 1;
                }
            }
            let mut tm = on_max + off_max;
            let mut denom = 1.0 - n_eq as f64 / n as f64;
            if denom == 0.0 {
                denom += 1.0;
            }
            tm /= denom;
            if tm == 0.0 {
                tm += 1.0;
            }
            let m = tmax / tm;
            gd.alpha /= m;
            gd.beta /= m;
            gd.gamma /= m;
            gd.t_ *= m;
            for v in tt.row_mut(j) {
                *v *= m;
            }
            for v in tau.row_mut(j) {
                *v *= m;
            }
            for v in tau_.row_mut(j) {
                *v *= m;
            }
            gd.t.iter_mut().for_each(|v| *v *= m);
            gd.tau.iter_mut().for_each(|v| *v *= m);
            gd.tau_.iter_mut().for_each(|v| *v *= m);
            alignment[j] = if m == 1.0 { f64::NAN } else { m };
        }
    }
    DynamicsFit { genes, fit_t: tt.transpose(), fit_tau: tau.transpose(), fit_tau_: tau_.transpose(), alignment_scaling: alignment }
}

/// Velocities of the dynamical model (scvelo.tl.velocity(mode="dynamical")).
/// Inputs are cells x genes moments and the per-gene fitted parameters
/// (unscaled beta as stored by [`recover_dynamics`]).
/// Returns (velocity of s, velocity of u).
#[allow(clippy::too_many_arguments)]
pub fn dynamical_velocity(
    mu: &RowMat,
    ms: &RowMat,
    fit: &DynamicsFit,
    conn: Option<&Csr>,
) -> (RowMat, RowMat) {
    let (n, g) = (ms.nrow, ms.ncol);
    let mut_ = mu.transpose();
    let mst = ms.transpose();
    let taut = fit.fit_tau.transpose();
    let tau_t = fit.fit_tau_.transpose();
    let cols: Vec<(Vec<f64>, Vec<f64>)> = (0..g)
        .into_par_iter()
        .map(|j| {
            let gd = &fit.genes[j];
            if !gd.recoverable || gd.alpha.is_nan() {
                return (vec![f64::NAN; n], vec![f64::NAN; n]);
            }
            let scaling = gd.scaling;
            let (alpha, beta, gamma, t_) = (gd.alpha, gd.beta * scaling, gd.gamma, gd.t_);
            let u: Vec<f64> = mut_.row(j).iter().map(|v| v / scaling).collect();
            let s = mst.row(j);
            let tau = taut.row(j);
            let tau_ = tau_t.row(j);
            let std_u = gd.std_u / scaling;
            let std_s = gd.std_s;
            let (u0_, s0_) = solution(t_, 0.0, 0.0, alpha, beta, gamma);
            let (u_inf, s_inf) = (alpha / beta, alpha / gamma);
            let mut res: Vec<Vec<f64>> = vec![vec![0.0; n]; 4];
            let mut dist = vec![0.0; n];
            let mut sgn = vec![0.0; n];
            for i in 0..n {
                let (ut, st) = solution(tau[i], 0.0, 0.0, alpha, beta, gamma);
                let (ut_, st_) = solution(tau_[i], u0_, s0_, 0.0, beta, gamma);
                let (du, du_) = ((u[i] - ut) / std_u, (u[i] - ut_) / std_u);
                let (ds, ds_) = ((s[i] - st) / std_s, (s[i] - st_) / std_s);
                let dx = du * du + ds * ds;
                let dx_ = du_ * du_ + ds_ * ds_;
                res[0][i] = dx_;
                res[1][i] = dx;
                res[2][i] = (u[i] / std_u).powi(2) + (s[i] / std_s).powi(2);
                res[3][i] = ((u[i] - u_inf) / std_u).powi(2) + ((s[i] - s_inf) / std_s).powi(2);
                let o = state_argmin(&[dx_, dx]);
                if o == 1 {
                    dist[i] = dx;
                    sgn[i] = sign(ds);
                } else {
                    dist[i] = dx_;
                    sgn[i] = sign(ds_);
                }
            }
            let res: Vec<Vec<f64>> = match conn {
                Some(c) => res.into_iter().map(|r| c.dot_vec(&r)).collect(),
                None => res,
            };
            let m2 = mean(&dist.iter().zip(&sgn).map(|(d, s)| s * d.sqrt()).collect::<Vec<_>>());
            let varx = mean(&dist) - m2 * m2;
            let norm = 1.0 / (2.0 * std::f64::consts::PI * varx.sqrt());
            let mut state = vec![0usize; n];
            let mut t = vec![0.0; n];
            for i in 0..n {
                // argmax over [l_off, l_on, l_zero, l_steady, 1e-6]; normalising
                // by the total (scVelo `normalized=True`) does not change it
                let vals = [
                    norm * (-0.5 * res[0][i] / varx).exp(),
                    norm * (-0.5 * res[1][i] / varx).exp(),
                    norm * (-0.5 * res[2][i] / varx).exp(),
                    norm * (-0.5 * res[3][i] / varx).exp(),
                    1e-6,
                ];
                let best = crate::stats::argmax(&vals);
                state[i] = best;
                t[i] = match best {
                    1 => tau[i],
                    0 => tau_[i] + t_,
                    _ => 0.0,
                };
            }
            let t = match conn {
                Some(c) => c.dot_vec(&t),
                None => t,
            };
            let mut vs = vec![0.0; n];
            let mut vu = vec![0.0; n];
            for i in 0..n {
                let o = state[i] < 2 && t[i] < t_;
                let o_ = state[i] < 2 && t[i] >= t_;
                let on = t[i] < t_;
                let (tv, av, u0, s0) = if on { (t[i], alpha, 0.0, 0.0) } else { (t[i] - t_, 0.0, u0_, s0_) };
                let (ut, st) = solution(tv, u0, s0, av, beta, gamma);
                let w = if o || o_ { 1.0 } else { 0.0 };
                let (ut, st) = (ut * w, st * w);
                let a = av * if o { 1.0 } else { 0.0 };
                let vt = ut * beta - st * gamma;
                let wt = (a - beta * ut) * scaling;
                vs[i] = np_max(vt, -s[i]);
                vu[i] = np_max(wt, -u[i] * scaling);
            }
            (vs, vu)
        })
        .collect();
    let mut vs = RowMat::zeros(g, n);
    let mut vu = RowMat::zeros(g, n);
    for (j, (a, b)) in cols.into_iter().enumerate() {
        vs.row_mut(j).copy_from_slice(&a);
        vu.row_mut(j).copy_from_slice(&b);
    }
    (vs.transpose(), vu.transpose())
}

/// Solution curves of a fitted gene for plotting: returns (t, u, s) on a grid
/// of `num` points along induction and repression.
pub fn gene_trajectory(gd: &GeneDynamics, num: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let beta = gd.beta * gd.scaling;
    let (u0_, s0_) = solution(gd.t_, 0.0, 0.0, gd.alpha, beta, gd.gamma);
    let tmax = gd.t.iter().copied().fold(gd.t_, np_max);
    let ts = linspace(0.0, tmax.max(gd.t_ * 1.0001), num);
    let mut uu = Vec::with_capacity(num);
    let mut ss = Vec::with_capacity(num);
    for &t in &ts {
        let (u, s) = if t < gd.t_ { solution(t, 0.0, 0.0, gd.alpha, beta, gd.gamma) } else { solution(t - gd.t_, u0_, s0_, 0.0, beta, gd.gamma) };
        // back to unscaled unspliced units
        uu.push(u * gd.scaling);
        ss.push(s);
    }
    (ts, uu, ss)
}
