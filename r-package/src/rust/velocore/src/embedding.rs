//! Projection of velocities into a low-dimensional embedding and the grid
//! smoothing used by grid / stream plots (scvelo.tl.velocity_embedding,
//! scvelo.pl.velocity_embedding_grid.compute_velocity_on_grid).

use crate::matrix::{Csr, RowMat};
use crate::neighbors::knn_query;
use crate::stats::percentile;
use rayon::prelude::*;

/// Expected displacement under the transition matrix:
/// ṽ_i = Σ_j π_ij δ̂_ij − mean(π_i·) Σ_j δ̂_ij with unit displacement vectors δ̂.
/// `t` must not be renormalised after removing the diagonal (scVelo semantics).
pub fn project(t: &Csr, x_emb: &RowMat) -> RowMat {
    let (n, d) = (x_emb.nrow, x_emb.ncol);
    let mut out = RowMat::zeros(n, d);
    out.data.par_chunks_mut(d).enumerate().for_each(|(i, vrow)| {
        let (idx, val) = t.row(i);
        let xi = x_emb.row(i);
        let mut sum_p = 0.0;
        let mut cnt = 0usize;
        let mut sum_dx = vec![0.0; d];
        let mut acc = vec![0.0; d];
        for (&j, &p) in idx.iter().zip(val) {
            if j == i || p == 0.0 {
                continue;
            }
            let xj = x_emb.row(j);
            let mut nrm = 0.0;
            for k in 0..d {
                nrm += (xj[k] - xi[k]) * (xj[k] - xi[k]);
            }
            let nrm = nrm.sqrt();
            for k in 0..d {
                let mut dx = (xj[k] - xi[k]) / nrm;
                if dx.is_nan() {
                    dx = 0.0;
                }
                acc[k] += p * dx;
                sum_dx[k] += dx;
            }
            sum_p += p;
            cnt += 1;
        }
        if cnt > 0 {
            let mp = sum_p / cnt as f64;
            for k in 0..d {
                vrow[k] = acc[k] - mp * sum_dx[k];
            }
        }
    });
    out
}

/// Reproduces matplotlib's quiver auto-scaling as used by scVelo
/// (`quiver_autoscale`): returns the factor such that dividing the
/// embedded velocities by `3 * factor` gives scVelo's arrow lengths.
/// Uses matplotlib's default figure geometry (6.4 x 4.8 in, axes 0.775 x 0.77).
pub fn quiver_autoscale(x_emb: &RowMat, v_emb: &RowMat) -> f64 {
    let n = x_emb.nrow;
    let sf = x_emb.data.iter().fold(0.0f64, |a, b| a.max(b.abs()));
    let (mut xmin, mut xmax, mut ymin, mut ymax) = (f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY);
    for i in 0..n {
        let (x, y) = (x_emb.get(i, 0) / sf, x_emb.get(i, 1) / sf);
        if x.is_finite() && y.is_finite() {
            xmin = xmin.min(x);
            xmax = xmax.max(x);
            ymin = ymin.min(y);
            ymax = ymax.max(y);
        }
    }
    // default 5% margins on the view limits
    let vlw = (xmax - xmin) * 1.1;
    let vlh = (ymax - ymin) * 1.1;
    let dpi = 100.0;
    let bbw = 0.775 * 6.4 * dpi;
    let bbh = 0.77 * 4.8 * dpi;
    let sx = bbw / if vlw > 0.0 { vlw } else { 1.0 };
    let sy = bbh / if vlh > 0.0 { vlh } else { 1.0 };
    let mut s = 0.0;
    let mut cnt = 0usize;
    for i in 0..n {
        let (u, v) = (v_emb.get(i, 0), v_emb.get(i, 1));
        let l = ((u * sx).powi(2) + (v * sy).powi(2)).sqrt();
        if l.is_finite() {
            s += l;
            cnt += 1;
        }
    }
    let amean = if cnt > 0 { s / cnt as f64 } else { 0.0 };
    let sn = (n as f64).sqrt().max(10.0);
    let scale = 1.8 * amean * sn / bbw;
    scale / sf
}

#[derive(Clone, Debug)]
pub struct GridOptions {
    pub density: f64,
    pub smooth: f64,
    pub n_neighbors: Option<usize>,
    pub min_mass: Option<f64>,
    pub autoscale: bool,
    pub adjust_for_stream: bool,
    pub cutoff_perc: f64,
}

impl Default for GridOptions {
    fn default() -> Self {
        GridOptions { density: 1.0, smooth: 0.5, n_neighbors: None, min_mass: None, autoscale: true, adjust_for_stream: false, cutoff_perc: 5.0 }
    }
}

#[derive(Clone, Debug)]
pub struct GridVelocity {
    /// Grid coordinates along x and y (for stream mode) or the kept grid points.
    pub xs: Vec<f64>,
    pub ys: Vec<f64>,
    /// Points (x, y) and their velocities (u, v). In stream mode these are the
    /// full ny x nx grid in row-major order (y outer), with NaN for masked points.
    pub points: RowMat,
    pub velocity: RowMat,
}

/// Gaussian-kernel average of cell velocities on a regular grid.
pub fn velocity_on_grid(x_emb: &RowMat, v_emb: &RowMat, opts: &GridOptions) -> GridVelocity {
    // keep valid cells
    let valid: Vec<usize> = (0..x_emb.nrow)
        .filter(|&i| {
            (x_emb.get(i, 0) + x_emb.get(i, 1) + v_emb.get(i, 0) + v_emb.get(i, 1)).is_finite()
        })
        .collect();
    let mut x = RowMat::zeros(valid.len(), 2);
    let mut v = RowMat::zeros(valid.len(), 2);
    for (k, &i) in valid.iter().enumerate() {
        x.set(k, 0, x_emb.get(i, 0));
        x.set(k, 1, x_emb.get(i, 1));
        v.set(k, 0, v_emb.get(i, 0));
        v.set(k, 1, v_emb.get(i, 1));
    }
    let n_obs = x.nrow;
    let npts = (50.0 * opts.density) as usize;
    let mut grs: Vec<Vec<f64>> = Vec::new();
    for d in 0..2 {
        let col = x.col(d);
        let m = col.iter().copied().fold(f64::INFINITY, f64::min);
        let mm = col.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let m2 = m - 0.01 * (mm - m).abs();
        let mm2 = mm + 0.01 * (mm - m2).abs();
        grs.push(crate::stats::linspace(m2, mm2, npts));
    }
    // meshgrid: x varies fastest
    let ng = npts * npts;
    let mut grid = RowMat::zeros(ng, 2);
    for yi in 0..npts {
        for xi in 0..npts {
            grid.set(yi * npts + xi, 0, grs[0][xi]);
            grid.set(yi * npts + xi, 1, grs[1][yi]);
        }
    }
    let k = opts.n_neighbors.unwrap_or((n_obs as f64 / 50.0) as usize).max(1);
    let (nidx, ndst) = knn_query(&x, &grid, k);
    let k = nidx.len() / ng;
    let scale = ((grs[0][1] - grs[0][0]) + (grs[1][1] - grs[1][0])) / 2.0 * opts.smooth;
    let pdf = |d: f64| (-(d * d) / (2.0 * scale * scale)).exp() / (scale * (2.0 * std::f64::consts::PI).sqrt());
    let mut vgrid = RowMat::zeros(ng, 2);
    let mut p_mass = vec![0.0; ng];
    let mut length = vec![0.0; ng];
    for gi in 0..ng {
        let (mut su, mut sv, mut w) = (0.0, 0.0, 0.0);
        let (mut au, mut av) = (0.0, 0.0);
        for q in 0..k {
            let j = nidx[gi * k + q];
            let wt = pdf(ndst[gi * k + q]);
            su += v.get(j, 0) * wt;
            sv += v.get(j, 1) * wt;
            w += wt;
            au += v.get(j, 0).abs();
            av += v.get(j, 1).abs();
        }
        p_mass[gi] = w;
        let denom = w.max(1.0);
        vgrid.set(gi, 0, su / denom);
        vgrid.set(gi, 1, sv / denom);
        length[gi] = au / k as f64 + av / k as f64;
    }
    let min_mass = opts.min_mass.unwrap_or(1.0);
    if opts.adjust_for_stream {
        let mass: Vec<f64> = (0..ng).map(|i| (vgrid.get(i, 0).powi(2) + vgrid.get(i, 1).powi(2)).sqrt()).collect();
        let mmax = mass.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let thr = 10f64.powf(min_mass - 6.0).min(mmax * 0.9);
        let lthr = percentile(&length, opts.cutoff_perc);
        // grid points are ordered y-outer / x-inner, i.e. the [y, x] layout
        // matplotlib's streamplot expects; only u is masked, as in scVelo
        let mut masked = vec![false; ng];
        for gi in 0..ng {
            if mass[gi] < thr || length[gi] < lthr {
                masked[gi] = true;
            }
        }
        for gi in 0..ng {
            if masked[gi] {
                vgrid.set(gi, 0, f64::NAN);
            }
        }
        GridVelocity { xs: grs[0].clone(), ys: grs[1].clone(), points: grid, velocity: vgrid }
    } else {
        let pmax = percentile(&p_mass, 99.0) / 100.0 * min_mass;
        let keep: Vec<usize> = (0..ng).filter(|&i| p_mass[i] > pmax).collect();
        let pts = grid.select_rows(&keep);
        let mut vel = vgrid.select_rows(&keep);
        if opts.autoscale && !keep.is_empty() {
            let f = 3.0 * quiver_autoscale(&pts, &vel);
            if f > 0.0 {
                vel.data.iter_mut().for_each(|e| *e /= f);
            }
        }
        GridVelocity { xs: grs[0].clone(), ys: grs[1].clone(), points: pts, velocity: vel }
    }
}
