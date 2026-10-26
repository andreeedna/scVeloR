//! Gene-shared latent time (scvelo.tl.latent_time) from per-gene fitted times.

use crate::matrix::{Csr, RowMat};
use crate::stats::{percentile, scale01};
use rayon::prelude::*;

/// root_time: re-root gene-wise times at a given cell (cells x genes).
pub fn root_time(t: &RowMat, root: usize) -> RowMat {
    let (n, g) = (t.nrow, t.ncol);
    let t_root: Vec<f64> = t.row(root).to_vec();
    let mut t_origin = vec![f64::NEG_INFINITY; g];
    for i in 0..n {
        let r = t.row(i);
        for j in 0..g {
            let after = if r[j] >= t_root[j] { r[j] - t_root[j] } else { 0.0 };
            if after > t_origin[j] {
                t_origin[j] = after;
            }
        }
    }
    let mut out = RowMat::zeros(n, g);
    out.data.par_chunks_mut(g.max(1)).enumerate().for_each(|(i, orow)| {
        let r = t.row(i);
        for j in 0..g {
            orow[j] = if r[j] >= t_root[j] { r[j] - t_root[j] } else { r[j] + t_origin[j] };
        }
    });
    out
}

/// compute_shared_time: per-cell quantiles across genes, choosing the two
/// quantiles whose distribution is closest to uniform.
pub fn shared_time(t: &RowMat) -> Vec<f64> {
    let (n, g) = (t.nrow, t.ncol);
    let tmin = t.data.iter().copied().fold(f64::INFINITY, f64::min);
    let percs = [15.0, 20.0, 25.0, 30.0, 35.0];
    let rows: Vec<[f64; 5]> = (0..n)
        .into_par_iter()
        .map(|i| {
            let mut r: Vec<f64> = t.row(i).iter().map(|v| v - tmin).collect();
            r.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            let mut out = [0.0; 5];
            for (k, &p) in percs.iter().enumerate() {
                out[k] = if g > 0 { crate::stats::percentile_sorted(&r, p) } else { f64::NAN };
            }
            out
        })
        .collect();
    let mut tx: Vec<Vec<f64>> = (0..5).map(|k| rows.iter().map(|r| r[k]).collect()).collect();
    let linx = crate::stats::linspace(0.0, 1.0, n);
    let mut mse = vec![0.0; 5];
    for k in 0..5 {
        let mut mx = tx[k].iter().copied().fold(f64::NEG_INFINITY, f64::max);
        if mx == 0.0 {
            mx = 1.0;
        }
        tx[k].iter_mut().for_each(|v| *v /= mx);
        let mut sorted = tx[k].clone();
        sorted.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        mse[k] = sorted.iter().zip(&linx).map(|(a, b)| (a - b) * (a - b)).sum();
    }
    let mut order: Vec<usize> = (0..5).collect();
    order.sort_by(|&a, &b| mse[a].partial_cmp(&mse[b]).unwrap_or(std::cmp::Ordering::Equal));
    let mut ts: Vec<f64> = (0..n).map(|i| tx[order[0]][i] + tx[order[1]][i]).collect();
    let mx = ts.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    ts.iter_mut().for_each(|v| *v /= mx);
    ts
}

/// Latent time from fitted gene times `t` (cells x genes, genes already
/// filtered by likelihood), candidate root cells (at most 4 are used, in the
/// given order) and optional end points; `conn` are the moment connectivities.
pub fn latent_time(t: &RowMat, roots: &[usize], fates: &[usize], conn: &Csr, min_confidence: f64) -> Vec<f64> {
    let n = t.nrow;
    // scVelo re-roots the already re-rooted matrix in successive iterations;
    // that behaviour is kept for parity.
    let mut tcur = t.clone();
    let mut acc = vec![0.0; n];
    let roots: Vec<usize> = roots.iter().copied().take(4).collect();
    for &r in &roots {
        tcur = root_time(&tcur, r);
        let st = shared_time(&tcur);
        for i in 0..n {
            acc[i] += st[i];
        }
    }
    let k = roots.len().max(1) as f64;
    let mut lt = scale01(&acc.iter().map(|v| v / k).collect::<Vec<_>>());
    if !fates.is_empty() {
        let fates: Vec<usize> = fates.iter().copied().take(4).collect();
        let mut acc2 = vec![0.0; n];
        for &f in &fates {
            tcur = root_time(&tcur, f);
            let st = shared_time(&tcur);
            for i in 0..n {
                acc2[i] += 1.0 - st[i];
            }
        }
        let m: Vec<f64> = acc2.iter().map(|v| v / fates.len() as f64).collect();
        let sm = scale01(&m);
        lt = scale01(&lt.iter().zip(&sm).map(|(a, b)| a + 0.2 * b).collect::<Vec<_>>());
    }
    // local coherence
    let tc = conn.dot_vec(&lt);
    let z = lt.iter().zip(&tc).map(|(a, b)| a * b).sum::<f64>() / tc.iter().map(|b| b * b).sum::<f64>();
    let tmax = lt.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let low: Vec<bool> = (0..n)
        .map(|i| {
            let conf = (1.0 - (lt[i] / tmax - tc[i] * z / tmax).abs()).powi(2);
            conf < min_confidence
        })
        .collect();
    let smoothed: Vec<f64> = (0..n)
        .map(|i| {
            let (idx, val) = conn.row(i);
            idx.iter().zip(val).filter(|(&j, _)| !low[j]).map(|(&j, &w)| w * lt[j]).sum()
        })
        .collect();
    scale01(&smoothed)
}

/// Velocity pseudotime from a diffusion map (eigen decomposition of the
/// symmetric velocity kernel): scanpy's DPT distance from `root`.
pub fn dpt_from_root(evals: &[f64], evecs: &RowMat, root: usize) -> Vec<f64> {
    let n = evecs.nrow;
    let mut row = vec![0.0; n];
    for (j, &l) in evals.iter().enumerate() {
        let w = if l < 0.9994 { l / (1.0 - l) } else { 1.0 };
        let er = evecs.get(root, j);
        for i in 0..n {
            let d = w * (er - evecs.get(i, j));
            row[i] += d * d;
        }
    }
    let mut row: Vec<f64> = row.iter().map(|v| v.sqrt()).collect();
    let mx = row.iter().copied().filter(|v| v.is_finite()).fold(f64::NEG_INFINITY, f64::max);
    row.iter_mut().for_each(|v| *v /= mx);
    let _ = percentile; // keep import used in all cfgs
    row
}
