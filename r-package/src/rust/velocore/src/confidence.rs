//! Velocity coherence scores (scvelo.tl.velocity_confidence).

use crate::matrix::{Csr, RowMat};
use crate::neighbors::Knn;
use rayon::prelude::*;

fn centered(v: &[f64]) -> Vec<f64> {
    let m = v.iter().sum::<f64>() / v.len() as f64;
    v.iter().map(|x| x - m).collect()
}

fn norm(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

/// Returns (velocity_length, velocity_confidence): the norm of the centred
/// velocity and its mean correlation with the neighbours' velocities.
pub fn velocity_confidence(v: &RowMat, knn: &Knn) -> (Vec<f64>, Vec<f64>) {
    let n = v.nrow;
    let vc: Vec<Vec<f64>> = (0..n).into_par_iter().map(|i| centered(v.row(i))).collect();
    let norms: Vec<f64> = vc.iter().map(|r| norm(r)).collect();
    let conf: Vec<f64> = (0..n)
        .into_par_iter()
        .map(|i| {
            let nb = knn.neighbors(i);
            let mut s = 0.0;
            for &j in nb {
                let d: f64 = vc[j].iter().zip(&vc[i]).map(|(a, b)| a * b).sum();
                s += d / (norms[j] * norms[i]);
            }
            (s / nb.len() as f64).max(0.0)
        })
        .collect();
    let length = norms.iter().map(|v| (v * 100.0).round() / 100.0).collect();
    (length, conf)
}

/// Correlation between the velocity and the expected displacement under the
/// transition matrix `t` (velocity_confidence_transition).
pub fn confidence_transition(t: &Csr, x: &RowMat, v: &RowMat) -> Vec<f64> {
    let tx = t.dot_dense(x);
    (0..x.nrow)
        .into_par_iter()
        .map(|i| {
            let dx: Vec<f64> = tx.row(i).iter().zip(x.row(i)).map(|(a, b)| a - b).collect();
            let dx = centered(&dx);
            let vi = centered(v.row(i));
            let mut nn = norm(&dx) * norm(&vi);
            if nn == 0.0 {
                nn = 1.0;
            }
            dx.iter().zip(&vi).map(|(a, b)| a * b).sum::<f64>() / nn
        })
        .collect()
}
