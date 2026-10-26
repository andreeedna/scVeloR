//! Nelder-Mead simplex, a line-by-line port of scipy.optimize's
//! `_minimize_neldermead` (non-adaptive, no bounds) so that the dynamical model
//! follows the same optimisation path as scVelo.

/// Minimise `f` from `x0` with at most `maxiter` iterations
/// (scipy counts the initial simplex as iteration 1).
pub fn nelder_mead(f: impl Fn(&[f64]) -> f64, x0: &[f64], maxiter: usize, xatol: f64, fatol: f64) -> Vec<f64> {
    let (rho, chi, psi, sigma) = (1.0, 2.0, 0.5, 0.5);
    let nonzdelt = 0.05;
    let zdelt = 0.00025;
    let n = x0.len();
    let mut sim: Vec<Vec<f64>> = Vec::with_capacity(n + 1);
    sim.push(x0.to_vec());
    for k in 0..n {
        let mut y = x0.to_vec();
        if y[k] != 0.0 {
            y[k] *= 1.0 + nonzdelt;
        } else {
            y[k] = zdelt;
        }
        sim.push(y);
    }
    let mut fsim: Vec<f64> = sim.iter().map(|x| f(x)).collect();
    sort_simplex(&mut sim, &mut fsim);

    let mut iterations = 1usize;
    while iterations < maxiter {
        let xdev = sim[1..]
            .iter()
            .flat_map(|x| x.iter().zip(&sim[0]).map(|(a, b)| (a - b).abs()))
            .fold(0.0f64, f64::max);
        let fdev = fsim[1..].iter().map(|v| (fsim[0] - v).abs()).fold(0.0f64, f64::max);
        let any_nan = fsim.iter().any(|v| v.is_nan()) || sim.iter().flatten().any(|v| v.is_nan());
        if !any_nan && xdev <= xatol && fdev <= fatol {
            break;
        }
        let mut xbar = vec![0.0; n];
        for x in &sim[..n] {
            for k in 0..n {
                xbar[k] += x[k];
            }
        }
        for v in xbar.iter_mut() {
            *v /= n as f64;
        }
        let worst = sim[n].clone();
        let xr: Vec<f64> = (0..n).map(|k| (1.0 + rho) * xbar[k] - rho * worst[k]).collect();
        let fxr = f(&xr);
        let mut doshrink = false;
        if fxr < fsim[0] {
            let xe: Vec<f64> = (0..n).map(|k| (1.0 + rho * chi) * xbar[k] - rho * chi * worst[k]).collect();
            let fxe = f(&xe);
            if fxe < fxr {
                sim[n] = xe;
                fsim[n] = fxe;
            } else {
                sim[n] = xr;
                fsim[n] = fxr;
            }
        } else if fxr < fsim[n - 1] {
            sim[n] = xr;
            fsim[n] = fxr;
        } else if fxr < fsim[n] {
            let xc: Vec<f64> = (0..n).map(|k| (1.0 + psi * rho) * xbar[k] - psi * rho * worst[k]).collect();
            let fxc = f(&xc);
            if fxc <= fxr {
                sim[n] = xc;
                fsim[n] = fxc;
            } else {
                doshrink = true;
            }
        } else {
            let xcc: Vec<f64> = (0..n).map(|k| (1.0 - psi) * xbar[k] + psi * worst[k]).collect();
            let fxcc = f(&xcc);
            if fxcc < fsim[n] {
                sim[n] = xcc;
                fsim[n] = fxcc;
            } else {
                doshrink = true;
            }
        }
        if doshrink {
            for j in 1..=n {
                let s0 = sim[0].clone();
                for k in 0..n {
                    sim[j][k] = s0[k] + sigma * (sim[j][k] - s0[k]);
                }
                fsim[j] = f(&sim[j]);
            }
        }
        iterations += 1;
        sort_simplex(&mut sim, &mut fsim);
    }
    sim.swap_remove(0)
}

/// Stable sort by function value (numpy argsort on <= 16 elements is an
/// insertion sort, hence stable). NaN sorts last like numpy.
fn sort_simplex(sim: &mut Vec<Vec<f64>>, fsim: &mut Vec<f64>) {
    let mut idx: Vec<usize> = (0..fsim.len()).collect();
    idx.sort_by(|&a, &b| {
        let (x, y) = (fsim[a], fsim[b]);
        match (x.is_nan(), y.is_nan()) {
            (true, true) => std::cmp::Ordering::Equal,
            (true, false) => std::cmp::Ordering::Greater,
            (false, true) => std::cmp::Ordering::Less,
            _ => x.partial_cmp(&y).unwrap(),
        }
    });
    let s: Vec<Vec<f64>> = idx.iter().map(|&i| sim[i].clone()).collect();
    let f: Vec<f64> = idx.iter().map(|&i| fsim[i]).collect();
    *sim = s;
    *fsim = f;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rosenbrock() {
        let f = |x: &[f64]| (1.0 - x[0]).powi(2) + 100.0 * (x[1] - x[0] * x[0]).powi(2);
        let x = nelder_mead(f, &[-1.2, 1.0], 2000, 1e-10, 1e-10);
        assert!((x[0] - 1.0).abs() < 1e-4 && (x[1] - 1.0).abs() < 1e-4);
    }
}
