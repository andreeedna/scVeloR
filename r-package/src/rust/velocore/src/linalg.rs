//! Tiny dense linear algebra for small matrices.

/// Cyclic Jacobi eigen decomposition of a small symmetric matrix.
/// Returns (eigenvalues, eigenvectors as columns: evecs[row][col]).
pub fn jacobi_eigen(m: &[Vec<f64>]) -> (Vec<f64>, Vec<Vec<f64>>) {
    let n = m.len();
    let mut a: Vec<Vec<f64>> = m.to_vec();
    let mut v = vec![vec![0.0; n]; n];
    for i in 0..n {
        v[i][i] = 1.0;
    }
    for _sweep in 0..100 {
        let off: f64 = (0..n)
            .flat_map(|p| (0..n).filter(move |&q| q != p).map(move |q| (p, q)))
            .map(|(p, q)| a[p][q] * a[p][q])
            .sum();
        if off < 1e-30 {
            break;
        }
        for p in 0..n {
            for q in (p + 1)..n {
                if a[p][q].abs() < 1e-300 {
                    continue;
                }
                let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let t = if theta == 0.0 { 1.0 } else { t };
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..n {
                    let akp = a[k][p];
                    let akq = a[k][q];
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for k in 0..n {
                    let apk = a[p][k];
                    let aqk = a[q][k];
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
                for k in 0..n {
                    let vkp = v[k][p];
                    let vkq = v[k][q];
                    v[k][p] = c * vkp - s * vkq;
                    v[k][q] = s * vkp + c * vkq;
                }
            }
        }
    }
    ((0..n).map(|i| a[i][i]).collect(), v)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn eig2() {
        let (e, _) = jacobi_eigen(&[vec![2.0, 1.0], vec![1.0, 2.0]]);
        let mut e = e;
        e.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!((e[0] - 1.0).abs() < 1e-12 && (e[1] - 3.0).abs() < 1e-12);
    }
}
