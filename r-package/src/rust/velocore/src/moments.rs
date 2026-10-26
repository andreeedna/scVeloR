//! First and second order moments over kNN neighbourhoods
//! (scvelo.pp.moments and second_order_moments).

use crate::matrix::{Csr, RowMat};
use rayon::prelude::*;

/// M = C · X for a sparse cells x genes matrix X, dense output.
pub fn first_order(conn: &Csr, x: &Csr) -> RowMat {
    conn.dot_sparse_dense_out(x)
}

/// M = C · X for a dense cells x genes matrix.
pub fn first_order_dense(conn: &Csr, x: &RowMat) -> RowMat {
    conn.dot_dense(x)
}

/// Element-wise product of two sparse matrices with the same shape.
pub fn hadamard(a: &Csr, b: &Csr) -> Csr {
    assert_eq!((a.nrow, a.ncol), (b.nrow, b.ncol));
    let rows: Vec<Vec<(usize, f64)>> = (0..a.nrow)
        .into_par_iter()
        .map(|i| {
            let (ai, av) = a.row(i);
            let (bi, bv) = b.row(i);
            let mut out = Vec::new();
            let (mut p, mut q) = (0, 0);
            while p < ai.len() && q < bi.len() {
                match ai[p].cmp(&bi[q]) {
                    std::cmp::Ordering::Less => p += 1,
                    std::cmp::Ordering::Greater => q += 1,
                    std::cmp::Ordering::Equal => {
                        let v = av[p] * bv[q];
                        if v != 0.0 {
                            out.push((ai[p], v));
                        }
                        p += 1;
                        q += 1;
                    }
                }
            }
            out
        })
        .collect();
    Csr::from_rows(a.ncol, rows)
}

/// Second order moments: (Mss, Mus) = (C·(s∘s), C·(s∘u)).
pub fn second_order(conn: &Csr, s: &Csr, u: &Csr) -> (RowMat, RowMat) {
    let ss = hadamard(s, s);
    let us = hadamard(s, u);
    (conn.dot_sparse_dense_out(&ss), conn.dot_sparse_dense_out(&us))
}

/// Dense variant of [`second_order`].
pub fn second_order_dense(conn: &Csr, s: &RowMat, u: &RowMat) -> (RowMat, RowMat) {
    let mut ss = s.clone();
    ss.data.par_iter_mut().for_each(|v| *v = *v * *v);
    let mut us = s.clone();
    us.data.par_iter_mut().zip(u.data.par_iter()).for_each(|(a, b)| *a *= *b);
    (conn.dot_dense(&ss), conn.dot_dense(&us))
}
