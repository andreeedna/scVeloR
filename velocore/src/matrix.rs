//! Minimal dense and sparse matrix types.
//!
//! Layout conventions (chosen so R / numpy buffers can be used without copies):
//! * [`RowMat`] is row-major: row `i` is contiguous. A genes x cells matrix
//!   coming from R (column-major) is exactly a cells x genes `RowMat`.
//! * [`Csr`] is compressed sparse row. An R `dgCMatrix` of genes x cells is the
//!   CSR representation of the cells x genes matrix.

use rayon::prelude::*;

#[derive(Clone, Debug, PartialEq)]
pub struct RowMat {
    pub nrow: usize,
    pub ncol: usize,
    pub data: Vec<f64>,
}

impl RowMat {
    pub fn zeros(nrow: usize, ncol: usize) -> Self {
        RowMat { nrow, ncol, data: vec![0.0; nrow * ncol] }
    }

    pub fn from_vec(nrow: usize, ncol: usize, data: Vec<f64>) -> Self {
        assert_eq!(data.len(), nrow * ncol, "data length does not match shape");
        RowMat { nrow, ncol, data }
    }

    #[inline]
    pub fn row(&self, i: usize) -> &[f64] {
        &self.data[i * self.ncol..(i + 1) * self.ncol]
    }

    #[inline]
    pub fn row_mut(&mut self, i: usize) -> &mut [f64] {
        let n = self.ncol;
        &mut self.data[i * n..(i + 1) * n]
    }

    #[inline]
    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.data[i * self.ncol + j]
    }

    #[inline]
    pub fn set(&mut self, i: usize, j: usize, v: f64) {
        self.data[i * self.ncol + j] = v;
    }

    /// Copy of column `j`.
    pub fn col(&self, j: usize) -> Vec<f64> {
        (0..self.nrow).map(|i| self.data[i * self.ncol + j]).collect()
    }

    pub fn transpose(&self) -> RowMat {
        let mut out = RowMat::zeros(self.ncol, self.nrow);
        let (nr, nc) = (self.nrow, self.ncol);
        out.data.par_chunks_mut(nr).enumerate().for_each(|(j, orow)| {
            for i in 0..nr {
                orow[i] = self.data[i * nc + j];
            }
        });
        out
    }

    /// Keep only the given columns, in the given order.
    pub fn select_cols(&self, cols: &[usize]) -> RowMat {
        let mut out = RowMat::zeros(self.nrow, cols.len());
        let nc = cols.len();
        out.data.par_chunks_mut(nc.max(1)).enumerate().for_each(|(i, orow)| {
            if nc == 0 {
                return;
            }
            let r = self.row(i);
            for (k, &c) in cols.iter().enumerate() {
                orow[k] = r[c];
            }
        });
        out
    }

    pub fn select_rows(&self, rows: &[usize]) -> RowMat {
        let mut out = RowMat::zeros(rows.len(), self.ncol);
        for (k, &r) in rows.iter().enumerate() {
            out.row_mut(k).copy_from_slice(self.row(r));
        }
        out
    }
}

/// Compressed sparse row matrix.
#[derive(Clone, Debug, PartialEq)]
pub struct Csr {
    pub nrow: usize,
    pub ncol: usize,
    pub indptr: Vec<usize>,
    pub indices: Vec<usize>,
    pub data: Vec<f64>,
}

impl Csr {
    pub fn new(nrow: usize, ncol: usize, indptr: Vec<usize>, indices: Vec<usize>, data: Vec<f64>) -> Self {
        assert_eq!(indptr.len(), nrow + 1, "indptr must have nrow + 1 entries");
        assert_eq!(indices.len(), data.len());
        assert_eq!(*indptr.last().unwrap(), data.len());
        Csr { nrow, ncol, indptr, indices, data }
    }

    /// Build from per-row (column index, value) lists. Rows are sorted by column.
    pub fn from_rows(ncol: usize, rows: Vec<Vec<(usize, f64)>>) -> Self {
        let nrow = rows.len();
        let nnz: usize = rows.iter().map(|r| r.len()).sum();
        let mut indptr = Vec::with_capacity(nrow + 1);
        let mut indices = Vec::with_capacity(nnz);
        let mut data = Vec::with_capacity(nnz);
        indptr.push(0);
        for mut r in rows {
            r.sort_unstable_by_key(|e| e.0);
            for (c, v) in r {
                indices.push(c);
                data.push(v);
            }
            indptr.push(indices.len());
        }
        Csr { nrow, ncol, indptr, indices, data }
    }

    #[inline]
    pub fn row(&self, i: usize) -> (&[usize], &[f64]) {
        let (a, b) = (self.indptr[i], self.indptr[i + 1]);
        (&self.indices[a..b], &self.data[a..b])
    }

    pub fn nnz(&self) -> usize {
        self.data.len()
    }

    /// y = A x
    pub fn dot_vec(&self, x: &[f64]) -> Vec<f64> {
        assert_eq!(x.len(), self.ncol);
        (0..self.nrow)
            .into_par_iter()
            .map(|i| {
                let (idx, val) = self.row(i);
                idx.iter().zip(val).map(|(&j, &v)| v * x[j]).sum()
            })
            .collect()
    }

    /// Y = A X for a dense row-major X (ncol x k).
    pub fn dot_dense(&self, x: &RowMat) -> RowMat {
        assert_eq!(x.nrow, self.ncol);
        let k = x.ncol;
        let mut out = RowMat::zeros(self.nrow, k);
        if k == 0 {
            return out;
        }
        out.data.par_chunks_mut(k).enumerate().for_each(|(i, orow)| {
            let (idx, val) = self.row(i);
            for (&j, &v) in idx.iter().zip(val) {
                let xr = x.row(j);
                for c in 0..k {
                    orow[c] += v * xr[c];
                }
            }
        });
        out
    }

    /// Y = A X for a sparse X (ncol x k), densified output.
    pub fn dot_sparse_dense_out(&self, x: &Csr) -> RowMat {
        assert_eq!(x.nrow, self.ncol);
        let k = x.ncol;
        let mut out = RowMat::zeros(self.nrow, k);
        if k == 0 {
            return out;
        }
        out.data.par_chunks_mut(k).enumerate().for_each(|(i, orow)| {
            let (idx, val) = self.row(i);
            for (&j, &v) in idx.iter().zip(val) {
                let (xi, xv) = x.row(j);
                for (&c, &w) in xi.iter().zip(xv) {
                    orow[c] += v * w;
                }
            }
        });
        out
    }

    pub fn transpose(&self) -> Csr {
        let mut counts = vec![0usize; self.ncol + 1];
        for &j in &self.indices {
            counts[j + 1] += 1;
        }
        for j in 0..self.ncol {
            counts[j + 1] += counts[j];
        }
        let indptr = counts.clone();
        let mut next = counts;
        let mut indices = vec![0usize; self.nnz()];
        let mut data = vec![0.0; self.nnz()];
        for i in 0..self.nrow {
            let (idx, val) = self.row(i);
            for (&j, &v) in idx.iter().zip(val) {
                let p = next[j];
                indices[p] = i;
                data[p] = v;
                next[j] += 1;
            }
        }
        Csr { nrow: self.ncol, ncol: self.nrow, indptr, indices, data }
    }

    /// Divide each row by the sum of absolute values of its entries
    /// (rows summing to zero are left untouched).
    pub fn normalize_rows_abs(&mut self) {
        for i in 0..self.nrow {
            let (a, b) = (self.indptr[i], self.indptr[i + 1]);
            let s: f64 = self.data[a..b].iter().map(|v| v.abs()).sum();
            if s > 0.0 {
                for v in &mut self.data[a..b] {
                    *v /= s;
                }
            }
        }
    }

    /// Drop stored entries equal to zero (or NaN).
    pub fn eliminate_zeros(&mut self) {
        let mut indptr = Vec::with_capacity(self.nrow + 1);
        let mut indices = Vec::with_capacity(self.nnz());
        let mut data = Vec::with_capacity(self.nnz());
        indptr.push(0);
        for i in 0..self.nrow {
            let (idx, val) = self.row(i);
            for (&j, &v) in idx.iter().zip(val) {
                if v != 0.0 && !v.is_nan() {
                    indices.push(j);
                    data.push(v);
                }
            }
            indptr.push(indices.len());
        }
        self.indptr = indptr;
        self.indices = indices;
        self.data = data;
    }

    /// Row-wise maximum over stored entries, 0 for empty rows (scipy semantics
    /// for matrices with implicit zeros).
    pub fn row_max(&self) -> Vec<f64> {
        (0..self.nrow)
            .map(|i| {
                let (idx, val) = self.row(i);
                let mut m = f64::NEG_INFINITY;
                for &v in val {
                    if v > m {
                        m = v;
                    }
                }
                if idx.len() < self.ncol {
                    m = m.max(0.0);
                }
                if m == f64::NEG_INFINITY {
                    0.0
                } else {
                    m
                }
            })
            .collect()
    }

    /// Map every stored value.
    pub fn map(&self, f: impl Fn(f64) -> f64 + Sync) -> Csr {
        let mut out = self.clone();
        out.data.par_iter_mut().for_each(|v| *v = f(*v));
        out
    }

    /// Triplets (row, col, value) with zero-based indices.
    pub fn to_triplets(&self) -> (Vec<usize>, Vec<usize>, Vec<f64>) {
        let mut r = Vec::with_capacity(self.nnz());
        for i in 0..self.nrow {
            for _ in self.indptr[i]..self.indptr[i + 1] {
                r.push(i);
            }
        }
        (r, self.indices.clone(), self.data.clone())
    }

    pub fn from_triplets(nrow: usize, ncol: usize, rows: &[usize], cols: &[usize], vals: &[f64]) -> Csr {
        let mut per_row: Vec<Vec<(usize, f64)>> = vec![Vec::new(); nrow];
        for ((&i, &j), &v) in rows.iter().zip(cols).zip(vals) {
            per_row[i].push((j, v));
        }
        // merge duplicates by summation
        let rows: Vec<Vec<(usize, f64)>> = per_row
            .into_iter()
            .map(|mut r| {
                r.sort_unstable_by_key(|e| e.0);
                let mut out: Vec<(usize, f64)> = Vec::with_capacity(r.len());
                for (j, v) in r {
                    if let Some(last) = out.last_mut() {
                        if last.0 == j {
                            last.1 += v;
                            continue;
                        }
                    }
                    out.push((j, v));
                }
                out
            })
            .collect();
        Csr::from_rows(ncol, rows)
    }

    pub fn to_dense(&self) -> RowMat {
        let mut out = RowMat::zeros(self.nrow, self.ncol);
        for i in 0..self.nrow {
            let (idx, val) = self.row(i);
            for (&j, &v) in idx.iter().zip(val) {
                out.data[i * self.ncol + j] += v;
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn csr_roundtrip() {
        let a = Csr::from_triplets(2, 3, &[0, 1, 1], &[2, 0, 2], &[1.0, 2.0, 3.0]);
        let t = a.transpose();
        assert_eq!(t.to_dense().data, vec![0.0, 2.0, 0.0, 0.0, 1.0, 3.0]);
        assert_eq!(a.dot_vec(&[1.0, 1.0, 1.0]), vec![1.0, 5.0]);
    }
}
