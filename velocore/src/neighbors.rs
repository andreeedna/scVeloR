//! Exact k-nearest-neighbour search and the neighbourhood graphs used by scVelo.

use crate::matrix::{Csr, RowMat};
use rayon::prelude::*;
use std::cmp::Ordering;
use std::collections::BinaryHeap;

/// kNN result. `indices[i]` / `distances[i]` hold the `k` nearest neighbours of
/// cell `i` **excluding the cell itself**, sorted by increasing distance.
/// This matches the `distances` graph produced by scanpy with
/// `n_neighbors = k + 1` (scanpy counts the cell itself as a neighbour).
#[derive(Clone, Debug)]
pub struct Knn {
    pub n: usize,
    pub k: usize,
    pub indices: Vec<usize>,
    pub distances: Vec<f64>,
}

impl Knn {
    #[inline]
    pub fn neighbors(&self, i: usize) -> &[usize] {
        &self.indices[i * self.k..(i + 1) * self.k]
    }

    #[inline]
    pub fn dists(&self, i: usize) -> &[f64] {
        &self.distances[i * self.k..(i + 1) * self.k]
    }

    /// Build from user-supplied neighbour indices (e.g. a Seurat graph).
    pub fn from_indices(n: usize, k: usize, indices: Vec<usize>, distances: Vec<f64>) -> Knn {
        assert_eq!(indices.len(), n * k);
        assert_eq!(distances.len(), n * k);
        Knn { n, k, indices, distances }
    }
}

#[derive(PartialEq)]
struct Cand(f64, usize);
impl Eq for Cand {}
impl PartialOrd for Cand {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Cand {
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.partial_cmp(&other.0).unwrap_or(Ordering::Equal).then(self.1.cmp(&other.1))
    }
}

/// Exact Euclidean kNN on the rows of `x` (cells x dims), excluding self.
pub fn knn(x: &RowMat, k: usize) -> Knn {
    let n = x.nrow;
    let k = k.min(n.saturating_sub(1));
    let norms: Vec<f64> = (0..n).map(|i| x.row(i).iter().map(|v| v * v).sum()).collect();
    let rows: Vec<(Vec<usize>, Vec<f64>)> = (0..n)
        .into_par_iter()
        .map(|i| {
            let xi = x.row(i);
            let mut heap: BinaryHeap<Cand> = BinaryHeap::with_capacity(k + 1);
            for j in 0..n {
                if j == i {
                    continue;
                }
                let xj = x.row(j);
                let mut dot = 0.0;
                for d in 0..xi.len() {
                    dot += xi[d] * xj[d];
                }
                let d2 = (norms[i] + norms[j] - 2.0 * dot).max(0.0);
                if heap.len() < k {
                    heap.push(Cand(d2, j));
                } else if let Some(top) = heap.peek() {
                    if d2 < top.0 {
                        heap.pop();
                        heap.push(Cand(d2, j));
                    }
                }
            }
            let mut v = heap.into_sorted_vec();
            v.truncate(k);
            // recompute exact distances to avoid cancellation error
            let idx: Vec<usize> = v.iter().map(|c| c.1).collect();
            let dist: Vec<f64> = idx
                .iter()
                .map(|&j| {
                    let xj = x.row(j);
                    xi.iter().zip(xj).map(|(a, b)| (a - b) * (a - b)).sum::<f64>().sqrt()
                })
                .collect();
            (idx, dist)
        })
        .collect();
    let mut indices = Vec::with_capacity(n * k);
    let mut distances = Vec::with_capacity(n * k);
    for (i, d) in rows {
        indices.extend(i);
        distances.extend(d);
    }
    Knn { n, k, indices, distances }
}

/// Neighbours of 2-d query points among reference points (used for grids).
/// Returns per query: (indices, distances) of the `k` nearest reference points.
pub fn knn_query(reference: &RowMat, queries: &RowMat, k: usize) -> (Vec<usize>, Vec<f64>) {
    let n = reference.nrow;
    let k = k.min(n).max(1);
    let res: Vec<(Vec<usize>, Vec<f64>)> = (0..queries.nrow)
        .into_par_iter()
        .map(|q| {
            let xq = queries.row(q);
            let mut heap: BinaryHeap<Cand> = BinaryHeap::with_capacity(k + 1);
            for j in 0..n {
                let xj = reference.row(j);
                let d2: f64 = xq.iter().zip(xj).map(|(a, b)| (a - b) * (a - b)).sum();
                if heap.len() < k {
                    heap.push(Cand(d2, j));
                } else if d2 < heap.peek().unwrap().0 {
                    heap.pop();
                    heap.push(Cand(d2, j));
                }
            }
            let v = heap.into_sorted_vec();
            (v.iter().map(|c| c.1).collect(), v.iter().map(|c| c.0.sqrt()).collect())
        })
        .collect();
    let mut idx = Vec::with_capacity(queries.nrow * k);
    let mut dst = Vec::with_capacity(queries.nrow * k);
    for (i, d) in res {
        idx.extend(i);
        dst.extend(d);
    }
    (idx, dst)
}

/// Row-normalised binary connectivities used for first/second order moments.
///
/// scVelo takes the (symmetric) UMAP connectivity graph, binarises it, adds the
/// diagonal and row-normalises. The support of the UMAP fuzzy union is the union
/// of the directed kNN graph and its transpose, so the result is the mean over
/// `{i} ∪ kNN(i) ∪ {j : i ∈ kNN(j)}`.
pub fn moment_connectivities(knn: &Knn) -> Csr {
    let n = knn.n;
    let mut adj: Vec<Vec<usize>> = (0..n).map(|i| {
        let mut v = Vec::with_capacity(knn.k * 2 + 1);
        v.push(i);
        v.extend_from_slice(knn.neighbors(i));
        v
    }).collect();
    for i in 0..n {
        for &j in knn.neighbors(i) {
            adj[j].push(i);
        }
    }
    let rows: Vec<Vec<(usize, f64)>> = adj
        .into_par_iter()
        .map(|mut v| {
            v.sort_unstable();
            v.dedup();
            let w = 1.0 / v.len() as f64;
            v.into_iter().map(|j| (j, w)).collect()
        })
        .collect();
    Csr::from_rows(n, rows)
}

/// Row-normalised binary connectivities on the directed kNN graph plus the
/// diagonal (scVelo `get_connectivities(adata, mode="distances")`).
/// With `recurse`, neighbours of neighbours are added with weight 0.5 before
/// clipping and normalising (`recurse_neighbors=True`).
pub fn knn_connectivities(knn: &Knn, recurse: bool) -> Csr {
    let n = knn.n;
    let rows: Vec<Vec<(usize, f64)>> = (0..n)
        .into_par_iter()
        .map(|i| {
            let mut v: Vec<(usize, f64)> = Vec::with_capacity(knn.k + 1);
            v.push((i, 1.0));
            for &j in knn.neighbors(i) {
                v.push((j, 1.0));
            }
            if recurse {
                // C + C·(0.5 C), clipped to [0, 1]
                let mut acc: std::collections::HashMap<usize, f64> = v.iter().copied().collect();
                let first: Vec<usize> = v.iter().map(|e| e.0).collect();
                for &j in &first {
                    acc.entry(j).and_modify(|w| *w += 0.5).or_insert(0.5);
                    for &l in knn.neighbors(j) {
                        *acc.entry(l).or_insert(0.0) += 0.5;
                    }
                }
                v = acc.into_iter().map(|(j, w)| (j, w.min(1.0))).collect();
            }
            let s: f64 = v.iter().map(|e| e.1).sum();
            v.into_iter().map(|(j, w)| (j, w / s)).collect()
        })
        .collect();
    Csr::from_rows(n, rows)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn knn_line() {
        let x = RowMat::from_vec(5, 1, vec![0.0, 1.0, 2.0, 3.0, 10.0]);
        let k = knn(&x, 2);
        assert_eq!(k.neighbors(0), &[1, 2]);
        assert_eq!(k.neighbors(4), &[3, 2]);
        let c = moment_connectivities(&k);
        // cell 4 is neighbour of nobody else except via its own list
        let (idx, w) = c.row(4);
        assert_eq!(idx, &[2, 3, 4]);
        assert!((w.iter().sum::<f64>() - 1.0).abs() < 1e-12);
    }
}
