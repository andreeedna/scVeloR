//! Trace the dynamical model of one gene: cargo run --release --example trace_gene -- <dir> <var index>
use velocore::dynamics::{fit_gene, DynamicsOptions};
use velocore::neighbors::{moment_connectivities, Knn};
use velocore::RowMat;

fn read(dir: &str, name: &str) -> Vec<f64> {
    let b = std::fs::read(format!("{dir}/{name}.bin")).unwrap();
    b.chunks_exact(8).map(|c| f64::from_le_bytes(c.try_into().unwrap())).collect()
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let (dir, j): (&str, usize) = (&a[1], a[2].parse().unwrap());
    let ms = read(dir, "Ms");
    let mu = read(dir, "Mu");
    let idx = read(dir, "knn_idx");
    let n = 3696;
    let g = ms.len() / n;
    let k = idx.len() / n;
    let knn = Knn::from_indices(n, k, idx.iter().map(|v| *v as usize).collect(), vec![1.0; n * k]);
    let conn = moment_connectivities(&knn);
    let ms = RowMat::from_vec(n, g, ms).col(j);
    let mu = RowMat::from_vec(n, g, mu).col(j);
    std::env::set_var("VELOCORE_TRACE", "1");
    let r = fit_gene(&mu, &ms, Some(&conn), &DynamicsOptions::default());
    eprintln!("result alpha {} beta {} gamma {} t_ {} scaling {} lik {}", r.alpha, r.beta, r.gamma, r.t_, r.scaling, r.likelihood);
}
