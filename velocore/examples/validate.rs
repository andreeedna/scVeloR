//! Compare velocore against arrays exported by validation/reference_scvelo.py.
//!
//! cargo run --release --example validate -- validation/out/pancreas

use std::collections::HashMap;
use std::time::Instant;
use velocore::dynamics::{dynamical_velocity, recover_dynamics, DynamicsOptions};
use velocore::embedding::{project, quiver_autoscale, velocity_on_grid, GridOptions};
use velocore::graph::{transition_matrix, velocity_graph, TransitionOptions, VelocityGraphOptions};
use velocore::neighbors::{knn, moment_connectivities, Knn};
use velocore::velocity::{deterministic, stochastic, SteadyStateOptions};
use velocore::{Csr, RowMat};

struct Data {
    dir: String,
    shapes: HashMap<String, Vec<usize>>,
}

impl Data {
    fn open(dir: &str) -> Data {
        let txt = std::fs::read_to_string(format!("{dir}/shapes.json")).expect("shapes.json");
        let mut shapes = HashMap::new();
        // tiny parser for {"name": [a, b], ...}
        let mut rest = txt.as_str();
        while let Some(q) = rest.find('"') {
            let r = &rest[q + 1..];
            let e = r.find('"').unwrap();
            let name = r[..e].to_string();
            let r2 = &r[e + 1..];
            let lb = r2.find('[').unwrap();
            let rb = r2.find(']').unwrap();
            let dims: Vec<usize> = r2[lb + 1..rb]
                .split(',')
                .filter(|s| !s.trim().is_empty())
                .map(|s| s.trim().parse().unwrap())
                .collect();
            shapes.insert(name, dims);
            rest = &r2[rb + 1..];
        }
        Data { dir: dir.to_string(), shapes }
    }
    fn has(&self, name: &str) -> bool {
        self.shapes.contains_key(name)
    }
    fn vec(&self, name: &str) -> Vec<f64> {
        let bytes = std::fs::read(format!("{}/{name}.bin", self.dir)).unwrap_or_else(|_| panic!("missing {name}"));
        bytes.chunks_exact(8).map(|c| f64::from_le_bytes(c.try_into().unwrap())).collect()
    }
    fn mat(&self, name: &str) -> RowMat {
        let s = &self.shapes[name];
        let (r, c) = (s[0], if s.len() > 1 { s[1] } else { 1 });
        RowMat::from_vec(r, c, self.vec(name))
    }
}

fn pearson(a: &[f64], b: &[f64]) -> f64 {
    let pairs: Vec<(f64, f64)> = a.iter().zip(b).filter(|(x, y)| x.is_finite() && y.is_finite()).map(|(x, y)| (*x, *y)).collect();
    let n = pairs.len() as f64;
    let ma = pairs.iter().map(|p| p.0).sum::<f64>() / n;
    let mb = pairs.iter().map(|p| p.1).sum::<f64>() / n;
    let (mut sab, mut saa, mut sbb) = (0.0, 0.0, 0.0);
    for (x, y) in pairs {
        sab += (x - ma) * (y - mb);
        saa += (x - ma) * (x - ma);
        sbb += (y - mb) * (y - mb);
    }
    sab / (saa * sbb).sqrt()
}

fn max_abs_diff(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).filter(|(x, y)| x.is_finite() && y.is_finite()).map(|(x, y)| (x - y).abs()).fold(0.0, f64::max)
}

fn rel_err(a: &[f64], b: &[f64]) -> f64 {
    let num: f64 = a.iter().zip(b).filter(|(x, y)| x.is_finite() && y.is_finite()).map(|(x, y)| (x - y).powi(2)).sum();
    let den: f64 = b.iter().filter(|y| y.is_finite()).map(|y| y * y).sum();
    (num / den).sqrt()
}

fn report(name: &str, ours: &[f64], theirs: &[f64]) {
    println!(
        "{name:<34} pearson={:.6}  rel_err={:.2e}  max_abs={:.2e}",
        pearson(ours, theirs),
        rel_err(ours, theirs),
        max_abs_diff(ours, theirs)
    );
}

fn csr_from_dense(m: &RowMat) -> Csr {
    let rows: Vec<Vec<(usize, f64)>> = (0..m.nrow)
        .map(|i| m.row(i).iter().enumerate().filter(|(_, v)| **v != 0.0).map(|(j, v)| (j, *v)).collect())
        .collect();
    Csr::from_rows(m.ncol, rows)
}

fn main() {
    let dir = std::env::args().nth(1).unwrap_or_else(|| "validation/out/pancreas".into());
    let d = Data::open(&dir);
    let s = d.mat("spliced");
    let u = d.mat("unspliced");
    let (n, g) = (s.nrow, s.ncol);
    println!("cells={n} genes={g}");

    // ---- neighbours: scanpy's graph (for parity) and our exact kNN (for info)
    let knn_idx = d.mat("knn_idx");
    let knn_dist = d.mat("knn_dist");
    let k = knn_idx.ncol;
    let sc_knn = Knn::from_indices(n, k, knn_idx.data.iter().map(|v| *v as usize).collect(), knn_dist.data.clone());
    let t0 = Instant::now();
    let ours_knn = knn(&d.mat("X_pca"), k);
    let mut overlap = 0usize;
    for i in 0..n {
        let a: std::collections::HashSet<usize> = sc_knn.neighbors(i).iter().copied().collect();
        overlap += ours_knn.neighbors(i).iter().filter(|j| a.contains(j)).count();
    }
    println!("exact kNN {:.2}s, overlap with scanpy (approximate) kNN: {:.4}", t0.elapsed().as_secs_f64(), overlap as f64 / (n * k) as f64);

    // ---- moments
    let t0 = Instant::now();
    let conn = moment_connectivities(&sc_knn);
    let s_sp = csr_from_dense(&s);
    let u_sp = csr_from_dense(&u);
    let ms = conn.dot_sparse_dense_out(&s_sp);
    let mu = conn.dot_sparse_dense_out(&u_sp);
    println!("moments: {:.3}s", t0.elapsed().as_secs_f64());
    report("Ms", &ms.data, &d.vec("Ms"));
    report("Mu", &mu.data, &d.vec("Mu"));

    // ---- stochastic velocity
    let t0 = Instant::now();
    let (mss, mus) = velocore::moments::second_order(&conn, &s_sp, &u_sp);
    let vr = stochastic(&ms, &mu, &mss, &mus, None, &SteadyStateOptions::default());
    println!("stochastic velocity: {:.3}s", t0.elapsed().as_secs_f64());
    report("velocity_gamma", &vr.gamma, &d.vec("velocity_gamma"));
    report("velocity_r2", &vr.r2, &d.vec("velocity_r2"));
    let vg_ref: Vec<bool> = d.vec("velocity_genes").iter().map(|v| *v > 0.5).collect();
    let agree = vr.velocity_genes.iter().zip(&vg_ref).filter(|(a, b)| a == b).count();
    println!(
        "velocity_genes: ours={} ref={} agree={}/{}",
        vr.velocity_genes.iter().filter(|b| **b).count(),
        vg_ref.iter().filter(|b| **b).count(),
        agree,
        g
    );
    report("velocity (all genes)", &vr.velocity.data, &d.vec("velocity"));

    // ---- velocity graph (stochastic -> sqrt transform)
    let vg_idx: Vec<usize> = (0..g).filter(|&j| vr.velocity_genes[j]).collect();
    let xs = ms.select_cols(&vg_idx);
    let vs = vr.velocity.select_cols(&vg_idx);
    let t0 = Instant::now();
    let vgraph = velocity_graph(&xs, &vs, &sc_knn, &VelocityGraphOptions { sqrt_transform: true, ..Default::default() });
    println!("velocity graph: {:.3}s  nnz={}", t0.elapsed().as_secs_f64(), vgraph.graph.nnz());
    {
        let mut refm: HashMap<(usize, usize), f64> = HashMap::new();
        for name in ["graph_pos", "graph_neg"] {
            let m = d.mat(name);
            for r in 0..m.nrow {
                refm.insert((m.get(r, 0) as usize, m.get(r, 1) as usize), m.get(r, 2));
            }
        }
        let (rows, cols, vals) = vgraph.graph.to_triplets();
        let mut a = Vec::new();
        let mut b = Vec::new();
        for ((r, c), v) in rows.iter().zip(&cols).zip(&vals) {
            a.push(*v);
            b.push(*refm.get(&(*r, *c)).unwrap_or(&0.0));
        }
        println!("graph entries: ours={} ref={}", vals.len(), refm.len());
        report("velocity_graph (shared support)", &a, &b);
        report("self_transition", &vgraph.self_transition, &d.vec("self_transition"));
    }

    // ---- embedding
    let x_umap = d.mat("X_umap");
    let t = transition_matrix(&vgraph, &TransitionOptions { use_negative_cosines: true, ..Default::default() });
    let mut v_emb = project(&t, &x_umap);
    let f = 3.0 * quiver_autoscale(&x_umap, &v_emb);
    v_emb.data.iter_mut().for_each(|v| *v /= f);
    report("velocity_umap", &v_emb.data, &d.vec("velocity_umap"));

    let gs = velocity_on_grid(&x_umap, &v_emb, &GridOptions { autoscale: false, adjust_for_stream: true, ..Default::default() });
    let su: Vec<f64> = (0..gs.points.nrow).map(|i| gs.velocity.get(i, 0)).collect();
    let sv: Vec<f64> = (0..gs.points.nrow).map(|i| gs.velocity.get(i, 1)).collect();
    report("stream grid u", &su, &d.vec("stream_u"));
    report("stream grid v", &sv, &d.vec("stream_v"));
    let nan_ours = su.iter().filter(|v| v.is_nan()).count();
    let nan_ref = d.vec("stream_u").iter().filter(|v| v.is_nan()).count();
    println!("stream masked cells: ours={nan_ours} ref={nan_ref}");
    let gq = velocity_on_grid(&x_umap, &v_emb, &GridOptions::default());
    println!("grid arrows: ours={} ref={}", gq.points.nrow, d.shapes["grid_X"][0]);
    if gq.points.nrow == d.shapes["grid_X"][0] {
        report("grid arrows V", &gq.velocity.data, &d.vec("grid_V"));
    }

    // ---- confidence
    let (len, conf) = velocore::confidence::velocity_confidence(&vs, &sc_knn);
    report("velocity_length", &len, &d.vec("velocity_length"));
    report("velocity_confidence", &conf, &d.vec("velocity_confidence"));
    let tconf = transition_matrix(&vgraph, &TransitionOptions::default());
    let ct = velocore::confidence::confidence_transition(&tconf, &xs, &vs);
    report("velocity_confidence_transition", &ct, &d.vec("velocity_confidence_transition"));

    // ---- dynamical model
    if d.has("fit_alpha") {
        let fit_alpha = d.vec("fit_alpha");
        let dyn_idx: Vec<usize> = (0..g).filter(|&j| !fit_alpha[j].is_nan()).collect();
        println!("dynamical genes (ref fitted): {}", dyn_idx.len());
        // scVelo keeps moments in float32; use its exact values to isolate
        // algorithmic differences from rounding (set VALIDATE_OWN_MOMENTS=1 to use ours)
        let own = std::env::var("VALIDATE_OWN_MOMENTS").is_ok();
        let (mu_d, ms_d) = if own {
            (mu.select_cols(&dyn_idx), ms.select_cols(&dyn_idx))
        } else {
            (d.mat("Mu").select_cols(&dyn_idx), d.mat("Ms").select_cols(&dyn_idx))
        };
        let t0 = Instant::now();
        let fit = recover_dynamics(&mu_d, &ms_d, Some(&conn), &DynamicsOptions::default(), Some(20.0));
        println!("recover_dynamics: {:.2}s for {} genes", t0.elapsed().as_secs_f64(), dyn_idx.len());
        let pick = |name: &str| -> Vec<f64> {
            let v = d.vec(name);
            dyn_idx.iter().map(|&j| v[j]).collect()
        };
        let gd = &fit.genes;
        report("fit_alpha", &gd.iter().map(|x| x.alpha).collect::<Vec<_>>(), &pick("fit_alpha"));
        report("fit_beta", &gd.iter().map(|x| x.beta).collect::<Vec<_>>(), &pick("fit_beta"));
        report("fit_gamma", &gd.iter().map(|x| x.gamma).collect::<Vec<_>>(), &pick("fit_gamma"));
        report("fit_t_", &gd.iter().map(|x| x.t_).collect::<Vec<_>>(), &pick("fit_t_"));
        report("fit_scaling", &gd.iter().map(|x| x.scaling).collect::<Vec<_>>(), &pick("fit_scaling"));
        report("fit_likelihood", &gd.iter().map(|x| x.likelihood).collect::<Vec<_>>(), &pick("fit_likelihood"));
        report("fit_pval_steady", &gd.iter().map(|x| x.pval_steady).collect::<Vec<_>>(), &pick("fit_pval_steady"));
        let exact = gd
            .iter()
            .zip(pick("fit_gamma"))
            .filter(|(a, b)| ((a.gamma - b) / b).abs() < 1e-6)
            .count();
        println!("genes with gamma identical to 1e-6: {exact}/{}", dyn_idx.len());
        let rg = pick("fit_gamma");
        let bad: Vec<usize> = (0..dyn_idx.len()).filter(|&k| ((gd[k].gamma - rg[k]) / rg[k]).abs() >= 1e-6).take(5).collect();
        for k in bad {
            println!("  mismatch gene col {} (var index {}): gamma ours={:.6} ref={:.6}", k, dyn_idx[k], gd[k].gamma, rg[k]);
        }
        let ref_t = d.mat("fit_t").select_cols(&dyn_idx);
        report("fit_t", &fit.fit_t.data, &ref_t.data);
        if let Some(dir64) = std::env::args().nth(2) {
            let d64 = Data::open(&dir64);
            let pick64 = |name: &str| -> Vec<f64> {
                let v = d64.vec(name);
                dyn_idx.iter().map(|&j| v[j]).collect()
            };
            println!("--- against scVelo run on float64 moments ({dir64})");
            for (name, ours) in [
                ("fit_alpha", gd.iter().map(|x| x.alpha).collect::<Vec<_>>()),
                ("fit_gamma", gd.iter().map(|x| x.gamma).collect::<Vec<_>>()),
                ("fit_t_", gd.iter().map(|x| x.t_).collect::<Vec<_>>()),
                ("fit_likelihood", gd.iter().map(|x| x.likelihood).collect::<Vec<_>>()),
            ] {
                report(&format!("{name} ours vs scVelo64"), &ours, &pick64(name));
                report(&format!("{name} scVelo32 vs scVelo64"), &pick(name), &pick64(name));
            }
            let g64 = pick64("fit_gamma");
            let same = gd.iter().zip(&g64).filter(|(a, b)| ((a.gamma - *b) / *b).abs() < 1e-6).count();
            let same32 = pick("fit_gamma").iter().zip(&g64).filter(|(a, b)| ((*a - *b) / *b).abs() < 1e-6).count();
            println!("identical gamma: ours vs scVelo64 {same}/{}, scVelo32 vs scVelo64 {same32}/{}", dyn_idx.len(), dyn_idx.len());
        }
        let (vdyn, _vdyn_u) = dynamical_velocity(&mu_d, &ms_d, &fit, Some(&conn));
        let ref_v = d.mat("velocity_dyn").select_cols(&dyn_idx);
        report("velocity (dynamical)", &vdyn.data, &ref_v.data);
    }
}
