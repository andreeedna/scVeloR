//! velocore: the numerical core of scVeloR.
//!
//! A from-scratch Rust implementation of the RNA velocity workflow popularised
//! by scVelo (Bergen et al., 2020): kNN moments, steady-state / stochastic /
//! dynamical models, velocity graph, transition matrices, embedding projection,
//! grid smoothing, streamlines and latent time. It has no I/O and no
//! dependency on a host language; R and Python bindings live in separate crates.

pub mod confidence;
pub mod dynamics;
pub mod embedding;
pub mod graph;
pub mod latent;
pub mod linalg;
pub mod matrix;
pub mod moments;
pub mod neighbors;
pub mod optim;
pub mod stats;
pub mod stream;
pub mod velocity;

pub use matrix::{Csr, RowMat};
