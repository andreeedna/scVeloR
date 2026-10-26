# Low-level bindings to the Rust core (src/rust). Not exported; see R/core.R.
# nolint start

#' @useDynLib scVelo2, .registration = TRUE
#' @import Matrix
#' @importFrom methods as is new
NULL

rs_knn <- function(x, n, d, k, n_threads) .Call(wrap__rs_knn, x, n, d, k, n_threads)

rs_connectivities <- function(idx, dist, n, k, mode, recurse) .Call(wrap__rs_connectivities, idx, dist, n, k, mode, recurse)

rs_moments <- function(idx, dist, n, k, s_p, s_i, s_x, u_p, u_i, u_x, ngenes, second_order) .Call(wrap__rs_moments, idx, dist, n, k, s_p, s_i, s_x, u_p, u_i, u_x, ngenes, second_order)

rs_velocity_steady <- function(ms, mu, mss, mus, n, g, stochastic, fit_offset, fit_offset2, perc_lo, perc_hi, min_r2, min_ratio, r2_adjusted, highly_variable, n_threads) .Call(wrap__rs_velocity_steady, ms, mu, mss, mus, n, g, stochastic, fit_offset, fit_offset2, perc_lo, perc_hi, min_r2, min_ratio, r2_adjusted, highly_variable, n_threads)

rs_velocity_graph <- function(x, v, n, g, idx, dist, k, n_recurse, sqrt_transform, max_neighs, seed, n_threads) .Call(wrap__rs_velocity_graph, x, v, n, g, idx, dist, k, n_recurse, sqrt_transform, max_neighs, seed, n_threads)

rs_transition_matrix <- function(i, j, x, n, self_transition, scale, self_transitions, use_negative_cosines, backward) .Call(wrap__rs_transition_matrix, i, j, x, n, self_transition, scale, self_transitions, use_negative_cosines, backward)

rs_velocity_embedding <- function(i, j, x, n, self_transition, scale, x_emb, d, autoscale) .Call(wrap__rs_velocity_embedding, i, j, x, n, self_transition, scale, x_emb, d, autoscale)

rs_velocity_grid <- function(x_emb, v_emb, n, density, smooth, n_neighbors, min_mass, autoscale, stream_mode, cutoff_perc) .Call(wrap__rs_velocity_grid, x_emb, v_emb, n, density, smooth, n_neighbors, min_mass, autoscale, stream_mode, cutoff_perc)

rs_streamlines <- function(xs, ys, u, v, density_x, density_y, min_length, max_length, direction) .Call(wrap__rs_streamlines, xs, ys, u, v, density_x, density_y, min_length, max_length, direction)

rs_velocity_confidence <- function(v, n, g, idx, dist, k) .Call(wrap__rs_velocity_confidence, v, n, g, idx, dist, k)

rs_confidence_transition <- function(i, j, tx, n, x, v, g) .Call(wrap__rs_confidence_transition, i, j, tx, n, x, v, g)

rs_recover_dynamics <- function(mu, ms, n, g, idx, dist, k, max_iter, fit_scaling, fit_steady_states, fit_connected_states, projection, t_max, n_threads) .Call(wrap__rs_recover_dynamics, mu, ms, n, g, idx, dist, k, max_iter, fit_scaling, fit_steady_states, fit_connected_states, projection, t_max, n_threads)

rs_dynamical_velocity <- function(mu, ms, n, g, alpha, beta, gamma, t_, scaling, std_u, std_s, fit_tau, fit_tau_, idx, dist, k, use_connectivities) .Call(wrap__rs_dynamical_velocity, mu, ms, n, g, alpha, beta, gamma, t_, scaling, std_u, std_s, fit_tau, fit_tau_, idx, dist, k, use_connectivities)

rs_latent_time <- function(fit_t, n, g, roots, fates, idx, dist, k, min_confidence) .Call(wrap__rs_latent_time, fit_t, n, g, roots, fates, idx, dist, k, min_confidence)

rs_pseudotime_kernel <- function(i, j, x, n) .Call(wrap__rs_pseudotime_kernel, i, j, x, n)

rs_dpt <- function(evals, evecs, n, root) .Call(wrap__rs_dpt, evals, evecs, n, root)

rs_gene_trajectory <- function(alpha, beta, gamma, t_, scaling, tmax, num) .Call(wrap__rs_gene_trajectory, alpha, beta, gamma, t_, scaling, tmax, num)

# nolint end
