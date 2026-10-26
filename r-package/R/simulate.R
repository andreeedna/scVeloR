#' Simulate spliced / unspliced counts from splicing kinetics
#'
#' Cells are sampled along a common time axis; each gene follows the
#' transcriptional induction / repression model with its own rates and
#' switching time, and counts are drawn from Poisson distributions. Useful to
#' test the workflow and to check that velocities and latent time recover the
#' known dynamics.
#'
#' @param n_cells,n_genes size of the data set.
#' @param t_max total time.
#' @param seed random seed.
#' @return a `velo` object with the true time in `obs$true_time`; ready for
#'   [velocity_filter_normalize()] with `n_top_genes = NULL`.
#' @export
simulate_splicing <- function(n_cells = 800, n_genes = 150, t_max = 20, seed = 1) {
  set.seed(seed)
  t <- sort(stats::runif(n_cells, 0, t_max))
  alpha <- stats::rlnorm(n_genes, log(15), 0.4)
  beta <- stats::rlnorm(n_genes, log(1), 0.3)
  gamma <- stats::rlnorm(n_genes, log(0.6), 0.3)
  switch <- stats::runif(n_genes, 0.3, 0.8) * t_max
  sol <- function(tau, u0, s0, a, b, g) {
    eu <- exp(-b * tau)
    es <- exp(-g * tau)
    u <- u0 * eu + a / b * (1 - eu)
    s <- s0 * es + a / g * (1 - es) + (a - u0 * b) / (g - b) * (es - eu)
    cbind(u, s)
  }
  u <- s <- matrix(0, n_genes, n_cells)
  for (j in seq_len(n_genes)) {
    on <- t < switch[j]
    x <- matrix(0, n_cells, 2)
    x[on, ] <- sol(t[on], 0, 0, alpha[j], beta[j], gamma[j])
    x0 <- sol(switch[j], 0, 0, alpha[j], beta[j], gamma[j])
    x[!on, ] <- sol(t[!on] - switch[j], x0[1], x0[2], 0, beta[j], gamma[j])
    u[j, ] <- stats::rpois(n_cells, x[, 1])
    s[j, ] <- stats::rpois(n_cells, x[, 2])
  }
  dn <- list(sprintf("gene%03d", seq_len(n_genes)), sprintf("cell%04d", seq_len(n_cells)))
  sp <- as_dgc(Matrix::Matrix(s, sparse = TRUE, dimnames = dn))
  up <- as_dgc(Matrix::Matrix(u, sparse = TRUE, dimnames = dn))
  v <- new_velo(sp, up)
  v$obs$true_time <- t
  v$var_true <- data.frame(alpha = alpha, beta = beta, gamma = gamma, switch = switch, row.names = dn[[1]])
  v
}
