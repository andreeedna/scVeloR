sim <- local({
  v <- simulate_splicing(n_cells = 600, n_genes = 120, seed = 7)
  v <- velocity_filter_normalize(v, min_shared_counts = 10, n_top_genes = NULL)
  velocity_moments(v, n_neighbors = 30, n_pcs = 20)
})

test_that("kNN graph and moments have the expected shape", {
  expect_equal(dim(sim$knn$idx), c(600, 29))
  expect_false(any(sim$knn$idx == seq_len(600)))
  expect_equal(dim(sim$layers$Ms), c(length(sim$genes), 600))
  expect_true(all(sim$layers$Ms >= 0))
})

test_that("stochastic velocities point forward in time", {
  v <- velocity(sim, mode = "stochastic")
  expect_gt(sum(v$var$velocity_genes), 10)
  v <- velocity_graph(v)
  # cells should mostly transition towards later true times
  tm <- transition_matrix(v)
  t <- v$obs$true_time
  expected_t <- as.numeric(tm %*% t)
  expect_gt(mean(expected_t > t), 0.6)
})

test_that("dynamical model recovers time ordering", {
  v <- velocity(sim, mode = "dynamical")
  expect_true(any(!is.na(v$var$fit_alpha)))
  ok <- !is.na(v$var$fit_gamma)
  # fitted degradation/splicing ratios correlate with the truth
  truth <- v$var_true[v$genes[ok], ]
  expect_gt(stats::cor(v$var$fit_gamma[ok] / v$var$fit_beta[ok], truth$gamma / truth$beta, method = "spearman"), 0.5)
  v <- velocity_graph(v)
  v <- latent_time(v)
  expect_gt(stats::cor(v$obs$latent_time, v$obs$true_time, method = "spearman"), 0.8)
})

test_that("embedding projection, grid and streamlines run", {
  v <- velocity_graph(velocity(sim))
  v$reductions$umap <- v$reductions$pca[, 1:2]
  v <- velocity_embedding(v, "umap")
  expect_equal(dim(v$velocity_embeddings$umap), c(600, 2))
  p <- plot_velocity_stream(v, color_by = "true_time")
  expect_s3_class(p, "ggplot")
  expect_s3_class(plot_velocity_grid(v), "ggplot")
  expect_s3_class(plot_phase(v, v$genes[1:2]), "ggplot")
})

test_that("cell barcodes are matched between velocyto and Cell Ranger", {
  loom <- c("s1:AAACCTGAGAGGGATAx", "s1:AAACCTGAGCCTTGATx")
  obj <- c("AAACCTGAGCCTTGAT-1", "AAACCTGAGAGGGATA-1", "TTTTTTTTTTTTTTTT-1")
  expect_equal(match_cells(loom, obj), c(2L, 1L, NA))
})
