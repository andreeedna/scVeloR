# Minimal velocyto loom: layers are genes x cells in HDF5 (cells x genes in R).
write_test_loom <- function(path, s, u, cells, genes) {
  rhdf5::h5createFile(path)
  for (g in c("layers", "row_attrs", "col_attrs")) rhdf5::h5createGroup(path, g)
  rhdf5::h5write(t(as.matrix(s)), path, "/matrix")
  rhdf5::h5write(t(as.matrix(s)), path, "/layers/spliced")
  rhdf5::h5write(t(as.matrix(u)), path, "/layers/unspliced")
  rhdf5::h5write(genes, path, "/row_attrs/Gene")
  rhdf5::h5write(cells, path, "/col_attrs/CellID")
  rhdf5::h5closeAll()
}

test_that("loom files are read in blocks with the right orientation", {
  v <- simulate_splicing(n_cells = 30, n_genes = 7, seed = 2)
  path <- tempfile(fileext = ".loom")
  write_test_loom(path, v$spliced, v$unspliced, v$cells, v$genes)
  d <- read_loom(path, block_size = 2)
  expect_equal(as.matrix(d$spliced), as.matrix(v$spliced), ignore_attr = TRUE)
  expect_equal(as.matrix(d$unspliced), as.matrix(v$unspliced), ignore_attr = TRUE)
  expect_equal(dimnames(d$spliced), list(v$genes, v$cells))
  expect_null(d$ambiguous)
  # the same file twice needs sample names
  expect_error(read_loom(c(path, path)), "sample_names")
  expect_equal(ncol(read_loom(c(path, path), sample_names = c("a", "b"))$spliced), 60)
})

test_that("cell matching stays fast and exact on large objects", {
  set.seed(1)
  bc <- unique(replicate(44000, paste(sample(c("A", "C", "G", "T"), 16, TRUE), collapse = "")))[1:40000]
  loom <- paste0("s1:", bc, "x")
  obj <- paste0(bc[c(40000:2, 1)], "-1")
  t <- system.time(m <- match_cells(loom, obj))[["elapsed"]]
  expect_equal(m, c(40000:2, 1))
  expect_lt(t, 5)
})

test_that("analysis steps leave the user's random state alone", {
  v <- suppressMessages(velocity_graph(velocity(velocity_moments(velocity_filter_normalize(
    simulate_splicing(n_cells = 300, n_genes = 60, seed = 5), min_shared_counts = 10, n_top_genes = NULL), n_pcs = 10))))
  set.seed(42)
  a <- runif(1)
  set.seed(42)
  v <- terminal_states(v)
  expect_equal(runif(1), a)
  # velocities projected on PC2/PC3 are not drawn on PC1/PC2
  v <- velocity_embedding(v, "pca", dims = 2:3)
  expect_error(velocity_frame(v, v, "pca", dims = 1:2), "not projected on PC_1")
  vf <- velocity_frame(v, v, "pca", dims = 2:3)
  expect_equal(vf$vx, unname(v$velocity_embeddings$pca[, "PC_2"]))
})
