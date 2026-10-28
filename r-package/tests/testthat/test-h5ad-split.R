# Minimal AnnData (>= 0.8 layout) written with rhdf5.
write_test_h5ad <- function(path, s, u, cells, genes, sample) {
  rhdf5::h5createFile(path)
  csr <- function(grp, m) {
    # m: genes x cells dgCMatrix == CSR of cells x genes
    rhdf5::h5createGroup(path, grp)
    rhdf5::h5write(m@x, path, paste0(grp, "/data"))
    rhdf5::h5write(m@i, path, paste0(grp, "/indices"))
    rhdf5::h5write(m@p, path, paste0(grp, "/indptr"))
    fid <- rhdf5::H5Fopen(path)
    gid <- rhdf5::H5Gopen(fid, grp)
    rhdf5::h5writeAttribute("csr_matrix", gid, "encoding-type")
    rhdf5::h5writeAttribute(c(ncol(m), nrow(m)), gid, "shape")
    rhdf5::H5Gclose(gid)
    rhdf5::H5Fclose(fid)
  }
  rhdf5::h5createGroup(path, "layers")
  csr("/layers/spliced", s)
  csr("/layers/unspliced", u)
  for (key in c("obs", "var")) rhdf5::h5createGroup(path, key)
  rhdf5::h5write(cells, path, "/obs/_index")
  rhdf5::h5createGroup(path, "/obs/sample")
  lv <- sort(unique(sample))
  rhdf5::h5write(match(sample, lv) - 1L, path, "/obs/sample/codes")
  rhdf5::h5write(lv, path, "/obs/sample/categories")
  rhdf5::h5write(genes, path, "/var/_index")
  rhdf5::h5createGroup(path, "obsm")
  rhdf5::h5write(t(matrix(seq_len(2 * length(cells)), ncol = 2)), path, "/obsm/X_umap")
  fid <- rhdf5::H5Fopen(path)
  for (key in c("obs", "var")) {
    gid <- rhdf5::H5Gopen(fid, key)
    rhdf5::h5writeAttribute("_index", gid, "_index")
    if (key == "obs") rhdf5::h5writeAttribute("sample", gid, "column-order")
    rhdf5::H5Gclose(gid)
  }
  rhdf5::H5Fclose(fid)
  rhdf5::h5closeAll()
}

test_that("h5ad files are read like loom files", {
  v <- simulate_splicing(n_cells = 60, n_genes = 20, seed = 3)
  sample <- rep(c("a", "b"), 30)
  path <- tempfile(fileext = ".h5ad")
  write_test_h5ad(path, v$spliced, v$unspliced, v$cells, v$genes, sample)
  d <- read_h5ad(path)
  expect_equal(as.matrix(d$spliced), as.matrix(v$spliced), ignore_attr = TRUE)
  expect_equal(as.matrix(d$unspliced), as.matrix(v$unspliced), ignore_attr = TRUE)
  expect_equal(colnames(d$spliced), v$cells)
  expect_equal(as.character(d$cell_attrs$sample), sample)
  expect_equal(dim(d$reductions$umap), c(60, 2))
  # two files as two samples
  d2 <- read_velocity_data(c(path, path), sample_names = c("s1", "s2"))
  expect_equal(ncol(d2$spliced), 120)
  expect_equal(as.vector(table(d2$cell_attrs$sample)), c(60, 60))
  # stand-alone object keeps metadata and embeddings
  sv <- add_velocity_data(NULL, path)
  expect_equal(names(sv$reductions), "umap")
  expect_true("sample" %in% colnames(cell_metadata(sv)))
})

test_that("per-sample analysis runs every step separately", {
  v <- simulate_splicing(n_cells = 600, n_genes = 100, seed = 11)
  v$meta <- data.frame(sample = rep(c("s1", "s2"), 300), phase = rep(c("early", "late"), each = 300),
                       row.names = v$cells)
  v <- split_velocity(v, "sample")
  v <- velocity_filter_normalize(v, min_shared_counts = 10, n_top_genes = NULL)
  v <- velocity_moments(v, n_neighbors = 20, n_pcs = 15)
  v <- velocity_graph(velocity(v))
  expect_equal(names(v$samples), c("s1", "s2"))
  expect_equal(length(v$samples$s1$cells), 300)
  # each sample has its own neighbour graph over its own cells
  expect_true(all(v$samples$s1$knn$idx <= 300))
  expect_false(is.null(v$samples$s2$graph))
  v$reductions$umap <- cbind(seq_along(v$cells) %% 37, seq_along(v$cells) %% 11) + 0.0
  rownames(v$reductions$umap) <- v$cells
  v <- velocity_embedding(v, "umap")
  expect_equal(dim(v$velocity_embeddings$umap), c(600, 2))
  expect_false(anyNA(v$velocity_embeddings$umap))
  v <- velocity_confidence(v)
  expect_false(anyNA(v$obs$velocity_confidence))
  p <- plot_velocity_stream(v, color_by = "sample")
  expect_s3_class(p, "ggplot")
  expect_equal(facet_vars(p), ".split")
  expect_s3_class(plot_phase(v, v$samples$s1$genes[1:2]), "ggplot")
  expect_error(plot_heatmap(v), "choose `sample`")
  rk <- rank_velocity_genes(v, "phase", n_genes = 2)
  expect_true(all(c("sample", "group", "gene") %in% colnames(rk)))
})

test_that("add_velocity works on facetted ggplots", {
  v <- velocity_graph(velocity(velocity_moments(velocity_filter_normalize(
    simulate_splicing(n_cells = 400, n_genes = 80, seed = 5), min_shared_counts = 10, n_top_genes = NULL))))
  v$reductions$umap <- v$reductions$pca[, 1:2]
  v <- velocity_embedding(v, "umap")
  df <- data.frame(x = v$reductions$umap[, 1], y = v$reductions$umap[, 2], grp = rep(c("a", "b"), 200), row.names = v$cells)
  p <- ggplot2::ggplot(df, ggplot2::aes(x, y)) + ggplot2::geom_point() + ggplot2::facet_wrap(~grp)
  for (type in c("stream", "grid", "arrows")) {
    q <- add_velocity(p, v, type = type)
    expect_gt(length(q$layers), 1)
    expect_true("grp" %in% colnames(q$layers[[2]]$data))
    expect_silent(ggplot2::ggplot_build(q))
  }
})
