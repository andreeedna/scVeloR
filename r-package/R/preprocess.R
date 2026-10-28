#' Filter genes, normalise counts and select highly variable genes
#'
#' Equivalent to `scvelo.pp.filter_and_normalize`: genes need at least
#' `min_shared_counts` counts summed over cells in which they are detected both
#' spliced and unspliced; spliced and unspliced counts are scaled to the median
#' library size (computed before gene filtering); the `n_top_genes` most
#' variable genes (Seurat flavour dispersion on spliced counts) are kept.
#'
#' @param object Seurat, SingleCellExperiment or velo object with velocity data.
#' @param min_shared_counts minimum shared counts per gene.
#' @param n_top_genes number of highly variable genes to keep (NULL keeps all).
#' @param genes optional character vector of genes to use instead of the
#'   dispersion-based selection (e.g. `Seurat::VariableFeatures(seu)`).
#' @param retain_genes genes kept regardless of the filters.
#' @return the updated object.
#' @export
velocity_filter_normalize <- function(object, min_shared_counts = 20, n_top_genes = 2000,
                                      genes = NULL, retain_genes = NULL) {
  if (is_split(object)) return(apply_split(object, velocity_filter_normalize, environment()))
  v <- get_velo(object)
  if (!is.null(v$params$normalize)) {
    message("Counts already normalised; skipping.")
    return(object)
  }
  s <- as_dgc(v$spliced)
  u <- as_dgc(v$unspliced)
  size_s <- Matrix::colSums(s)
  size_u <- Matrix::colSums(u)
  v$obs$initial_size_spliced <- size_s
  v$obs$initial_size_unspliced <- size_u

  keep <- rep(TRUE, nrow(s))
  if (!is.null(min_shared_counts) && min_shared_counts > 0) {
    both <- (s > 0) & (u > 0)
    shared <- Matrix::rowSums(s * both) + Matrix::rowSums(u * both)
    keep <- shared >= min_shared_counts
  }
  if (!is.null(retain_genes)) keep <- keep | rownames(s) %in% retain_genes
  s <- s[keep, , drop = FALSE]
  u <- u[keep, , drop = FALSE]
  a <- if (!is.null(v$ambiguous)) as_dgc(v$ambiguous)[keep, , drop = FALSE] else NULL

  scale_cols <- function(m, size) {
    f <- size / max(stats::median(size), 1)
    f[f == 0] <- 1
    m %*% Matrix::Diagonal(x = 1 / f)
  }
  s <- as_dgc(scale_cols(s, size_s))
  u <- as_dgc(scale_cols(u, size_u))
  dimnames(s) <- dimnames(u) <- list(rownames(v$spliced)[keep], v$cells)

  hv <- NULL
  if (!is.null(genes)) {
    hv <- rownames(s) %in% genes
    if (!any(hv)) stop("None of the given `genes` passed filtering.", call. = FALSE)
  } else if (!is.null(n_top_genes) && n_top_genes < nrow(s)) {
    hv <- dispersion_hvg(s, n_top_genes)
  }
  if (!is.null(hv)) {
    if (!is.null(retain_genes)) hv <- hv | rownames(s) %in% retain_genes
    s <- s[hv, , drop = FALSE]
    u <- u[hv, , drop = FALSE]
    if (!is.null(a)) a <- a[hv, , drop = FALSE]
  }
  v$spliced <- s
  v$unspliced <- u
  v$ambiguous_raw_total <- if (!is.null(v$ambiguous)) Matrix::colSums(v$ambiguous) else NULL
  v$ambiguous <- a
  v$genes <- rownames(s)
  v$var <- data.frame(row.names = v$genes)
  v$layers <- list()
  v$params$normalize <- list(min_shared_counts = min_shared_counts, n_top_genes = n_top_genes,
                             n_genes_filtered = sum(keep))
  message(sprintf("Kept %d of %d genes (%d passed the shared-count filter).",
                  length(v$genes), length(keep), sum(keep)))
  set_velo(object, v)
}

# Seurat-flavour dispersion (scvelo.pp.filter_genes_dispersion, log=TRUE, n_bins=20)
dispersion_hvg <- function(x, n_top_genes, n_bins = 20) {
  n <- ncol(x)
  mu <- Matrix::rowMeans(x)
  mu2 <- Matrix::rowMeans(x * x)
  vr <- (mu2 - mu^2) * (n / (n - 1))
  mu[mu == 0] <- 1e-12
  disp <- vr / mu
  disp[disp == 0] <- NA
  disp <- log(disp)
  lmu <- log1p(mu)
  bins <- cut(lmu, breaks = n_bins)
  bmean <- as.vector(tapply(disp, bins, mean, na.rm = TRUE))
  bsd <- as.vector(tapply(disp, bins, stats::sd, na.rm = TRUE))
  one <- is.na(bsd)
  bsd[one] <- bmean[one]
  bmean[one] <- 0
  norm <- (disp - bmean[as.integer(bins)]) / bsd[as.integer(bins)]
  norm[is.na(norm)] <- 0
  norm <- as.vector(norm)
  cut_off <- sort(norm, decreasing = TRUE)[n_top_genes]
  norm >= cut_off
}

#' Neighbour graph and first/second order moments
#'
#' Computes the kNN graph in PCA space (exact, multi-threaded) and the
#' neighbourhood means of spliced (`Ms`) and unspliced (`Mu`) counts, as
#' `scvelo.pp.moments`.
#'
#' @param object Seurat, SingleCellExperiment or velo object.
#' @param n_neighbors number of neighbours (counting the cell itself, as scanpy).
#' @param n_pcs number of principal components.
#' @param reduction reduction to build neighbours on. `"auto"` uses the host
#'   object's `pca` when present (so neighbours agree with your clustering),
#'   otherwise a PCA of log-normalised spliced counts is computed.
#' @param knn optional precomputed neighbour indices (cells x k, 1-based,
#'   excluding the cell itself).
#' @return the updated object.
#' @export
velocity_moments <- function(object, n_neighbors = 30, n_pcs = 30, reduction = "auto", knn = NULL) {
  if (is_split(object)) return(apply_split(object, velocity_moments, environment()))
  v <- get_velo(object)
  if (is.null(v$params$normalize)) {
    object <- velocity_filter_normalize(object)
    v <- get_velo(object)
  }
  n <- length(v$cells)
  if (!is.null(knn)) {
    stopifnot(nrow(knn) == n)
    v$knn <- list(idx = knn, dist = matrix(1, nrow(knn), ncol(knn)))
  } else {
    emb <- NULL
    if (identical(reduction, "auto")) {
      emb <- tryCatch(get_embedding(object, "pca", v), error = function(e) NULL)
      if (is.null(emb)) {
        v$reductions$pca <- compute_pca(v$spliced, n_pcs)
        emb <- v$reductions$pca
        v$params$moments$reduction <- "pca (computed on log spliced counts)"
      } else {
        v$params$moments$reduction <- "pca (from object)"
      }
    } else {
      emb <- get_embedding(object, reduction, v)
      v$params$moments$reduction <- reduction
    }
    emb <- emb[, seq_len(min(n_pcs, ncol(emb))), drop = FALSE]
    res <- rs_knn(as.double(t(emb)), n, ncol(emb), as.integer(n_neighbors - 1L), default_threads())
    k <- res$k
    v$knn <- list(idx = matrix(res$idx + 1L, n, k, byrow = TRUE), dist = matrix(res$dist, n, k, byrow = TRUE))
  }
  ka <- knn_args(v)
  s <- as_dgc(v$spliced)
  u <- as_dgc(v$unspliced)
  mo <- rs_moments(ka$idx, ka$dist, ka$n, ka$k, s@p, s@i, s@x, u@p, u@i, u@x, nrow(s), FALSE)
  dn <- list(v$genes, v$cells)
  v$layers$Ms <- matrix(mo$Ms, nrow(s), n, dimnames = dn)
  v$layers$Mu <- matrix(mo$Mu, nrow(s), n, dimnames = dn)
  v$params$moments$n_neighbors <- n_neighbors
  v$params$moments$n_pcs <- n_pcs
  set_velo(object, v)
}

compute_pca <- function(x, n_pcs) {
  lx <- log1p(Matrix::t(x))
  n_pcs <- min(n_pcs, min(dim(lx)) - 1)
  pc <- irlba::prcomp_irlba(lx, n = n_pcs, center = TRUE, scale. = FALSE)
  emb <- pc$x
  rownames(emb) <- colnames(x)
  colnames(emb) <- paste0("PC_", seq_len(ncol(emb)))
  emb
}
