#' Root and end points of the velocity-inferred Markov process
#'
#' Stationary states of the backward (roots) and forward (end points)
#' transition matrices, smoothed over the kNN graph (`scvelo.tl.terminal_states`).
#'
#' @param object Seurat, SingleCellExperiment or velo object.
#' @param eps eigenvalue tolerance: eigenvectors with eigenvalue >= 1 - eps are used.
#' @param n_eig number of eigenpairs computed.
#' @return the updated object with `root_cells` and `end_points` in the cell metadata.
#' @export
terminal_states <- function(object, eps = 1e-3, n_eig = 10) {
  if (is_split(object)) return(apply_split(object, terminal_states, environment()))
  v <- get_velo(object)
  if (is.null(v$graph)) {
    object <- velocity_graph(object)
    v <- get_velo(object)
  }
  ka <- knn_args(v)
  n <- length(v$cells)
  cd <- rs_connectivities(ka$idx, ka$dist, ka$n, ka$k, "distances", FALSE)
  conn <- triplets_to_sparse(cd, n)
  score <- function(backward) {
    g <- v$graph
    tr <- rs_transition_matrix(g$i, g$j, g$x, g$n, g$self_transition, 10, FALSE, FALSE, backward)
    tm <- triplets_to_sparse(tr, n)
    ev <- stationary_vectors(Matrix::t(tm), eps = eps, k = n_eig)
    s <- as.numeric(Matrix::rowSums(conn %*% ev))
    s <- pmin(pmax(s, 0), stats::quantile(s, 0.98, names = FALSE))
    scale01(s)
  }
  v$obs$root_cells <- score(TRUE)
  v$obs$end_points <- score(FALSE)
  v$params$terminal_states <- list(eps = eps, graph_mode = v$params$velocity$mode)
  set_velo(object, v)
}

# eigenvectors of A (= Tᵀ) with eigenvalue >= 1 - eps, absolute, clipped to the
# 2-98 percentile range and max-normalised (scvelo.tl.terminal_states.eigs)
stationary_vectors <- function(a, eps = 1e-3, k = 10) {
  k <- min(k, nrow(a) - 2)
  set.seed(0)
  e <- RSpectra::eigs(a, k = k, which = "LR", opts = list(maxitr = 5000, tol = 1e-10))
  ord <- order(Re(e$values), decreasing = TRUE)
  vals <- Re(e$values)[ord]
  vecs <- Re(e$vectors)[, ord, drop = FALSE]
  keep <- vals >= 1 - eps
  if (!any(keep)) keep <- seq_along(vals) == 1
  vecs <- abs(vecs[, keep, drop = FALSE])
  for (j in seq_len(ncol(vecs))) {
    q <- stats::quantile(vecs[, j], c(0.02, 0.98), names = FALSE)
    x <- vecs[, j]
    x[x < q[1]] <- 0
    x <- pmin(pmax(x, 0), q[2])
    vecs[, j] <- x / max(x)
  }
  vecs
}

scale01 <- function(x) {
  ok <- is.finite(x)
  if (!any(ok)) return(x)
  x <- x - min(x[ok])
  mx <- max(x[ok])
  if (mx != 0) x / mx else x
}

#' Velocity pseudotime
#'
#' Diffusion pseudotime on the velocity graph, measured from the root cells
#' and (inversely) from the end points (`scvelo.tl.velocity_pseudotime`).
#'
#' @param object Seurat, SingleCellExperiment or velo object.
#' @param n_dcs number of diffusion components.
#' @param root,end optional cell names / indices to use as root and end cell.
#' @return the updated object with `velocity_pseudotime` in the cell metadata.
#' @export
velocity_pseudotime <- function(object, n_dcs = 10, root = NULL, end = NULL) {
  if (is_split(object)) return(apply_split(object, velocity_pseudotime, environment()))
  v <- get_velo(object)
  if (is.null(v$obs$root_cells) || !identical(v$params$terminal_states$graph_mode, v$params$velocity$mode)) {
    object <- terminal_states(object)
    v <- get_velo(object)
  }
  n <- length(v$cells)
  g <- v$graph
  kt <- rs_pseudotime_kernel(g$i, g$j, g$x, n)
  k <- triplets_to_sparse(kt, n)
  k <- (k + Matrix::t(k)) / 2
  e <- RSpectra::eigs_sym(k, k = min(n_dcs, n - 2), which = "LM", opts = list(maxitr = 5000, tol = 1e-10))
  ord <- order(e$values, decreasing = TRUE)
  evals <- e$values[ord]
  evecs <- e$vectors[, ord, drop = FALSE]
  ka <- knn_args(v)
  conn <- triplets_to_sparse(rs_connectivities(ka$idx, ka$dist, ka$n, ka$k, "moments", FALSE), n)
  pick <- function(x, score) {
    if (!is.null(x)) return(if (is.character(x)) match(x, v$cells) else as.integer(x))
    which.max(scale01(as.numeric(conn %*% score)))
  }
  iroot <- pick(root, v$obs$root_cells)
  iend <- pick(end, v$obs$end_points)
  d_root <- rs_dpt(evals, as.double(t(evecs)), n, iroot)
  d_end <- 1 - rs_dpt(evals, as.double(t(evecs)), n, iend)
  pt <- ifelse(is.finite(d_root), d_root, 0) + ifelse(is.finite(d_end), d_end, 0)
  both <- is.finite(d_root) & is.finite(d_end)
  pt[both] <- pt[both] / 2
  pt <- scale01(pt)
  pt[!is.finite(d_root) & !is.finite(d_end)] <- NA
  v$obs$velocity_pseudotime <- pt
  v$reductions$velocity_diffmap <- evecs
  set_velo(object, v)
}

#' Gene-shared latent time
#'
#' Couples gene-wise fitted times of the dynamical model to a universal
#' latent time (`scvelo.tl.latent_time`).
#'
#' @param object Seurat, SingleCellExperiment or velo object after
#'   `velocity(mode = "dynamical")`.
#' @param min_likelihood minimum likelihood of genes used.
#' @param min_confidence local coherence threshold.
#' @param root_cells optional per-cell root scores (default: from
#'   [terminal_states()] on the dynamical velocity graph).
#' @param end_points optional per-cell end point scores to add as a prior.
#' @param t_max optional total duration to rescale latent time to.
#' @return the updated object with `latent_time` in the cell metadata.
#' @export
latent_time <- function(object, min_likelihood = 0.1, min_confidence = 0.75,
                        root_cells = NULL, end_points = NULL, t_max = NULL) {
  if (is_split(object)) return(apply_split(object, latent_time, environment()))
  v <- get_velo(object)
  if (is.null(v$layers$fit_t)) stop("Run `velocity(mode = \"dynamical\")` first.", call. = FALSE)
  if (is.null(v$graph)) {
    object <- velocity_graph(object)
    v <- get_velo(object)
  }
  if (is.null(root_cells)) {
    if (is.null(v$obs$root_cells) || !identical(v$params$terminal_states$graph_mode, v$params$velocity$mode)) {
      object <- terminal_states(object)
      v <- get_velo(object)
    }
    root_cells <- v$obs$root_cells
  }
  ft <- v$layers$fit_t
  ok <- !apply(is.na(ft), 1, any)
  if (!is.null(min_likelihood)) ok <- ok & !is.na(v$var$fit_likelihood) & v$var$fit_likelihood >= min_likelihood
  if (!any(ok)) stop("No genes pass `min_likelihood`.", call. = FALSE)
  ft <- ft[ok, , drop = FALSE]
  t_sum <- colSums(ft)
  ord <- order(t_sum)
  rv <- root_cells[ord]
  rv[is.na(rv)] <- 0
  roots <- ord[rv > 1 - 1e-3]
  if (!length(roots)) {
    warning("No root cells detected; using the cells with the smallest fitted times.", call. = FALSE)
    roots <- ord
  }
  fates <- integer(0)
  if (!is.null(end_points)) {
    ordf <- order(t_sum, decreasing = TRUE)
    ev <- end_points[ordf]
    ev[is.na(ev)] <- 0
    fates <- ordf[ev > 1 - 1e-3]
  }
  ka <- knn_args(v)
  lt <- rs_latent_time(as.double(ft), length(v$cells), nrow(ft), as.integer(roots[seq_len(min(4, length(roots)))]),
                       as.integer(fates[seq_len(min(4, length(fates)))]), ka$idx, ka$dist, ka$k, min_confidence)
  if (!is.null(t_max)) lt <- lt * t_max
  v$obs$latent_time <- lt
  set_velo(object, v)
}

#' Genes with group-specific velocities
#'
#' Ranks velocity genes per group by a Welch t-statistic of their velocities
#' (group vs. rest), after scVelo's quality filters (R², spliced/unspliced
#' Spearman correlation, likelihood for the dynamical model).
#'
#' @param object Seurat, SingleCellExperiment or velo object.
#' @param group_by cell metadata column with groups (e.g. clusters).
#' @param n_genes number of genes reported per group.
#' @param min_r2,min_corr,min_likelihood quality filters.
#' @return a data frame (group, rank, gene, score).
#' @export
rank_velocity_genes <- function(object, group_by, n_genes = 10, min_r2 = 0.1, min_corr = 0.1,
                                min_likelihood = 0.1) {
  if (is_split(object)) {
    res <- collect_split(object, rank_velocity_genes, environment())
    return(do.call(rbind, lapply(names(res), function(s) if (!is.null(res[[s]])) cbind(sample = s, res[[s]]))))
  }
  v <- get_velo(object)
  md <- cell_metadata(object, v)
  if (!group_by %in% colnames(md)) stop("Unknown `group_by` column: ", group_by, call. = FALSE)
  groups <- as.factor(md[[group_by]])
  vel <- v$layers$velocity
  keep <- !apply(is.na(vel), 1, any) & v$var$velocity_genes
  sp <- vapply(seq_len(nrow(vel)), function(i) {
    if (!keep[i]) return(NA_real_)
    suppressWarnings(stats::cor(v$layers$Ms[i, ], v$layers$Mu[i, ], method = "spearman"))
  }, numeric(1))
  keep <- keep & !is.na(sp) & sp > min_corr
  if (!is.null(v$var$velocity_r2) && v$params$velocity$mode %in% c("stochastic", "deterministic")) {
    keep <- keep & v$var$velocity_r2 > min_r2
  }
  if (!is.null(v$var$fit_likelihood) && identical(v$params$velocity$mode, "dynamical")) {
    keep <- keep & !is.na(v$var$fit_likelihood) & v$var$fit_likelihood > min_likelihood
  }
  x <- vel[keep, , drop = FALSE]
  out <- lapply(levels(groups), function(gr) {
    m <- groups == gr
    m[is.na(m)] <- FALSE
    if (sum(m) < 2 || sum(!m) < 2) return(NULL)
    a <- x[, m, drop = FALSE]
    b <- x[, !m, drop = FALSE]
    score <- (rowMeans(a) - rowMeans(b)) /
      sqrt(apply(a, 1, stats::var) / ncol(a) + apply(b, 1, stats::var) / ncol(b))
    score[!is.finite(score)] <- 0
    top <- order(score, decreasing = TRUE)[seq_len(min(n_genes, length(score)))]
    data.frame(group = gr, rank = seq_along(top), gene = rownames(x)[top], score = score[top],
               row.names = NULL, stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

#' Top genes of the dynamical model by likelihood
#' @param object Seurat, SingleCellExperiment or velo object.
#' @param n number of genes.
#' @return character vector of gene names.
#' @export
top_likelihood_genes <- function(object, n = 10) {
  if (is_split(object)) return(collect_split(object, top_likelihood_genes, environment()))
  v <- get_velo(object)
  if (is.null(v$var$fit_likelihood)) stop("Run the dynamical model first.", call. = FALSE)
  l <- v$var$fit_likelihood
  v$genes[order(l, decreasing = TRUE, na.last = NA)][seq_len(min(n, sum(!is.na(l))))]
}

#' Run the complete velocity workflow
#'
#' Convenience wrapper: attach loom counts, filter / normalise, compute
#' moments, estimate velocities, build the velocity graph, project on an
#' embedding and compute confidence, terminal states, pseudotime and (for the
#' dynamical model) latent time.
#'
#' @param object Seurat or SingleCellExperiment object (or NULL for a stand-alone
#'   analysis of the loom file).
#' @param loom output of [read_loom()] or loom path(s); can be omitted when
#'   velocity data was already attached.
#' @param mode velocity model.
#' @param reduction embedding for the projection (e.g. `"umap"`).
#' @param n_top_genes,min_shared_counts,n_neighbors,n_pcs preprocessing settings.
#' @param genes optional genes to use instead of dispersion-based selection.
#' @param split_by optional cell metadata column (e.g. `"sample"`): run the
#'   whole workflow separately for each sample (see [split_velocity()]).
#' @param ... passed to [recover_dynamics()].
#' @return the updated object.
#' @export
run_velocity <- function(object = NULL, loom = NULL, mode = c("stochastic", "dynamical", "deterministic"),
                         reduction = "umap", n_top_genes = 2000, min_shared_counts = 20,
                         n_neighbors = 30, n_pcs = 30, genes = NULL, split_by = NULL, ...) {
  mode <- match.arg(mode)
  if (!is.null(loom)) object <- add_velocity_data(object, loom)
  if (!is.null(split_by)) object <- split_velocity(object, split_by)
  object <- velocity_filter_normalize(object, min_shared_counts = min_shared_counts,
                                      n_top_genes = n_top_genes, genes = genes)
  object <- velocity_moments(object, n_neighbors = n_neighbors, n_pcs = n_pcs)
  object <- velocity(object, mode = mode, ...)
  object <- velocity_graph(object)
  v <- get_velo(object)
  has_emb <- tryCatch({ get_embedding(object, reduction, v); TRUE }, error = function(e) FALSE)
  if (!has_emb && reduction == "umap" && !is_split(object) && requireNamespace("uwot", quietly = TRUE)) {
    pcs <- v$reductions$pca %||% get_embedding(object, "pca", v)
    v$reductions$umap <- uwot::umap(pcs, n_neighbors = n_neighbors, min_dist = 0.5)
    rownames(v$reductions$umap) <- v$cells
    colnames(v$reductions$umap) <- c("UMAP_1", "UMAP_2")
    object <- set_velo(object, v)
    has_emb <- TRUE
  }
  if (has_emb) object <- velocity_embedding(object, reduction = reduction)
  object <- velocity_confidence(object)
  object <- terminal_states(object)
  object <- velocity_pseudotime(object)
  if (mode == "dynamical") object <- latent_time(object)
  object
}
