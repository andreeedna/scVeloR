#' Estimate RNA velocity
#'
#' Gene-specific velocities from the steady-state (`"deterministic"`),
#' second-order (`"stochastic"`, default) or full splicing-kinetics model
#' (`"dynamical"`, which runs [recover_dynamics()] first if needed).
#'
#' @param object Seurat, SingleCellExperiment or velo object.
#' @param mode `"stochastic"`, `"deterministic"` or `"dynamical"`.
#' @param fit_offset,fit_offset2 fit intercepts for first / second order moments.
#' @param perc extreme quantiles used for the steady-state fit.
#' @param min_r2 minimum R² of the steady-state fit for velocity genes.
#' @param min_ratio minimum steady-state ratio for velocity genes.
#' @param min_likelihood minimum likelihood for velocity genes (dynamical).
#' @param r2_adjusted compute R² on the full-data fit rather than the quantile fit.
#' @param ... passed to [recover_dynamics()] for `mode = "dynamical"`.
#' @return the updated object; velocities are in `get_velo(x)$layers$velocity`,
#'   gene parameters in `get_velo(x)$var`.
#' @export
velocity <- function(object, mode = c("stochastic", "deterministic", "dynamical"),
                     fit_offset = FALSE, fit_offset2 = FALSE, perc = c(5, 95),
                     min_r2 = 0.01, min_ratio = 0.01, min_likelihood = 0.001,
                     r2_adjusted = FALSE, ...) {
  mode <- match.arg(mode)
  v <- get_velo(object)
  if (is.null(v$layers$Ms)) {
    object <- velocity_moments(object)
    v <- get_velo(object)
  }
  n <- length(v$cells)
  g <- length(v$genes)
  if (mode == "dynamical") {
    if (is.null(v$var$fit_alpha)) {
      object <- recover_dynamics(object, ...)
      v <- get_velo(object)
    }
    return(set_velo(object, dynamical_velocity(v, min_r2 = min_r2, min_likelihood = min_likelihood)))
  }
  mss <- mus <- numeric(0)
  if (mode == "stochastic") {
    ka <- knn_args(v)
    s <- as_dgc(v$spliced)
    u <- as_dgc(v$unspliced)
    mo <- rs_moments(ka$idx, ka$dist, ka$n, ka$k, s@p, s@i, s@x, u@p, u@i, u@x, nrow(s), TRUE)
    mss <- mo$Mss
    mus <- mo$Mus
  }
  res <- rs_velocity_steady(
    as.double(v$layers$Ms), as.double(v$layers$Mu), mss, mus, n, g,
    mode == "stochastic", fit_offset, fit_offset2,
    as.double(perc[1]), as.double(perc[length(perc)]),
    min_r2, min_ratio, r2_adjusted, integer(0), default_threads()
  )
  dn <- list(v$genes, v$cells)
  v$layers$velocity <- matrix(res$velocity, g, n, dimnames = dn)
  v$layers$variance_velocity <- if (length(res$variance_velocity)) matrix(res$variance_velocity, g, n, dimnames = dn) else NULL
  v$layers$velocity_u <- NULL
  v$var$velocity_gamma <- res$gamma
  v$var$velocity_offset <- res$offset
  v$var$velocity_r2 <- res$r2
  v$var$velocity_qreg_ratio <- res$qreg_ratio
  v$var$velocity_genes <- as.logical(res$velocity_genes)
  if (!is.na(res$low_signal_r2)) {
    warning(sprintf("Very low signal in splicing dynamics: the R² threshold was lowered to %.4f. Interpret with care.",
                    res$low_signal_r2), call. = FALSE)
  }
  if (sum(v$var$velocity_genes) < 10) {
    warning("Fewer than 10 velocity genes; consider lowering `min_r2`.", call. = FALSE)
  }
  v$params$velocity <- list(mode = mode, fit_offset = fit_offset, perc = perc)
  v$graph <- NULL
  v$velocity_embeddings <- list()
  message(sprintf("%s velocity: %d velocity genes.", mode, sum(v$var$velocity_genes)))
  set_velo(object, v)
}

#' Recover the full splicing kinetics (dynamical model)
#'
#' Fits transcription, splicing and degradation rates, switching time and
#' cell-wise latent times per gene by expectation-maximisation, following
#' `scvelo.tl.recover_dynamics`, in parallel across genes.
#'
#' @param object Seurat, SingleCellExperiment or velo object.
#' @param genes genes to fit; default: velocity genes of the steady-state model.
#' @param n_top_genes fit only the `n_top_genes` most expressed of `genes`.
#' @param max_iter EM iterations (scVelo default 10). Larger values give better
#'   fits at higher cost.
#' @param fit_scaling,fit_steady_states,fit_connected_states see scVelo.
#' @param projection use orthogonal projection for the final time assignment.
#' @param t_max total time span for aligning gene times (NULL: no alignment).
#' @return the updated object (parameters in `var$fit_*`, times in layers).
#' @export
recover_dynamics <- function(object, genes = NULL, n_top_genes = NULL, max_iter = 10,
                             fit_scaling = TRUE, fit_steady_states = TRUE,
                             fit_connected_states = TRUE, projection = TRUE, t_max = 20) {
  v <- get_velo(object)
  if (is.null(v$layers$Ms)) {
    object <- velocity_moments(object)
    v <- get_velo(object)
  }
  n <- length(v$cells)
  # steady-state velocity genes with adjusted R², as scvelo does when no
  # velocity genes were computed before
  det <- rs_velocity_steady(as.double(v$layers$Ms), as.double(v$layers$Mu), numeric(0), numeric(0),
                            n, length(v$genes), FALSE, FALSE, FALSE, 5, 95, 0.01, 0.01, TRUE,
                            integer(0), default_threads())
  v$var$fit_r2 <- det$r2
  if (is.null(genes)) {
    genes <- if (!is.null(v$var$velocity_genes)) v$genes[v$var$velocity_genes] else v$genes[det$velocity_genes == 1L]
  }
  genes <- intersect(genes, v$genes)
  if (!length(genes)) stop("No genes to fit.", call. = FALSE)
  if (!is.null(n_top_genes) && length(genes) > n_top_genes) {
    tot <- rowSums(v$layers$Ms[genes, , drop = FALSE])
    genes <- genes[order(tot, decreasing = TRUE)][seq_len(n_top_genes)]
  }
  gi <- match(genes, v$genes)
  ka <- knn_args(v)
  t0 <- Sys.time()
  res <- rs_recover_dynamics(
    as.double(v$layers$Mu[gi, , drop = FALSE]), as.double(v$layers$Ms[gi, , drop = FALSE]),
    n, length(gi), ka$idx, ka$dist, ka$k, as.integer(max_iter), fit_scaling, fit_steady_states,
    fit_connected_states, projection, if (is.null(t_max)) NA_real_ else as.double(t_max), default_threads()
  )
  message(sprintf("Recovered dynamics of %d genes in %.1f s.", length(gi),
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  nm <- c("alpha", "beta", "gamma", "t_", "scaling", "std_u", "std_s", "likelihood", "variance",
          "pval_steady", "steady_u", "steady_s", "alignment_scaling")
  for (p in nm) {
    col <- rep(NA_real_, length(v$genes))
    col[gi] <- res[[p]]
    v$var[[paste0("fit_", p)]] <- col
  }
  full <- function(x) {
    m <- matrix(NA_real_, length(v$genes), n, dimnames = list(v$genes, v$cells))
    m[gi, ] <- matrix(x, length(gi), n)
    m
  }
  v$layers$fit_t <- full(res$fit_t)
  v$layers$fit_tau <- full(res$fit_tau)
  v$layers$fit_tau_ <- full(res$fit_tau_)
  v$params$dynamics <- list(max_iter = max_iter, t_max = t_max, genes = genes)
  set_velo(object, v)
}

dynamical_velocity <- function(v, min_r2 = 0.01, min_likelihood = 0.001) {
  gi <- which(!is.na(v$var$fit_alpha))
  n <- length(v$cells)
  ka <- knn_args(v)
  vv <- v$var
  res <- rs_dynamical_velocity(
    as.double(v$layers$Mu[gi, , drop = FALSE]), as.double(v$layers$Ms[gi, , drop = FALSE]), n, length(gi),
    vv$fit_alpha[gi], vv$fit_beta[gi], vv$fit_gamma[gi], vv$fit_t_[gi], vv$fit_scaling[gi],
    vv$fit_std_u[gi], vv$fit_std_s[gi],
    as.double(v$layers$fit_tau[gi, , drop = FALSE]), as.double(v$layers$fit_tau_[gi, , drop = FALSE]),
    ka$idx, ka$dist, ka$k, TRUE
  )
  full <- function(x) {
    m <- matrix(NA_real_, length(v$genes), n, dimnames = list(v$genes, v$cells))
    m[gi, ] <- matrix(x, length(gi), n)
    m
  }
  v$layers$velocity <- full(res$velocity)
  v$layers$velocity_u <- full(res$velocity_u)
  v$layers$variance_velocity <- NULL
  vg <- !is.na(vv$fit_likelihood) & vv$fit_likelihood > min_likelihood
  if (!is.null(min_r2)) vg <- vg & !is.na(vv$fit_r2) & vv$fit_r2 > min_r2
  q <- stats::quantile(vv$fit_scaling, c(0.1, 0.9), na.rm = TRUE, names = FALSE, type = 7)
  vg <- vg & vv$fit_scaling > min(q[1], 0.03) & vv$fit_scaling < max(q[2], 3)
  vg[is.na(vg)] <- FALSE
  v$var$velocity_genes <- vg
  v$params$velocity <- list(mode = "dynamical")
  v$graph <- NULL
  v$velocity_embeddings <- list()
  message(sprintf("dynamical velocity: %d velocity genes.", sum(vg)))
  v
}

#' Velocity graph
#'
#' Cosine correlations between each cell's velocity and its displacements to
#' neighbouring cells (neighbours of neighbours by default), as
#' `scvelo.tl.velocity_graph`.
#'
#' @param object Seurat, SingleCellExperiment or velo object.
#' @param n_recurse_neighbors neighbourhood depth.
#' @param sqrt_transform variance-stabilise velocities and displacements;
#'   default: only for the stochastic model (scVelo behaviour).
#' @param max_neighbors subsample neighbourhoods larger than this.
#' @param genes genes to use (default: velocity genes).
#' @return the updated object.
#' @export
velocity_graph <- function(object, n_recurse_neighbors = 2, sqrt_transform = NULL,
                           max_neighbors = NULL, genes = NULL) {
  v <- get_velo(object)
  if (is.null(v$layers$velocity)) {
    object <- velocity(object)
    v <- get_velo(object)
  }
  if (is.null(genes)) genes <- v$genes[v$var$velocity_genes]
  gi <- match(intersect(genes, v$genes), v$genes)
  vel <- v$layers$velocity[gi, , drop = FALSE]
  ok <- !apply(is.na(vel), 1, any)
  gi <- gi[ok]
  if (is.null(sqrt_transform)) sqrt_transform <- identical(v$params$velocity$mode, "stochastic")
  ka <- knn_args(v)
  n <- length(v$cells)
  res <- rs_velocity_graph(
    as.double(v$layers$Ms[gi, , drop = FALSE]), as.double(v$layers$velocity[gi, , drop = FALSE]),
    n, length(gi), ka$idx, ka$dist, ka$k, as.integer(n_recurse_neighbors), sqrt_transform,
    as.integer(max_neighbors %||% 0L), 0L, default_threads()
  )
  v$graph <- c(res$graph, list(self_transition = res$self_transition, n = n))
  v$obs$velocity_self_transition <- res$self_transition
  v$velocity_embeddings <- list()
  v$params$graph <- list(n_recurse_neighbors = n_recurse_neighbors, sqrt_transform = sqrt_transform, n_genes = length(gi))
  set_velo(object, v)
}

#' Velocity graph as a sparse matrix
#' @param object Seurat, SingleCellExperiment or velo object.
#' @return signed dgCMatrix of cosine correlations (cells x cells).
#' @export
get_velocity_graph <- function(object) {
  v <- get_velo(object)
  g <- vgraph_args(v)
  triplets_to_sparse(g, g$n, v$cells)
}

#' Cell-to-cell transition probabilities
#'
#' @param object Seurat, SingleCellExperiment or velo object.
#' @param scale kernel width (scVelo default 10).
#' @param self_transitions allow self transitions.
#' @param use_negative_cosines include negatively correlated transitions.
#' @param backward transpose (pull back in time).
#' @return row-stochastic dgCMatrix.
#' @export
transition_matrix <- function(object, scale = 10, self_transitions = TRUE,
                              use_negative_cosines = FALSE, backward = FALSE) {
  v <- get_velo(object)
  g <- vgraph_args(v)
  tr <- rs_transition_matrix(g$i, g$j, g$x, g$n, g$self_transition, scale, self_transitions,
                             use_negative_cosines, backward)
  triplets_to_sparse(tr, g$n, v$cells)
}

#' Project velocities onto an embedding
#'
#' @param object Seurat, SingleCellExperiment or velo object.
#' @param reduction embedding name (e.g. `"umap"`, `"tsne"`, `"pca"`).
#' @param dims embedding dimensions to use.
#' @param scale kernel width of the transition matrix.
#' @param autoscale rescale arrows as scVelo does.
#' @return the updated object (`get_velo(x)$velocity_embeddings[[reduction]]`).
#' @export
velocity_embedding <- function(object, reduction = "umap", dims = 1:2, scale = 10, autoscale = TRUE) {
  v <- get_velo(object)
  if (is.null(v$graph)) {
    object <- velocity_graph(object)
    v <- get_velo(object)
  }
  emb <- get_embedding(object, reduction, v)[, dims, drop = FALSE]
  g <- v$graph
  ve <- rs_velocity_embedding(g$i, g$j, g$x, g$n, g$self_transition, scale,
                              as.double(t(emb)), ncol(emb), autoscale)
  ve <- matrix(ve, nrow(emb), ncol(emb), byrow = TRUE, dimnames = list(v$cells, colnames(emb)))
  v$velocity_embeddings[[reduction]] <- ve
  set_velo(object, v)
}

#' Velocity length and confidence
#'
#' `velocity_length`: norm of the velocity vector; `velocity_confidence`: mean
#' correlation with neighbouring velocities (local coherence);
#' `velocity_confidence_transition`: agreement between velocity and the
#' expected displacement under the transition matrix.
#' @param object Seurat, SingleCellExperiment or velo object.
#' @return the updated object with new cell metadata columns.
#' @export
velocity_confidence <- function(object) {
  v <- get_velo(object)
  if (is.null(v$graph)) {
    object <- velocity_graph(object)
    v <- get_velo(object)
  }
  gi <- which(v$var$velocity_genes & !apply(is.na(v$layers$velocity), 1, any))
  ka <- knn_args(v)
  n <- length(v$cells)
  vel <- as.double(v$layers$velocity[gi, , drop = FALSE])
  res <- rs_velocity_confidence(vel, n, length(gi), ka$idx, ka$dist, ka$k)
  v$obs$velocity_length <- res$length
  v$obs$velocity_confidence <- res$confidence
  g <- v$graph
  tr <- rs_transition_matrix(g$i, g$j, g$x, g$n, g$self_transition, 10, TRUE, FALSE, FALSE)
  v$obs$velocity_confidence_transition <- rs_confidence_transition(
    tr$i, tr$j, tr$x, n, as.double(v$layers$Ms[gi, , drop = FALSE]), vel, length(gi))
  set_velo(object, v)
}
