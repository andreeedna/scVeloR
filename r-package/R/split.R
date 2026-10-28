#' Analyse each sample separately
#'
#' Splits the velocity data by a cell metadata column (e.g. `"sample"`). Every
#' later step ([velocity_filter_normalize()], [velocity_moments()],
#' [velocity()], [velocity_graph()], [velocity_embedding()], [latent_time()],
#' ...) then runs independently on each sample — own genes, neighbours, model
#' fits and velocity graph — while the per-cell results are recombined on the
#' shared embedding, so that plots can be split by sample.
#'
#' Call it right after [add_velocity_data()] (or use its `split_by` argument),
#' before normalisation.
#'
#' @param object Seurat, SingleCellExperiment or velo object with velocity data.
#' @param split_by cell metadata column defining the samples.
#' @param min_cells samples with fewer cells are left out (with a warning).
#' @return the updated object.
#' @export
split_velocity <- function(object, split_by, min_cells = 50) {
  v <- get_velo(object)
  if (!is.null(v$params$normalize) && is.null(v$samples)) {
    stop("Split before `velocity_filter_normalize()`: re-run `add_velocity_data()` and then `split_velocity()`.", call. = FALSE)
  }
  md <- cell_metadata(object, v)
  if (!split_by %in% colnames(md)) stop("Unknown `split_by` column: ", split_by, call. = FALSE)
  grp <- as.character(md[[split_by]])
  lv <- if (is.factor(md[[split_by]])) levels(droplevels(md[[split_by]])) else sort(unique(stats::na.omit(grp)))
  samples <- list()
  for (s in lv) {
    cells <- v$cells[!is.na(grp) & grp == s]
    if (length(cells) < min_cells) {
      warning(sprintf("Sample '%s' has only %d cells and is left out.", s, length(cells)), call. = FALSE)
      next
    }
    sub <- new_velo(v$spliced[, cells, drop = FALSE], v$unspliced[, cells, drop = FALSE],
                    if (!is.null(v$ambiguous)) v$ambiguous[, cells, drop = FALSE] else NULL)
    sub$meta <- md[cells, , drop = FALSE]
    sub$params$sample <- s
    samples[[s]] <- sub
  }
  if (!length(samples)) stop("No sample has enough cells.", call. = FALSE)
  v$samples <- samples
  v$params$split_by <- split_by
  v$obs <- data.frame(row.names = v$cells)
  v$velocity_embeddings <- list()
  message(sprintf("Velocity will be analysed separately for %d samples of '%s': %s.",
                  length(samples), split_by, paste(names(samples), collapse = ", ")))
  set_velo(object, v, sync_obs = FALSE)
}

#' Velocity data of one sample
#'
#' @param object Seurat, SingleCellExperiment or velo object analysed with
#'   [split_velocity()].
#' @param sample sample name.
#' @return the sample's `velo` object.
#' @export
get_sample_velo <- function(object, sample) {
  v <- get_velo(object)
  if (is.null(v$samples)) stop("The analysis is not split by sample.", call. = FALSE)
  if (!sample %in% names(v$samples)) {
    stop(sprintf("Unknown sample '%s'. Available: %s", sample, paste(names(v$samples), collapse = ", ")), call. = FALSE)
  }
  host_sync(object, v$samples[[sample]])
}

is_split <- function(object) {
  v <- tryCatch(get_velo(object), error = function(e) NULL)
  !is.null(v) && !is.null(v$samples)
}

# Copy the host's embeddings and metadata into a sample velo object, so that
# step functions applied to it see e.g. the Seurat UMAP / PCA.
host_sync <- function(object, sub) {
  cells <- sub$cells
  if (is_seurat(object)) {
    for (r in names(object@reductions)) {
      emb <- SeuratObject::Embeddings(object, reduction = r)
      if (all(cells %in% rownames(emb))) sub$reductions[[r]] <- emb[cells, , drop = FALSE]
    }
  } else if (is_sce(object)) {
    for (r in SingleCellExperiment::reducedDimNames(object)) {
      emb <- SingleCellExperiment::reducedDim(object, r)
      rownames(emb) <- colnames(object)
      sub$reductions[[tolower(r)]] <- emb[cells, , drop = FALSE]
    }
  } else if (inherits(object, "velo")) {
    for (r in names(object$reductions)) {
      emb <- object$reductions[[r]]
      if (!is.null(rownames(emb)) && all(cells %in% rownames(emb))) sub$reductions[[r]] <- emb[cells, , drop = FALSE]
    }
  }
  md <- cell_metadata(object)
  sub$meta <- md[cells, setdiff(colnames(md), colnames(sub$obs)), drop = FALSE]
  sub
}

# Run a step function on every sample and recombine the per-cell results.
apply_split <- function(object, fun, env, ...) {
  v <- get_velo(object)
  argn <- setdiff(names(formals(fun)), c("object", "..."))
  args <- mget(argn, envir = env)
  for (s in names(v$samples)) {
    message(sprintf("-- sample %s", s))
    sub <- host_sync(object, v$samples[[s]])
    sub <- do.call(fun, c(list(sub), args, list(...)))
    v$samples[[s]] <- sub
  }
  set_velo(object, merge_split(v))
}

# Like apply_split for functions returning values (data frames, matrices, ...).
collect_split <- function(object, fun, env, ...) {
  v <- get_velo(object)
  argn <- setdiff(names(formals(fun)), c("object", "..."))
  args <- mget(argn, envir = env)
  lapply(stats::setNames(names(v$samples), names(v$samples)), function(s) {
    do.call(fun, c(list(host_sync(object, v$samples[[s]])), args, list(...)))
  })
}

merge_split <- function(v) {
  obs <- data.frame(row.names = v$cells)
  emb <- list()
  for (s in names(v$samples)) {
    sub <- v$samples[[s]]
    idx <- match(sub$cells, v$cells)
    for (col in colnames(sub$obs)) {
      if (is.null(obs[[col]])) obs[[col]] <- if (is.numeric(sub$obs[[col]])) NA_real_ else NA
      obs[[col]][idx] <- sub$obs[[col]]
    }
    for (r in names(sub$velocity_embeddings)) {
      ve <- sub$velocity_embeddings[[r]]
      if (is.null(emb[[r]])) {
        emb[[r]] <- matrix(NA_real_, length(v$cells), ncol(ve), dimnames = list(v$cells, colnames(ve)))
      }
      emb[[r]][idx, ] <- ve
    }
  }
  v$obs <- obs
  v$velocity_embeddings <- emb
  first <- v$samples[[1]]
  v$params$velocity <- first$params$velocity
  v$params$normalize <- first$params$normalize
  v
}
