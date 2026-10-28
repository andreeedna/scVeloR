#' The velocity data object
#'
#' All scVeloR results live in a light-weight `velo` object (an S3 list). When
#' you work with Seurat or SingleCellExperiment objects you normally never touch
#' it: every function accepts the Seurat / SCE object, reads the `velo` object
#' stored inside it (`seu@misc$scVeloR` or `S4Vectors::metadata(sce)$scVeloR`),
#' updates it and writes cell-level results (latent time, confidence, ...) to
#' the object's cell metadata.
#'
#' Fields:
#' * `cells`, `genes`: identifiers (genes are the genes kept after filtering).
#' * `spliced`, `unspliced`: normalised counts, sparse genes x cells.
#' * `var`: per-gene data frame (filters, model parameters).
#' * `obs`: per-cell data frame (results).
#' * `layers`: dense genes x cells matrices (`Ms`, `Mu`, `velocity`, ...).
#' * `knn`: neighbour graph (`idx` n x k, 1-based; `dist`).
#' * `graph`: velocity graph (triplets plus self-transition probabilities).
#' * `reductions`, `velocity_embeddings`: cells x dims matrices.
#' * `params`: settings used by each step.
#'
#' @name velo
NULL

new_velo <- function(spliced, unspliced, ambiguous = NULL) {
  stopifnot(identical(dim(spliced), dim(unspliced)))
  v <- list(
    cells = colnames(spliced),
    genes = rownames(spliced),
    spliced = spliced,
    unspliced = unspliced,
    ambiguous = ambiguous,
    var = data.frame(row.names = rownames(spliced)),
    obs = data.frame(row.names = colnames(spliced)),
    layers = list(),
    knn = NULL,
    graph = NULL,
    reductions = list(),
    velocity_embeddings = list(),
    params = list()
  )
  class(v) <- "velo"
  v
}

#' @method print velo
#' @export
print.velo <- function(x, ...) {
  cat(sprintf("<velo> %d cells x %d genes\n", length(x$cells), length(x$genes)))
  if (!is.null(x$samples)) {
    cat(sprintf("  split by '%s', analysed separately:\n", x$params$split_by))
    for (s in names(x$samples)) {
      sv <- x$samples[[s]]
      cat(sprintf("    %-12s %5d cells, %4d genes, %s\n", s, length(sv$cells), length(sv$genes),
                  if (is.null(sv$params$velocity$mode)) "no velocity yet"
                  else sprintf("%s, %d velocity genes", sv$params$velocity$mode, sum(sv$var$velocity_genes, na.rm = TRUE))))
    }
  }
  if (length(x$layers)) cat("  layers:     ", paste(names(x$layers), collapse = ", "), "\n")
  if (!is.null(x$params$velocity$mode) && is.null(x$samples)) {
    cat(sprintf("  velocity:    %s model, %d velocity genes\n", x$params$velocity$mode, sum(x$var$velocity_genes, na.rm = TRUE)))
  }
  if (!is.null(x$graph)) cat(sprintf("  graph:       %d edges\n", length(x$graph$x)))
  if (length(x$reductions)) cat("  reductions: ", paste(names(x$reductions), collapse = ", "), "\n")
  if (length(x$velocity_embeddings)) cat("  embedded:   ", paste(names(x$velocity_embeddings), collapse = ", "), "\n")
  if (ncol(x$obs)) cat("  cell results:", paste(colnames(x$obs), collapse = ", "), "\n")
  invisible(x)
}

is_seurat <- function(x) inherits(x, "Seurat")
is_sce <- function(x) methods::is(x, "SingleCellExperiment")

#' Extract the velo object from Seurat / SCE / velo
#' @param object a Seurat, SingleCellExperiment or velo object
#' @return the `velo` object
#' @export
get_velo <- function(object) {
  if (inherits(object, "velo")) return(object)
  v <- if (is_seurat(object)) {
    # scVelo2: objects saved before the rename
    object@misc$scVeloR %||% object@misc$scVelo2
  } else if (is_sce(object)) {
    md <- S4Vectors::metadata(object)
    md$scVeloR %||% md$scVelo2
  } else {
    stop("`object` must be a Seurat, SingleCellExperiment or velo object.", call. = FALSE)
  }
  if (is.null(v)) {
    stop("No velocity data found. Run `add_velocity_data()` (or `run_velocity()`) first.", call. = FALSE)
  }
  v
}

#' Store a velo object inside a Seurat / SCE object
#' @param object a Seurat, SingleCellExperiment or velo object
#' @param v the velo object
#' @param sync_obs copy per-cell results to the object's cell metadata
#' @return the updated object
#' @export
set_velo <- function(object, v, sync_obs = TRUE) {
  if (inherits(object, "velo")) return(v)
  if (is_seurat(object)) {
    object@misc$scVeloR <- v
    object@misc$scVelo2 <- NULL
    if (sync_obs && ncol(v$obs)) {
      md <- object@meta.data
      idx <- match(rownames(md), v$cells)
      for (col in colnames(v$obs)) md[[col]] <- v$obs[[col]][idx]
      object@meta.data <- md
    }
  } else if (is_sce(object)) {
    md <- S4Vectors::metadata(object)
    md$scVeloR <- v
    md$scVelo2 <- NULL
    S4Vectors::metadata(object) <- md
    if (sync_obs && ncol(v$obs)) {
      cd <- SummarizedExperiment::colData(object)
      idx <- match(colnames(object), v$cells)
      for (col in colnames(v$obs)) cd[[col]] <- v$obs[[col]][idx]
      SummarizedExperiment::colData(object) <- cd
    }
  }
  object
}

# Cell metadata of the host object, aligned to the velo cells.
cell_metadata <- function(object, v = get_velo(object)) {
  md <- if (is_seurat(object)) {
    object@meta.data
  } else if (is_sce(object)) {
    as.data.frame(SummarizedExperiment::colData(object))
  } else {
    data.frame(row.names = v$cells)
  }
  md <- md[match(v$cells, rownames(md)), , drop = FALSE]
  rownames(md) <- v$cells
  # metadata read from loom / h5ad files fills in columns the host lacks
  if (!is.null(v$meta)) {
    for (col in setdiff(colnames(v$meta), colnames(md))) md[[col]] <- v$meta[[col]]
  }
  for (col in colnames(v$obs)) md[[col]] <- v$obs[[col]]
  md
}

# Embedding (cells x dims) aligned to velo cells, from the host object or velo.
get_embedding <- function(object, reduction, v = get_velo(object)) {
  emb <- NULL
  if (!is.null(v$reductions[[reduction]])) emb <- v$reductions[[reduction]]
  if (is.null(emb) && is_seurat(object) && reduction %in% names(object@reductions)) {
    emb <- SeuratObject::Embeddings(object, reduction = reduction)
  }
  if (is.null(emb) && is_sce(object)) {
    rd <- SingleCellExperiment::reducedDimNames(object)
    hit <- rd[tolower(rd) == tolower(reduction)]
    if (length(hit)) emb <- SingleCellExperiment::reducedDim(object, hit[1])
  }
  if (is.null(emb)) {
    stop(sprintf("Reduction '%s' not found. Available: %s", reduction,
                 paste(available_reductions(object, v), collapse = ", ")), call. = FALSE)
  }
  if (!is.null(rownames(emb))) {
    miss <- setdiff(v$cells, rownames(emb))
    if (length(miss)) stop(sprintf("%d cells have no coordinates in '%s'.", length(miss), reduction), call. = FALSE)
    emb <- emb[v$cells, , drop = FALSE]
  }
  emb
}

available_reductions <- function(object, v) {
  r <- names(v$reductions)
  if (is_seurat(object)) r <- c(r, names(object@reductions))
  if (is_sce(object)) r <- c(r, SingleCellExperiment::reducedDimNames(object))
  unique(r)
}

default_threads <- function() {
  n <- getOption("scVeloR.threads", 0L)
  as.integer(n)
}

# --- conversions to the Rust core -------------------------------------------

knn_args <- function(v) {
  if (is.null(v$knn)) stop("No neighbour graph. Run `velocity_moments()` first.", call. = FALSE)
  list(idx = as.integer(t(v$knn$idx) - 1L), dist = as.double(t(v$knn$dist)),
       n = nrow(v$knn$idx), k = ncol(v$knn$idx))
}

as_dgc <- function(m) {
  if (!methods::is(m, "dgCMatrix")) m <- methods::as(methods::as(m, "CsparseMatrix"), "generalMatrix")
  if (!methods::is(m, "dgCMatrix")) m <- methods::as(m, "dgCMatrix")
  m
}

triplets_to_sparse <- function(tr, n, names = NULL) {
  Matrix::sparseMatrix(i = tr$i, j = tr$j, x = tr$x, dims = c(n, n), dimnames = list(names, names))
}

vgraph_args <- function(v) {
  if (is.null(v$graph)) stop("No velocity graph. Run `velocity_graph()` first.", call. = FALSE)
  v$graph
}

`%||%` <- function(a, b) if (is.null(a)) b else a
