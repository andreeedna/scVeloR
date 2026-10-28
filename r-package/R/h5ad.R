#' Read spliced / unspliced counts from AnnData `.h5ad` files
#'
#' Reads count layers written by scVelo / scanpy / anndata (any version from 0.6
#' on; sparse CSR / CSC or dense), the cell (`obs`) and gene (`var`) tables and
#' the embeddings stored in `obsm` (e.g. `X_umap`), without Python.
#'
#' @param files path(s) to `.h5ad` files.
#' @param sample_names optional prefixes added to cell names (one per file);
#'   also stored in `cell_attrs$sample`.
#' @param layers layers to read (missing ones are skipped).
#' @param obsm embeddings to read from `obsm` (`NULL`: all two-dimensional ones).
#' @return a list like [read_loom()]: `spliced`, `unspliced`, `ambiguous`
#'   (dgCMatrix, genes x cells), `cell_attrs`, `gene_attrs` and `reductions`
#'   (cells x dims matrices named without the `X_` prefix, e.g. `umap`).
#' @export
#' @examples
#' \dontrun{
#' dat <- read_h5ad("adata.h5ad")
#' v <- add_velocity_data(NULL, dat)   # stand-alone analysis with the stored UMAP
#' }
read_h5ad <- function(files, sample_names = NULL, layers = c("spliced", "unspliced", "ambiguous"), obsm = NULL) {
  parts <- lapply(files, read_one_h5ad, layers = layers, obsm = obsm)
  combine_parts(parts, layers, sample_names)
}

#' Read velocity counts from `.loom` or `.h5ad` files
#'
#' Dispatches on the file extension; several files (also mixed formats) are
#' combined as samples.
#' @param files path(s) to `.loom` / `.h5ad` files.
#' @param sample_names optional per-file prefixes for cell names.
#' @param ... passed to [read_loom()] / [read_h5ad()].
#' @return see [read_loom()].
#' @export
read_velocity_data <- function(files, sample_names = NULL, ...) {
  ext <- tolower(tools::file_ext(files))
  bad <- !ext %in% c("loom", "h5ad")
  if (any(bad)) stop("Unsupported file type: ", paste(files[bad], collapse = ", "), call. = FALSE)
  layers <- c("spliced", "unspliced", "ambiguous")
  parts <- lapply(seq_along(files), function(k) {
    if (ext[k] == "loom") read_one_loom(files[k], layers, NULL) else read_one_h5ad(files[k], layers, NULL)
  })
  combine_parts(parts, layers, sample_names)
}

read_one_h5ad <- function(file, layers, obsm) {
  if (!file.exists(file)) stop("File not found: ", file, call. = FALSE)
  on.exit(rhdf5::h5closeAll(), add = TRUE)
  ls <- rhdf5::h5ls(file)
  full <- sub("^//", "/", paste(ls$group, ls$name, sep = "/"))
  obs <- h5ad_dataframe(file, "obs", ls, full)
  var <- h5ad_dataframe(file, "var", ls, full)
  cells <- make.unique(rownames(obs))
  genes <- make.unique(rownames(var))
  rownames(obs) <- cells
  rownames(var) <- genes
  out <- list(cell_attrs = obs, gene_attrs = var)
  for (l in layers) {
    path <- paste0("/layers/", l)
    if (!path %in% full) next
    m <- h5ad_matrix(file, path, ls, full, length(cells), length(genes))
    dimnames(m) <- list(genes, cells)
    out[[l]] <- m
  }
  if (is.null(out$spliced) || is.null(out$unspliced)) {
    stop("No 'spliced' / 'unspliced' layers in ", file, call. = FALSE)
  }
  red <- ls$name[ls$group == "/obsm"]
  if (!is.null(obsm)) red <- intersect(red, obsm)
  out$reductions <- list()
  for (r in red) {
    x <- tryCatch(rhdf5::h5read(file, paste0("/obsm/", r)), error = function(e) NULL)
    if (is.null(x) || !is.matrix(x) || !is.numeric(x)) next
    if (nrow(x) != length(cells)) x <- t(x)
    if (nrow(x) != length(cells)) next
    nm <- sub("^X_", "", r)
    rownames(x) <- cells
    colnames(x) <- paste0(nm, "_", seq_len(ncol(x)))
    out$reductions[[nm]] <- x
  }
  out
}

# Sparse (CSR / CSC, new or legacy encoding) or dense matrix -> genes x cells dgCMatrix.
h5ad_matrix <- function(file, path, ls, full, n_cells, n_genes) {
  info <- ls[full == path, ]
  if (info$otype == "H5I_DATASET") {
    # dense cells x genes in HDF5 -> rhdf5 returns genes x cells
    x <- rhdf5::h5read(file, path)
    if (!identical(dim(x), c(n_genes, n_cells))) x <- t(x)
    return(as_dgc(Matrix::Matrix(x, sparse = TRUE)))
  }
  att <- rhdf5::h5readAttributes(file, path)
  fmt <- att[["encoding-type"]] %||% att[["h5sparse_format"]] %||% "csr"
  fmt <- sub("_matrix$", "", fmt)
  data <- as.numeric(rhdf5::h5read(file, paste0(path, "/data")))
  indices <- as.integer(rhdf5::h5read(file, paste0(path, "/indices")))
  indptr <- as.numeric(rhdf5::h5read(file, paste0(path, "/indptr")))
  if (fmt == "csr") {
    # CSR of cells x genes is exactly the CSC layout of genes x cells
    methods::new("dgCMatrix", i = indices, p = as.integer(indptr), x = data, Dim = c(as.integer(n_genes), as.integer(n_cells)))
  } else {
    # CSC of cells x genes: columns are genes
    j <- rep.int(seq_len(n_genes), diff(indptr))
    Matrix::sparseMatrix(i = j, j = indices + 1L, x = data, dims = c(n_genes, n_cells))
  }
}

# obs / var table: new layout (group of columns) or legacy compound dataset.
h5ad_dataframe <- function(file, key, ls, full) {
  path <- paste0("/", key)
  info <- ls[full == path, ]
  if (!nrow(info)) stop("No /", key, " in ", file, call. = FALSE)
  if (info$otype == "H5I_DATASET") {
    df <- as.data.frame(rhdf5::h5read(file, path), stringsAsFactors = FALSE)
    df[] <- lapply(df, as.vector)
    idx <- intersect(c("index", "_index", "obs_names", "var_names"), colnames(df))[1]
    rn <- as.character(df[[idx]])
    df[[idx]] <- NULL
    # legacy categoricals: integer codes with categories in uns/<col>_categories
    for (col in colnames(df)) {
      cat_path <- paste0("/uns/", col, "_categories")
      if (cat_path %in% full && is.numeric(df[[col]])) {
        lv <- as.character(rhdf5::h5read(file, cat_path))
        codes <- df[[col]]
        codes[codes < 0] <- NA
        df[[col]] <- factor(lv[codes + 1L], levels = lv)
      }
    }
    rownames(df) <- rn
    return(df)
  }
  att <- rhdf5::h5readAttributes(file, path)
  idx_name <- att[["_index"]] %||% "_index"
  rn <- as.character(rhdf5::h5read(file, paste0(path, "/", idx_name)))
  cols <- att[["column-order"]]
  if (is.null(cols)) cols <- setdiff(ls$name[ls$group == path], c(idx_name, "__categories"))
  df <- data.frame(row.names = rn)
  for (col in cols) {
    cp <- paste0(path, "/", col)
    ci <- ls[full == cp, ]
    if (!nrow(ci)) next
    val <- tryCatch({
      if (ci$otype == "H5I_GROUP") {
        # anndata >= 0.8 categorical
        codes <- as.integer(rhdf5::h5read(file, paste0(cp, "/codes")))
        lv <- as.character(rhdf5::h5read(file, paste0(cp, "/categories")))
        codes[codes < 0] <- NA
        factor(lv[codes + 1L], levels = lv)
      } else {
        x <- as.vector(rhdf5::h5read(file, cp))
        cat_path <- paste0(path, "/__categories/", col)
        if (cat_path %in% full) {
          # anndata 0.7 categorical
          lv <- as.character(rhdf5::h5read(file, cat_path))
          x[x < 0] <- NA
          factor(lv[x + 1L], levels = lv)
        } else {
          x
        }
      }
    }, error = function(e) NULL)
    if (!is.null(val) && length(val) == length(rn)) df[[col]] <- val
  }
  df
}
