#' Read a velocyto loom file
#'
#' Reads the spliced, unspliced (and ambiguous) count layers of one or more
#' `.loom` files produced by velocyto, in memory-friendly blocks, as sparse
#' genes x cells matrices.
#'
#' @param files path(s) to `.loom` files.
#' @param sample_names optional prefixes added to cell names (one per file),
#'   useful when merging samples.
#' @param layers layers to read.
#' @param block_size number of genes read per block.
#' @return a list with `spliced`, `unspliced`, `ambiguous` (dgCMatrix, genes x
#'   cells), `cell_attrs` and `gene_attrs` data frames.
#' @export
#' @examples
#' \dontrun{
#' loom <- read_loom("sample1.loom")
#' }
read_loom <- function(files, sample_names = NULL, layers = c("spliced", "unspliced", "ambiguous"),
                      block_size = NULL) {
  if (!is.null(sample_names)) stopifnot(length(sample_names) == length(files))
  parts <- lapply(seq_along(files), function(k) {
    res <- read_one_loom(files[k], layers, block_size)
    if (!is.null(sample_names)) {
      new <- paste0(sample_names[k], "_", colnames(res[[layers[1]]]))
      for (l in layers) if (!is.null(res[[l]])) colnames(res[[l]]) <- new
      rownames(res$cell_attrs) <- new
    }
    res
  })
  if (length(parts) == 1) return(parts[[1]])
  genes <- Reduce(intersect, lapply(parts, function(p) rownames(p[[layers[1]]])))
  out <- list()
  for (l in layers) {
    if (is.null(parts[[1]][[l]])) next
    out[[l]] <- do.call(cbind, lapply(parts, function(p) p[[l]][genes, , drop = FALSE]))
  }
  out$cell_attrs <- do.call(rbind, lapply(parts, `[[`, "cell_attrs"))
  out$gene_attrs <- parts[[1]]$gene_attrs[match(genes, rownames(parts[[1]][[layers[1]]])), , drop = FALSE]
  out
}

read_one_loom <- function(file, layers, block_size) {
  if (!file.exists(file)) stop("File not found: ", file, call. = FALSE)
  ls <- rhdf5::h5ls(file)
  paths <- paste(ls$group, ls$name, sep = "/")
  paths <- sub("^//", "/", paths)
  ra <- ls$name[ls$group == "/row_attrs"]
  ca <- ls$name[ls$group == "/col_attrs"]
  gene_key <- intersect(c("Gene", "gene_names", "var_names", "Accession"), ra)[1]
  cell_key <- intersect(c("CellID", "obs_names", "cell_names"), ca)[1]
  if (is.na(gene_key) || is.na(cell_key)) stop("Cannot find gene / cell identifiers in ", file, call. = FALSE)
  genes <- as.character(rhdf5::h5read(file, paste0("/row_attrs/", gene_key)))
  cells <- as.character(rhdf5::h5read(file, paste0("/col_attrs/", cell_key)))
  genes <- make.unique(genes)
  gene_attrs <- data.frame(row.names = genes)
  for (a in ra) {
    val <- tryCatch(rhdf5::h5read(file, paste0("/row_attrs/", a)), error = function(e) NULL)
    if (is.null(dim(val)) && length(val) == length(genes)) gene_attrs[[a]] <- as.vector(val)
  }
  cell_attrs <- data.frame(row.names = cells)
  for (a in ca) {
    val <- tryCatch(rhdf5::h5read(file, paste0("/col_attrs/", a)), error = function(e) NULL)
    if (is.null(dim(val)) && length(val) == length(cells)) cell_attrs[[a]] <- as.vector(val)
  }
  n_genes <- length(genes)
  n_cells <- length(cells)
  if (is.null(block_size)) block_size <- max(64L, as.integer(4e7 / max(1, n_cells)))
  out <- list(cell_attrs = cell_attrs, gene_attrs = gene_attrs)
  for (l in layers) {
    path <- paste0("/layers/", l)
    if (!path %in% paths) {
      if (l == "spliced" && "/matrix" %in% paths) path <- "/matrix" else next
    }
    blocks <- split(seq_len(n_genes), ceiling(seq_len(n_genes) / block_size))
    mats <- lapply(blocks, function(b) {
      # HDF5 stores genes x cells (row-major); rhdf5 returns cells x genes
      m <- rhdf5::h5read(file, path, index = list(NULL, b))
      if (!identical(dim(m), c(n_cells, length(b)))) m <- t(m)
      Matrix::t(methods::as(Matrix::Matrix(m, sparse = TRUE), "CsparseMatrix"))
    })
    m <- do.call(rbind, mats)
    m <- as_dgc(m)
    dimnames(m) <- list(genes, cells)
    out[[l]] <- m
  }
  rhdf5::h5closeAll()
  out
}

#' Match velocyto cell IDs to the cell names of an object
#'
#' velocyto names cells `sample:BARCODEx`, while Seurat / Cell Ranger use
#' `BARCODE-1` (often with a sample prefix). Cells are first matched exactly,
#' then by their nucleotide barcode, and finally by barcode plus sample label.
#'
#' @param loom_cells cell IDs from the loom file(s).
#' @param object_cells cell names of the Seurat / SCE object.
#' @return integer vector: for each object cell, the index of the loom cell
#'   (NA if absent).
#' @export
match_cells <- function(loom_cells, object_cells) {
  m <- match(object_cells, loom_cells)
  if (!anyNA(m)) return(m)
  core <- function(x) {
    bc <- regmatches(x, regexpr("[ACGTN]{8,}", x))
    out <- rep(NA_character_, length(x))
    out[grepl("[ACGTN]{8,}", x)] <- bc
    out
  }
  lc <- core(loom_cells)
  oc <- core(object_cells)
  if (!anyDuplicated(stats::na.omit(lc)) && !anyDuplicated(stats::na.omit(oc))) {
    m2 <- match(oc, lc)
    m[is.na(m)] <- m2[is.na(m)]
    return(m)
  }
  # barcode collisions across samples: also require the sample label to match
  sample_of <- function(x) sub("[:_-]?[ACGTN]{8,}.*$", "", x)
  ls <- sample_of(loom_cells)
  os <- sample_of(object_cells)
  key_l <- paste(ls, lc)
  key_o <- paste(os, oc)
  m2 <- match(key_o, key_l)
  m[is.na(m)] <- m2[is.na(m)]
  m
}

#' Attach spliced / unspliced counts to a Seurat or SingleCellExperiment object
#'
#' @param object a Seurat or SingleCellExperiment object (or NULL to build a
#'   stand-alone velo object).
#' @param loom output of [read_loom()] or path(s) to loom files.
#' @param min_match minimum fraction of object cells that must be found in the
#'   loom file(s).
#' @return the object with a `velo` object attached (cells not found in the
#'   loom are excluded from velocity analysis).
#' @export
add_velocity_data <- function(object = NULL, loom, min_match = 0.5) {
  if (is.character(loom)) loom <- read_loom(loom)
  if (is.null(object)) {
    v <- new_velo(loom$spliced, loom$unspliced, loom$ambiguous)
    v$obs <- data.frame(row.names = v$cells)
    return(v)
  }
  obj_cells <- if (is_seurat(object)) colnames(object) else if (is_sce(object)) colnames(object) else stop("Unsupported object.")
  m <- match_cells(colnames(loom$spliced), obj_cells)
  frac <- mean(!is.na(m))
  if (frac < min_match) {
    stop(sprintf(paste0("Only %.1f%% of the cells were found in the loom file(s). ",
                        "Check the cell names (e.g. %s vs %s) or set `sample_names` in read_loom()."),
                 100 * frac, obj_cells[1], colnames(loom$spliced)[1]), call. = FALSE)
  }
  if (frac < 1) message(sprintf("%d of %d cells not found in the loom file(s); they are left out of velocity analysis.",
                                sum(is.na(m)), length(m)))
  keep <- !is.na(m)
  pick <- function(x) {
    if (is.null(x)) return(NULL)
    x <- x[, m[keep], drop = FALSE]
    colnames(x) <- obj_cells[keep]
    x
  }
  v <- new_velo(pick(loom$spliced), pick(loom$unspliced), pick(loom$ambiguous))
  set_velo(object, v, sync_obs = FALSE)
}
