#' Read a velocyto loom file
#'
#' Reads the spliced, unspliced (and ambiguous) count layers of one or more
#' `.loom` files produced by velocyto, in memory-friendly blocks, as sparse
#' genes x cells matrices.
#'
#' @param files path(s) to `.loom` files.
#' @param sample_names optional prefixes added to cell names (one per file),
#'   useful when merging samples; also stored in `cell_attrs$sample`.
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
  parts <- lapply(files, read_one_loom, layers = layers, block_size = block_size)
  combine_parts(parts, layers, sample_names)
}

# Prefix cell names per sample, record the sample, and bind several files.
combine_parts <- function(parts, layers, sample_names) {
  if (!is.null(sample_names)) stopifnot(length(sample_names) == length(parts))
  layers <- intersect(layers, unique(unlist(lapply(parts, names))))
  for (k in seq_along(parts)) {
    res <- parts[[k]]
    if (!is.null(sample_names)) {
      new <- paste0(sample_names[k], "_", colnames(res[[layers[1]]]))
      for (l in layers) if (!is.null(res[[l]])) colnames(res[[l]]) <- new
      rownames(res$cell_attrs) <- new
      res$cell_attrs$sample <- rep(sample_names[k], length(new))
      for (r in names(res$reductions)) rownames(res$reductions[[r]]) <- new
    }
    parts[[k]] <- res
  }
  if (length(parts) == 1) return(parts[[1]])
  if (anyDuplicated(unlist(lapply(parts, function(p) colnames(p[[layers[1]]]))))) {
    stop("The same cell names occur in several files: pass `sample_names` (one per file).", call. = FALSE)
  }
  genes <- Reduce(intersect, lapply(parts, function(p) rownames(p[[layers[1]]])))
  out <- list()
  for (l in layers) {
    if (any(vapply(parts, function(p) is.null(p[[l]]), logical(1)))) next
    out[[l]] <- do.call(cbind, lapply(parts, function(p) p[[l]][genes, , drop = FALSE]))
  }
  cols <- Reduce(intersect, lapply(parts, function(p) colnames(p$cell_attrs)))
  out$cell_attrs <- do.call(rbind, lapply(parts, function(p) p$cell_attrs[, cols, drop = FALSE]))
  out$gene_attrs <- parts[[1]]$gene_attrs[match(genes, rownames(parts[[1]][[layers[1]]])), , drop = FALSE]
  red <- Reduce(intersect, lapply(parts, function(p) names(p$reductions)))
  if (length(red)) {
    out$reductions <- lapply(stats::setNames(red, red), function(r) do.call(rbind, lapply(parts, function(p) p$reductions[[r]])))
  }
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
#' `BARCODE-1`, often with a sample prefix added when merging samples
#' (`s1_BARCODE-1`). Cells are matched exactly first, then by their nucleotide
#' barcode; when the same barcode occurs in several samples, the remaining
#' parts of the names (sample labels) decide. Matches are one-to-one.
#'
#' @param loom_cells cell IDs from the loom file(s).
#' @param object_cells cell names of the Seurat / SCE object.
#' @return integer vector: for each object cell, the index of the loom cell
#'   (NA if absent or ambiguous).
#' @export
match_cells <- function(loom_cells, object_cells) {
  m <- match(object_cells, loom_cells)
  if (!anyNA(m)) return(m)
  pat <- "[ACGTN]{10,}"
  core <- function(x) {
    out <- rep(NA_character_, length(x))
    hit <- regexpr(pat, x)
    out[hit > 0] <- regmatches(x, hit)
    out
  }
  # sample labels: the name without the barcode, split into tokens; tokens
  # shared by every name (e.g. the GEM suffix "1", the trailing "x") carry no
  # information and are dropped; only the names in `need` get their tokens
  tokens <- function(x, need) {
    out <- vector("list", length(x))
    if (!length(need)) return(out)
    tk <- strsplit(sub(pat, " ", x), "[^A-Za-z0-9]+")
    flat <- unlist(tk)
    ids <- rep(seq_along(tk), lengths(tk))
    first <- unique(tk[[1]][nzchar(tk[[1]])])
    common <- first[vapply(first, function(t) length(unique(ids[flat == t])) == length(x), logical(1))]
    out[need] <- lapply(tk[need], function(t) setdiff(t[nzchar(t)], common))
    out
  }
  lc <- core(loom_cells)
  oc <- core(object_cells)
  # barcodes found once on each side (the usual case) are matched in one go
  once <- function(x) !is.na(x) & !(duplicated(x) | duplicated(x, fromLast = TRUE))
  oc_once <- once(oc)
  direct <- is.na(m) & oc_once
  m[direct] <- match(oc[direct], ifelse(once(lc), lc, NA), incomparables = NA)
  # the rest share a barcode across samples: sample labels decide
  todo <- which(is.na(m) & !is.na(oc) & oc %in% lc)
  by_bc <- split(seq_along(lc), lc)[unique(oc[todo])]
  lt <- tokens(loom_cells, unlist(by_bc))
  ot <- tokens(object_cells, todo)
  for (i in todo) {
    cand <- by_bc[[oc[i]]]
    if (length(cand) > 1 || !oc_once[i]) {
      score <- vapply(cand, function(j) length(intersect(lt[[j]], ot[[i]])), numeric(1))
      best <- cand[score == max(score)]
      cand <- if (max(score) > 0 && length(best) == 1) best else integer(0)
    }
    if (length(cand) == 1) m[i] <- cand
  }
  dup <- duplicated(m, incomparables = NA) | duplicated(m, fromLast = TRUE, incomparables = NA)
  m[dup] <- NA_integer_
  m
}

#' Attach spliced / unspliced counts to a Seurat or SingleCellExperiment object
#'
#' @param object a Seurat or SingleCellExperiment object (or NULL to build a
#'   stand-alone velo object).
#' @param loom output of [read_loom()] / [read_h5ad()], or path(s) to `.loom`
#'   or `.h5ad` files.
#'   Several samples: `read_loom(c("a.loom", "b.loom"), sample_names = c("a", "b"))`
#'   with the same prefixes used in the object (e.g. Seurat's `add.cell.ids`).
#' @param min_match minimum fraction of object cells that must be found in the
#'   loom file(s).
#' @param rename optional function mapping loom cell IDs to the object's cell
#'   names, for naming schemes that cannot be matched automatically.
#' @param split_by optional cell metadata column (e.g. `"sample"`): analyse
#'   each sample separately, see [split_velocity()].
#' @return the object with a `velo` object attached (cells not found in the
#'   loom are excluded from velocity analysis).
#' @export
add_velocity_data <- function(object = NULL, loom, min_match = 0.5, rename = NULL, split_by = NULL) {
  if (is.character(loom)) loom <- read_velocity_data(loom)
  if (!is.null(rename)) {
    new <- rename(colnames(loom$spliced))
    for (l in c("spliced", "unspliced", "ambiguous")) if (!is.null(loom[[l]])) colnames(loom[[l]]) <- new
  }
  if (is.null(object)) {
    v <- new_velo(loom$spliced, loom$unspliced, loom$ambiguous)
    v$meta <- loom$cell_attrs[v$cells, , drop = FALSE]
    v$reductions <- loom$reductions %||% list()
    if (!is.null(split_by)) v <- split_velocity(v, split_by)
    return(v)
  }
  obj_cells <- if (is_seurat(object)) colnames(object) else if (is_sce(object)) colnames(object) else stop("Unsupported object.")
  m <- match_cells(colnames(loom$spliced), obj_cells)
  frac <- mean(!is.na(m))
  if (frac < min_match) {
    stop(sprintf(paste0("Only %.1f%% of the cells were found in the loom file(s). ",
                        "Check the cell names (e.g. %s vs %s): for several samples use the same ",
                        "prefixes in read_loom(sample_names = ...), or pass `rename`."),
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
  if (!is.null(loom$cell_attrs)) {
    meta <- loom$cell_attrs[m[keep], , drop = FALSE]
    rownames(meta) <- obj_cells[keep]
    v$meta <- meta
  }
  object <- set_velo(object, v, sync_obs = FALSE)
  if (!is.null(split_by)) object <- split_velocity(object, split_by)
  object
}
