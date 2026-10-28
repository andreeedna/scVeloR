#' Add velocity streamlines or arrows to an existing plot
#'
#' Overlays the velocity field on any ggplot of an embedding made from the
#' same object — in particular Seurat's [Seurat::DimPlot()] and
#' [Seurat::FeaturePlot()], including facetted (`split.by`) plots and
#' patchworks (several features / groupings). The flow is computed separately
#' for the cells shown in each panel.
#'
#' @param plot a ggplot or patchwork object.
#' @param object the Seurat / SingleCellExperiment / velo object the plot was
#'   made from, with velocities projected on `reduction`.
#' @param reduction embedding used by the plot.
#' @param type `"stream"` (streamlines), `"grid"` (arrows on a grid) or
#'   `"arrows"` (one arrow per cell).
#' @param dims embedding dimensions used by the plot.
#' @param color colour of lines and arrows.
#' @param ... settings of the velocity layer: for streams `density`, `smooth`,
#'   `min_mass`, `cutoff_perc`, `max_length`, `linewidth`, `arrow_size`; for
#'   grids `density`, `smooth`, `min_mass`, `arrow_length`, `arrow_size`; for
#'   arrows `n_arrows`, `arrow_length`, `arrow_size`.
#' @return the plot with the velocity layer added.
#' @export
#' @examples
#' \dontrun{
#' DimPlot(seu, group.by = "celltype", split.by = "sample") |> add_velocity(seu)
#' FeaturePlot(seu, "latent_time") |> add_velocity(seu, type = "grid")
#' }
add_velocity <- function(plot, object, reduction = "umap", type = c("stream", "grid", "arrows"),
                         dims = 1:2, color = "grey10", ...) {
  type <- match.arg(type)
  v <- ensure_embedding(object, reduction)
  vf <- velocity_frame(object, v, reduction, dims)
  add_velocity_frame(plot, vf, type, color, list(...))
}

add_velocity_frame <- function(plot, vf, type, color, args) {
  if (inherits(plot, "patchwork")) {
    for (i in seq_along(plot)) {
      if (inherits(plot[[i]], "ggplot")) plot[[i]] <- add_velocity_frame(plot[[i]], vf, type, color, args)
    }
    return(plot)
  }
  if (!inherits(plot, "ggplot")) stop("`plot` must be a ggplot or patchwork object.", call. = FALSE)
  pd <- plot$data
  cells <- if (is.data.frame(pd) && !is.null(rownames(pd))) intersect(rownames(pd), rownames(vf)) else character(0)
  if (!length(cells)) cells <- rownames(vf)
  df <- vf[cells, , drop = FALSE]
  by <- facet_vars(plot)
  by <- by[by %in% colnames(pd)]
  for (b in by) df[[b]] <- pd[cells, b]
  plot + velocity_field_layers(df, type, by, args, color)
}

facet_vars <- function(plot) {
  f <- plot$facet
  if (inherits(f, "FacetWrap")) return(names(f$params$facets))
  if (inherits(f, "FacetGrid")) return(c(names(f$params$rows), names(f$params$cols)))
  character(0)
}

#' DimPlot with RNA velocity
#'
#' [Seurat::DimPlot()] with velocity streamlines (or arrows) on top. All
#' DimPlot arguments work (`group.by`, `split.by`, `cols`, `label`, `pt.size`,
#' `shuffle`, `ncol`, ...); with `split.by` the flow is computed per panel.
#' When the analysis was split by sample ([split_velocity()]) the plot is split
#' by that column unless `split.by` is given.
#'
#' @param object Seurat object with velocities projected on `reduction`.
#' @param reduction,dims embedding and dimensions.
#' @param group.by,split.by as in [Seurat::DimPlot()].
#' @param ... further [Seurat::DimPlot()] arguments.
#' @param velocity `"stream"`, `"grid"` or `"arrows"`.
#' @param velocity.args list of velocity layer settings (see [add_velocity()]).
#' @param velocity.color colour of the velocity layer.
#' @return a ggplot / patchwork object.
#' @export
#' @examples
#' \dontrun{
#' VelocityDimPlot(seu, group.by = "seurat_clusters", label = TRUE)
#' VelocityDimPlot(seu, split.by = "sample", velocity = "grid", ncol = 2)
#' }
VelocityDimPlot <- function(object, reduction = "umap", dims = c(1, 2), group.by = NULL, split.by = NULL, ...,
                            velocity = c("stream", "grid", "arrows"), velocity.args = list(),
                            velocity.color = "grey10") {
  velocity <- match.arg(velocity)
  check_seurat(object)
  v <- ensure_embedding(object, reduction)
  if (is.null(split.by)) split.by <- v$params$split_by
  if (!setequal(colnames(object), v$cells)) object <- object[, v$cells]
  p <- Seurat::DimPlot(object, reduction = reduction, dims = dims, group.by = group.by, split.by = split.by, ...)
  vf <- velocity_frame(object, v, reduction, dims)
  add_velocity_frame(p, vf, velocity, velocity.color, velocity.args)
}

#' FeaturePlot with RNA velocity
#'
#' [Seurat::FeaturePlot()] with velocity streamlines (or arrows) on top, for
#' genes or metadata such as `latent_time`, `velocity_pseudotime` or
#' `velocity_confidence`. All FeaturePlot arguments work (`split.by`, `order`,
#' `min.cutoff`, `cols`, `ncol`, ...).
#'
#' @inheritParams VelocityDimPlot
#' @param features features to plot.
#' @param ... further [Seurat::FeaturePlot()] arguments.
#' @return a ggplot / patchwork object.
#' @export
#' @examples
#' \dontrun{
#' VelocityFeaturePlot(seu, c("latent_time", "Pcsk2"), split.by = "sample")
#' }
VelocityFeaturePlot <- function(object, features, reduction = "umap", dims = c(1, 2), split.by = NULL, ...,
                                velocity = c("stream", "grid", "arrows"), velocity.args = list(),
                                velocity.color = "grey10") {
  velocity <- match.arg(velocity)
  check_seurat(object)
  v <- ensure_embedding(object, reduction)
  if (is.null(split.by)) split.by <- v$params$split_by
  if (!setequal(colnames(object), v$cells)) object <- object[, v$cells]
  p <- Seurat::FeaturePlot(object, features = features, reduction = reduction, dims = dims, split.by = split.by, ...)
  vf <- velocity_frame(object, v, reduction, dims)
  add_velocity_frame(p, vf, velocity, velocity.color, velocity.args)
}

check_seurat <- function(object) {
  if (!is_seurat(object)) stop("This function needs a Seurat object; use plot_velocity_stream() otherwise.", call. = FALSE)
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Package 'Seurat' is required.", call. = FALSE)
}
