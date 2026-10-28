#' ggplot2 theme for embeddings
#'
#' A clean theme with no grid and light axes, close to scVelo's figures.
#' @param base_size base font size.
#' @param axes show axis titles and a light frame.
#' @return a ggplot2 theme.
#' @export
theme_velo <- function(base_size = 11, axes = TRUE) {
  th <- ggplot2::theme_void(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0, margin = ggplot2::margin(b = 6)),
      plot.subtitle = ggplot2::element_text(colour = "grey35", margin = ggplot2::margin(b = 6)),
      legend.title = ggplot2::element_text(size = base_size * 0.85, colour = "grey25"),
      legend.text = ggplot2::element_text(size = base_size * 0.8),
      legend.key.size = ggplot2::unit(0.9, "lines"),
      plot.margin = ggplot2::margin(8, 8, 8, 8),
      strip.text = ggplot2::element_text(face = "italic", size = base_size, margin = ggplot2::margin(b = 4))
    )
  if (axes) {
    th <- th + ggplot2::theme(
      axis.title = ggplot2::element_text(size = base_size * 0.8, colour = "grey40"),
      axis.title.y = ggplot2::element_text(angle = 90, size = base_size * 0.8, colour = "grey40"),
      panel.border = ggplot2::element_rect(fill = NA, colour = "grey80", linewidth = 0.4)
    )
  }
  th
}

#' Default categorical palette
#'
#' scanpy / scVelo's default 20-colour palette, extended for more groups.
#' @param n number of colours.
#' @return character vector of colours.
#' @export
velo_palette <- function(n) {
  base <- c("#1f77b4", "#ff7f0e", "#279e68", "#d62728", "#aa40fc", "#8c564b", "#e377c2",
            "#b5bd61", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
            "#c49c94", "#f7b6d2", "#dbdb8d", "#9edae5", "#ad494a", "#8c6d31")
  if (n <= length(base)) return(base[seq_len(n)])
  grDevices::hcl.colors(n, "Dark 3")
}

default_point_size <- function(n) min(3, max(0.25, sqrt(120000 / n) / 4.5))

# Data frame of cells with embedding coordinates and the colour variable.
cells_frame <- function(object, v, reduction, color_by, dims = 1:2) {
  emb <- get_embedding(object, reduction, v)[, dims, drop = FALSE]
  df <- data.frame(x = emb[, 1], y = emb[, 2], row.names = v$cells)
  if (!is.null(color_by)) {
    md <- cell_metadata(object, v)
    if (color_by %in% colnames(md)) {
      df$color <- md[[color_by]]
    } else if (color_by %in% v$genes) {
      df$color <- v$layers$Ms[color_by, ]
    } else if (is_seurat(object) && color_by %in% rownames(object)) {
      df$color <- as.numeric(SeuratObject::FetchData(object, vars = color_by, cells = v$cells)[, 1])
    } else {
      stop(sprintf("`color_by = \"%s\"` is neither a metadata column nor a gene.", color_by), call. = FALSE)
    }
    if (is.character(df$color) || is.logical(df$color)) df$color <- factor(df$color)
  }
  list(df = df, labs = colnames(emb))
}

guess_color_by <- function(object, v) {
  md <- cell_metadata(object, v)
  cand <- c("clusters", "seurat_clusters", "cell_type", "celltype", "cluster", "ident")
  hit <- cand[cand %in% colnames(md)]
  if (length(hit)) return(hit[1])
  if (is_seurat(object)) {
    md$ident <- SeuratObject::Idents(object)[v$cells]
    return(NULL)
  }
  NULL
}

cell_layer <- function(p, df, color_by, point_size, point_alpha, palette, colors_continuous,
                       legend, base_size, by = NULL) {
  if (is.null(color_by) || is.null(df$color)) {
    return(p + ggplot2::geom_point(data = df, ggplot2::aes(x = .data$x, y = .data$y),
                                   colour = "grey70", size = point_size, alpha = point_alpha, stroke = 0))
  }
  if (is.factor(df$color)) {
    lv <- levels(droplevels(df$color))
    pal <- palette %||% velo_palette(length(lv))
    if (is.null(names(pal))) pal <- stats::setNames(pal[seq_along(lv)], lv)
    p <- p + ggplot2::geom_point(data = df, ggplot2::aes(x = .data$x, y = .data$y, colour = .data$color),
                                 size = point_size, alpha = point_alpha, stroke = 0) +
      ggplot2::scale_colour_manual(values = pal, name = color_by, na.value = "grey85")
    if (identical(legend, "on data")) {
      keys <- c(list(color = df$color), lapply(stats::setNames(by, by), function(b) df[[b]]))
      cent <- stats::aggregate(df[, c("x", "y")], keys, stats::median)
      p <- p + ggplot2::geom_label(data = cent, ggplot2::aes(x = .data$x, y = .data$y, label = .data$color),
                                   fontface = "bold", size = base_size / 3.4, linewidth = 0,
                                   fill = grDevices::adjustcolor("white", 0.65),
                                   label.padding = ggplot2::unit(0.12, "lines"), colour = "grey10") +
        ggplot2::guides(colour = "none")
    } else {
      p <- p + ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1)))
    }
    if (identical(legend, "none")) p <- p + ggplot2::guides(colour = "none")
    p
  } else {
    o <- order(df$color, na.last = FALSE)
    p + ggplot2::geom_point(data = df[o, ], ggplot2::aes(x = .data$x, y = .data$y, colour = .data$color),
                            size = point_size, alpha = point_alpha, stroke = 0) +
      ggplot2::scale_colour_gradientn(colours = colors_continuous, name = color_by, na.value = "grey85") +
      ggplot2::guides(colour = ggplot2::guide_colourbar(barheight = ggplot2::unit(5, "lines"),
                                                        barwidth = ggplot2::unit(0.5, "lines")))
  }
}

continuous_colors <- function(name = "viridis") {
  switch(name,
         gnuplot = c("#000000", "#3A0083", "#7B05C2", "#B71C9A", "#E3471F", "#F9A100", "#FFFF00"),
         spectral = rev(c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")),
         grDevices::hcl.colors(9, name))
}

ensure_embedding <- function(object, reduction) {
  v <- get_velo(object)
  if (is.null(v$velocity_embeddings[[reduction]])) {
    object <- velocity_embedding(object, reduction = reduction)
    v <- get_velo(object)
  }
  v
}

#' Velocity streamlines on an embedding
#'
#' Streamlines of the projected velocity field over the cells (the classic
#' scVelo figure), drawn with ggplot2. Line width follows the local speed of
#' the field.
#'
#' @param object Seurat, SingleCellExperiment or velo object.
#' @param reduction embedding, e.g. `"umap"`.
#' @param color_by metadata column or gene used to colour cells (default:
#'   clusters if found).
#' @param density streamline density (scVelo default 2).
#' @param smooth grid smoothing (kernel width relative to grid spacing).
#' @param min_mass minimum density of cells for showing the flow.
#' @param cutoff_perc hide the flow where velocities are in the lowest percentile.
#' @param max_length maximum streamline length (axes units).
#' @param linewidth line width multiplier.
#' @param line_color streamline colour.
#' @param arrow_size arrow head size multiplier.
#' @param point_size,point_alpha cell point size and opacity.
#' @param palette named or unnamed vector of colours for groups.
#' @param continuous_palette palette name for continuous variables
#'   (`"viridis"`, `"gnuplot"`, `"spectral"`, or any [grDevices::hcl.colors()] name).
#' @param legend `"on data"`, `"right"` or `"none"`.
#' @param title plot title.
#' @param base_size base font size.
#' @param split_by cell metadata column to facet by (e.g. `"sample"`); the
#'   flow is computed separately in each panel. Defaults to the column the
#'   analysis was split by ([split_velocity()]), if any.
#' @param ncol number of facet columns.
#' @return a ggplot object.
#' @export
plot_velocity_stream <- function(object, reduction = "umap", color_by = NULL, density = 2, smooth = 0.5,
                                 min_mass = 1, cutoff_perc = 5, max_length = 4, linewidth = 1,
                                 line_color = "grey10", arrow_size = 1, point_size = NULL, point_alpha = 0.35,
                                 palette = NULL, continuous_palette = "viridis", legend = "on data",
                                 title = NULL, base_size = 11, split_by = NULL, ncol = NULL) {
  field_plot(object, "stream", reduction, color_by, split_by, ncol, point_size, point_alpha, palette,
             continuous_palette, legend, title, base_size, line_color,
             list(density = density, smooth = smooth, min_mass = min_mass, cutoff_perc = cutoff_perc,
                  max_length = max_length, linewidth = linewidth, arrow_size = arrow_size))
}

#' Velocity arrows on a grid
#'
#' @inheritParams plot_velocity_stream
#' @param arrow_length arrow length multiplier; `NULL` scales arrows so that
#'   the longest ones span about one grid cell.
#' @param arrow_color arrow colour.
#' @return a ggplot object.
#' @export
plot_velocity_grid <- function(object, reduction = "umap", color_by = NULL, density = 1, smooth = 0.5,
                               min_mass = 1, arrow_length = NULL, arrow_size = 1, arrow_color = "grey10",
                               point_size = NULL, point_alpha = 0.35, palette = NULL,
                               continuous_palette = "viridis", legend = "on data", title = NULL, base_size = 11,
                               split_by = NULL, ncol = NULL) {
  field_plot(object, "grid", reduction, color_by, split_by, ncol, point_size, point_alpha, palette,
             continuous_palette, legend, title, base_size, arrow_color,
             list(density = density, smooth = smooth, min_mass = min_mass, arrow_length = arrow_length,
                  arrow_size = arrow_size))
}

#' Velocity arrows for individual cells
#'
#' @inheritParams plot_velocity_stream
#' @param n_arrows maximum number of cells with an arrow per panel (random subset).
#' @param arrow_length arrow length multiplier; `NULL` scales the longest
#'   arrows to about 4% of the embedding width.
#' @param arrow_color arrow colour.
#' @return a ggplot object.
#' @export
plot_velocity_embedding <- function(object, reduction = "umap", color_by = NULL, n_arrows = 1500,
                                    arrow_length = NULL, arrow_size = 1, arrow_color = "grey15",
                                    point_size = NULL, point_alpha = 0.35, palette = NULL,
                                    continuous_palette = "viridis", legend = "on data", title = NULL,
                                    base_size = 11, split_by = NULL, ncol = NULL) {
  field_plot(object, "arrows", reduction, color_by, split_by, ncol, point_size, point_alpha, palette,
             continuous_palette, legend, title, base_size, arrow_color,
             list(n_arrows = n_arrows, arrow_length = arrow_length, arrow_size = arrow_size))
}

# Shared body of the three embedding plots. With `split_by` (default: the
# column the analysis was split by) the flow is computed per panel.
field_plot <- function(object, type, reduction, color_by, split_by, ncol, point_size, point_alpha, palette,
                       continuous_palette, legend, title, base_size, color, args) {
  v <- ensure_embedding(object, reduction)
  if (is.null(color_by)) color_by <- guess_color_by(object, v)
  if (is.null(split_by)) split_by <- v$params$split_by
  cf <- cells_frame(object, v, reduction, color_by)
  df <- cf$df
  vf <- velocity_frame(object, v, reduction)
  df$vx <- vf$vx
  df$vy <- vf$vy
  by <- NULL
  if (!is.null(split_by)) {
    md <- cell_metadata(object, v)
    if (!split_by %in% colnames(md)) stop("Unknown `split_by` column: ", split_by, call. = FALSE)
    df$.split <- md[[split_by]]
    df <- df[!is.na(df$.split), , drop = FALSE]
    by <- ".split"
  }
  if (is.null(point_size)) point_size <- default_point_size(nrow(df) / max(1, length(unique(df$.split))))
  p <- ggplot2::ggplot()
  p <- cell_layer(p, df, color_by, point_size, point_alpha, palette, continuous_colors(continuous_palette),
                  legend, base_size, by = by)
  p <- p + velocity_field_layers(df, type, by, args, color) +
    ggplot2::labs(x = cf$labs[1], y = cf$labs[2], title = title %||% color_by) +
    theme_velo(base_size)
  if (!is.null(by)) p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$.split), ncol = ncol)
  p
}

#' Cells on an embedding coloured by a variable
#'
#' @inheritParams plot_velocity_stream
#' @param color_by metadata column (e.g. `"latent_time"`, `"velocity_confidence"`)
#'   or gene.
#' @return a ggplot object.
#' @export
plot_cells <- function(object, color_by, reduction = "umap", point_size = NULL, point_alpha = 0.9,
                       palette = NULL, continuous_palette = "viridis", legend = "right", title = NULL,
                       base_size = 11) {
  v <- get_velo(object)
  cf <- cells_frame(object, v, reduction, color_by)
  if (is.null(point_size)) point_size <- default_point_size(nrow(cf$df))
  p <- ggplot2::ggplot()
  p <- cell_layer(p, cf$df, color_by, point_size, point_alpha, palette, continuous_colors(continuous_palette),
                  legend, base_size)
  p + ggplot2::labs(x = cf$labs[1], y = cf$labs[2], title = title %||% color_by) + theme_velo(base_size)
}

#' Phase portraits (spliced vs. unspliced)
#'
#' Unspliced against spliced moments per gene, with the steady-state ratio
#' (dashed) and, when available, the fitted trajectory of the dynamical model
#' (solid).
#'
#' @param object Seurat, SingleCellExperiment or velo object.
#' @param genes genes to show.
#' @param color_by metadata column or gene used to colour cells.
#' @param ncol number of facet columns.
#' @param point_size,point_alpha point aesthetics.
#' @param palette colours for groups.
#' @param continuous_palette palette for continuous colour variables.
#' @param base_size base font size.
#' @return a ggplot object.
#' @export
plot_phase <- function(object, genes, color_by = NULL, ncol = NULL, point_size = NULL, point_alpha = 0.6,
                       palette = NULL, continuous_palette = "viridis", base_size = 11) {
  v <- get_velo(object)
  if (is.null(color_by)) color_by <- guess_color_by(object, v)
  md <- cell_metadata(object, v)
  color_of <- function(vv) {
    if (is.null(color_by)) return(NULL)
    if (color_by %in% colnames(md)) md[vv$cells, color_by] else vv$layers$Ms[color_by, ]
  }
  split <- !is.null(v$samples)
  parts <- if (split) {
    lapply(names(v$samples), function(s) phase_frames(v$samples[[s]], genes, color_of(v$samples[[s]]), s))
  } else {
    list(phase_frames(v, genes, color_of(v), NULL))
  }
  pts <- do.call(rbind, lapply(parts, `[[`, "pts"))
  if (is.null(pts) || !nrow(pts)) stop("None of the genes are in the velocity data.", call. = FALSE)
  lines <- do.call(rbind, lapply(parts, `[[`, "lines"))
  curves <- do.call(rbind, lapply(parts, `[[`, "curves"))
  genes <- intersect(genes, unique(as.character(pts$gene)))
  col <- pts$color
  if (is.character(col) || is.logical(col)) col <- factor(col)
  if (all(is.na(col))) col <- NULL
  pts$color <- col
  lvl <- function(d) {
    if (is.null(d)) return(NULL)
    d$gene <- factor(d$gene, levels = genes)
    if (split) d$sample <- factor(d$sample, levels = names(v$samples))
    d
  }
  pts <- lvl(pts)
  lines <- lvl(lines)
  curves <- lvl(curves)
  if (is.null(point_size)) point_size <- default_point_size(nrow(pts) / length(genes)) * 0.9
  p <- ggplot2::ggplot(pts, ggplot2::aes(x = .data$s, y = .data$u))
  if (is.null(col)) {
    p <- p + ggplot2::geom_point(colour = "grey50", size = point_size, alpha = point_alpha, stroke = 0)
  } else if (is.factor(col)) {
    lv <- levels(droplevels(col))
    pal <- palette %||% velo_palette(length(lv))
    if (is.null(names(pal))) pal <- stats::setNames(pal[seq_along(lv)], lv)
    p <- p + ggplot2::geom_point(ggplot2::aes(colour = .data$color), size = point_size, alpha = point_alpha, stroke = 0) +
      ggplot2::scale_colour_manual(values = pal, name = color_by) +
      ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1)))
  } else {
    p <- p + ggplot2::geom_point(ggplot2::aes(colour = .data$color), size = point_size, alpha = point_alpha, stroke = 0) +
      ggplot2::scale_colour_gradientn(colours = continuous_colors(continuous_palette), name = color_by)
  }
  if (!is.null(lines)) {
    p <- p + ggplot2::geom_line(data = lines, colour = "#7B2CBF", linetype = "dashed", linewidth = 0.6)
  }
  if (!is.null(curves)) {
    p <- p + ggplot2::geom_path(data = curves, colour = "#7B2CBF", linewidth = 0.9)
  }
  facet <- if (split) {
    ggplot2::facet_wrap(ggplot2::vars(.data$sample, .data$gene), scales = "free", ncol = ncol %||% length(genes),
                        labeller = ggplot2::label_wrap_gen(multi_line = FALSE))
  } else {
    ggplot2::facet_wrap(~gene, scales = "free", ncol = ncol)
  }
  p + facet +
    ggplot2::labs(x = "spliced", y = "unspliced") +
    theme_velo(base_size) +
    ggplot2::theme(axis.text = ggplot2::element_text(size = base_size * 0.7, colour = "grey45"),
                   panel.spacing = ggplot2::unit(1, "lines"))
}

# Points, steady-state lines and model curves of the phase portraits of one
# velo object (one sample when the analysis is split).
phase_frames <- function(v, genes, col, sample) {
  genes <- intersect(genes, v$genes)
  if (!length(genes)) return(list())
  tag <- function(d) {
    if (!is.null(d) && !is.null(sample)) d$sample <- rep(sample, nrow(d))
    d
  }
  pts <- do.call(rbind, lapply(genes, function(g) {
    data.frame(gene = g, s = v$layers$Ms[g, ], u = v$layers$Mu[g, ], color = if (is.null(col)) NA else col)
  }))
  dyn <- identical(v$params$velocity$mode, "dynamical") && !is.null(v$var$fit_gamma)
  lines <- do.call(rbind, lapply(genes, function(g) {
    i <- match(g, v$genes)
    gamma <- if (dyn) v$var$fit_gamma[i] / v$var$fit_beta[i] else v$var$velocity_gamma[i]
    off <- if (dyn) 0 else v$var$velocity_offset[i] %||% 0
    if (is.null(gamma) || length(gamma) == 0 || is.na(gamma)) return(NULL)
    smax <- max(v$layers$Ms[g, ])
    data.frame(gene = g, s = c(0, smax), u = off + gamma * c(0, smax))
  }))
  curves <- NULL
  if (!is.null(v$var$fit_alpha)) {
    curves <- do.call(rbind, lapply(genes, function(g) {
      i <- match(g, v$genes)
      if (is.na(v$var$fit_alpha[i])) return(NULL)
      tmax <- max(v$layers$fit_t[i, ], na.rm = TRUE)
      tr <- rs_gene_trajectory(v$var$fit_alpha[i], v$var$fit_beta[i], v$var$fit_gamma[i], v$var$fit_t_[i],
                               v$var$fit_scaling[i], tmax, 300L)
      data.frame(gene = g, s = tr$s, u = tr$u)
    }))
  }
  list(pts = tag(pts), lines = tag(lines), curves = tag(curves))
}

#' Expression, velocity and phase portrait of a gene
#'
#' @inheritParams plot_velocity_stream
#' @param gene a gene name.
#' @param sample sample to show when the analysis is split by sample.
#' @return a patchwork object (or a list of ggplots if patchwork is missing).
#' @export
plot_velocity_gene <- function(object, gene, reduction = "umap", color_by = NULL, base_size = 11, sample = NULL) {
  if (is_split(object)) object <- pick_sample(object, sample)
  v <- get_velo(object)
  if (!gene %in% v$genes) stop("Gene not in velocity data: ", gene, call. = FALSE)
  emb <- get_embedding(object, reduction, v)
  df <- data.frame(x = emb[, 1], y = emb[, 2], expr = v$layers$Ms[gene, ], vel = v$layers$velocity[gene, ])
  ps <- default_point_size(nrow(df))
  lim <- stats::quantile(abs(df$vel), 0.98, na.rm = TRUE, names = FALSE)
  df$vel_c <- pmax(pmin(df$vel, lim), -lim)
  p1 <- ggplot2::ggplot(df[order(df$expr), ], ggplot2::aes(.data$x, .data$y, colour = .data$expr)) +
    ggplot2::geom_point(size = ps, stroke = 0) +
    ggplot2::scale_colour_gradientn(colours = continuous_colors("viridis"), name = NULL) +
    ggplot2::labs(title = "expression (Ms)", x = colnames(emb)[1], y = colnames(emb)[2]) + theme_velo(base_size)
  p2 <- ggplot2::ggplot(df[order(abs(df$vel_c)), ], ggplot2::aes(.data$x, .data$y, colour = .data$vel_c)) +
    ggplot2::geom_point(size = ps, stroke = 0) +
    ggplot2::scale_colour_gradient2(low = "#2166AC", mid = "grey92", high = "#B2182B", midpoint = 0, name = NULL) +
    ggplot2::labs(title = "velocity", x = colnames(emb)[1], y = colnames(emb)[2]) + theme_velo(base_size)
  p3 <- plot_phase(object, gene, color_by = color_by, base_size = base_size) +
    ggplot2::theme(legend.position = "none") + ggplot2::labs(title = "phase portrait")
  if (requireNamespace("patchwork", quietly = TRUE)) {
    patchwork::wrap_plots(p1, p2, p3, nrow = 1) + patchwork::plot_annotation(
      title = gene, theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold.italic", size = base_size * 1.3)))
  } else {
    list(p1, p2, p3)
  }
}

#' Gene expression cascade along latent time
#'
#' Heatmap of smoothed, scaled spliced moments of the given genes, with cells
#' ordered by `sort_by` and genes by the time of their peak.
#'
#' @param object Seurat, SingleCellExperiment or velo object.
#' @param genes genes (default: top 50 dynamical genes by likelihood).
#' @param sort_by cell metadata column ordering the cells.
#' @param group_by optional metadata column shown as a colour bar.
#' @param n_bins number of bins along the ordering.
#' @param palette colours for groups.
#' @param base_size base font size.
#' @param sample sample to show when the analysis is split by sample.
#' @return a ggplot object.
#' @export
plot_heatmap <- function(object, genes = NULL, sort_by = "latent_time", group_by = NULL, n_bins = 100,
                         palette = NULL, base_size = 10, sample = NULL) {
  if (is_split(object)) object <- pick_sample(object, sample)
  v <- get_velo(object)
  md <- cell_metadata(object, v)
  if (!sort_by %in% colnames(md)) stop("Unknown `sort_by` column: ", sort_by, call. = FALSE)
  if (is.null(genes)) genes <- top_likelihood_genes(object, 50)
  genes <- intersect(genes, v$genes)
  tt <- md[[sort_by]]
  ok <- is.finite(tt)
  bins <- cut(rank(tt[ok], ties.method = "first"), breaks = n_bins, labels = FALSE)
  x <- v$layers$Ms[genes, ok, drop = FALSE]
  xb <- t(apply(x, 1, function(r) tapply(r, bins, mean)))
  xb <- t(apply(xb, 1, function(r) {
    r <- stats::filter(r, rep(1 / 5, 5), sides = 2)
    r[is.na(r)] <- NA
    r
  }))
  xb <- t(apply(xb, 1, function(r) {
    rr <- range(r, na.rm = TRUE)
    if (diff(rr) > 0) (r - rr[1]) / diff(rr) else r * 0
  }))
  peak <- apply(xb, 1, function(r) which.max(r))
  ord <- genes[order(peak)]
  hm <- data.frame(gene = factor(rep(genes, ncol(xb)), levels = rev(ord)), bin = rep(seq_len(ncol(xb)), each = length(genes)),
                   value = as.vector(xb))
  p <- ggplot2::ggplot(hm, ggplot2::aes(x = .data$bin, y = .data$gene, fill = .data$value)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_gradientn(colours = continuous_colors("viridis"), na.value = "white", name = "scaled\nexpression") +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::labs(x = sort_by, y = NULL) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(face = "italic", size = base_size * 0.8))
  if (!is.null(group_by) && group_by %in% colnames(md)) {
    grp <- md[[group_by]][ok]
    dom <- tapply(as.character(grp), bins, function(z) names(which.max(table(z))))
    bar <- data.frame(bin = as.integer(names(dom)), group = factor(dom, levels = levels(factor(md[[group_by]]))))
    lv <- levels(bar$group)
    pal <- palette %||% stats::setNames(velo_palette(length(lv)), lv)
    p <- p + ggplot2::geom_tile(data = bar, ggplot2::aes(x = .data$bin, y = length(genes) + 1.2, fill = NULL, colour = NULL),
                                fill = pal[as.character(bar$group)], height = 0.9, inherit.aes = FALSE) +
      ggplot2::coord_cartesian(clip = "off")
  }
  p
}

#' Proportions of spliced, unspliced and ambiguous counts
#'
#' @param object Seurat, SingleCellExperiment or velo object.
#' @param group_by metadata column for per-group bars.
#' @param base_size base font size.
#' @return a ggplot object.
#' @export
plot_proportions <- function(object, group_by = NULL, base_size = 11) {
  v <- get_velo(object)
  if (is.null(v$obs$initial_size_spliced)) stop("Run `velocity_filter_normalize()` first.", call. = FALSE)
  md <- cell_metadata(object, v)
  if (is.null(group_by)) group_by <- guess_color_by(object, v)
  tot <- data.frame(spliced = v$obs$initial_size_spliced, unspliced = v$obs$initial_size_unspliced)
  if (!is.null(v$ambiguous_raw_total)) tot$ambiguous <- v$ambiguous_raw_total
  grp <- if (!is.null(group_by)) as.character(md[[group_by]]) else rep("all cells", nrow(tot))
  agg <- stats::aggregate(tot, list(group = grp), sum)
  agg <- rbind(data.frame(group = "all cells", t(colSums(tot))), agg)
  long <- do.call(rbind, lapply(setdiff(colnames(agg), "group"), function(k) {
    data.frame(group = agg$group, type = k, frac = agg[[k]] / rowSums(agg[, -1, drop = FALSE]))
  }))
  long$group <- factor(long$group, levels = rev(unique(agg$group)))
  long$type <- factor(long$type, levels = c("ambiguous", "unspliced", "spliced")[c("ambiguous", "unspliced", "spliced") %in% long$type])
  ggplot2::ggplot(long, ggplot2::aes(x = .data$frac, y = .data$group, fill = .data$type)) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(.data$frac > 0.06, sprintf("%.0f%%", 100 * .data$frac), "")),
                       position = ggplot2::position_stack(vjust = 0.5), size = base_size / 3.6, colour = "white") +
    ggplot2::scale_fill_manual(values = c(spliced = "#3A6EA5", unspliced = "#E07A5F", ambiguous = "#B0B7BF"),
                               breaks = c("spliced", "unspliced", "ambiguous"), name = NULL) +
    ggplot2::scale_x_continuous(labels = function(x) paste0(100 * x, "%"), expand = c(0, 0)) +
    ggplot2::labs(x = "fraction of counts", y = NULL) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "top")
}

utils::globalVariables(".data")

pick_sample <- function(object, sample) {
  v <- get_velo(object)
  if (is.null(sample)) {
    stop(sprintf("The analysis is split by '%s': choose `sample` (one of %s).", v$params$split_by,
                 paste(names(v$samples), collapse = ", ")), call. = FALSE)
  }
  get_sample_velo(object, sample)
}
