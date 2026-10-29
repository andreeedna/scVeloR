# Velocity field layers (streamlines, grid arrows, cell arrows) as plain
# ggplot2 layers, computed separately for each panel/group so they can be
# added to facetted plots (split by sample) or to Seurat's DimPlot/FeaturePlot.

velocity_layer_defaults <- function(type) {
  switch(type,
         stream = list(density = 2, smooth = 0.5, min_mass = 1, cutoff_perc = 5, max_length = 4,
                       linewidth = 1, arrow_size = 1),
         grid = list(density = 1, smooth = 0.5, min_mass = 1, arrow_length = NULL, arrow_size = 1),
         arrows = list(n_arrows = 1500, arrow_length = NULL, arrow_size = 1),
         stop("`type` must be one of \"stream\", \"grid\", \"arrows\".", call. = FALSE))
}

# df: data frame with x, y, vx, vy (+ grouping columns `by`); returns a list of
# ggplot layers whose data carry the grouping columns.
velocity_field_layers <- function(df, type = "stream", by = NULL, args = list(), color = "grey10") {
  type <- match.arg(type, c("stream", "grid", "arrows"))
  a <- utils::modifyList(velocity_layer_defaults(type), args)
  ok <- is.finite(df$x) & is.finite(df$y) & is.finite(df$vx) & is.finite(df$vy)
  df <- df[ok, , drop = FALSE]
  groups <- if (length(by)) split(df, df[by], drop = TRUE) else list(all = df)
  tag <- function(out, g) {
    for (b in by) out[[b]] <- rep(g[[b]][1], nrow(out))
    out
  }
  res <- lapply(groups, function(g) {
    if (nrow(g) < 10) return(NULL)
    xe <- as.double(t(as.matrix(g[, c("x", "y")])))
    ve <- as.double(t(as.matrix(g[, c("vx", "vy")])))
    if (type == "stream") {
      grid <- rs_velocity_grid(xe, ve, nrow(g), 1, a$smooth, 0L, as.double(a$min_mass), FALSE, TRUE, a$cutoff_perc)
      gv <- matrix(grid$velocity, ncol = 2, byrow = TRUE)
      if (all(is.na(gv[, 1]))) return(NULL)
      st <- rs_streamlines(grid$xs, grid$ys, gv[, 1], gv[, 2], a$density, a$density, 0.1, a$max_length, "both")
      if (!length(st$line)) return(NULL)
      list(lines = tag(data.frame(line = st$line, x = st$x, y = st$y, speed = st$speed), g),
           arrows = tag(data.frame(line = st$arrow_line, x0 = st$arrow_x0, y0 = st$arrow_y0,
                                   x1 = st$arrow_x1, y1 = st$arrow_y1), g),
           lmax = max(sqrt(rowSums(gv^2)), na.rm = TRUE))
    } else if (type == "grid") {
      grid <- rs_velocity_grid(xe, ve, nrow(g), a$density, a$smooth, 0L, as.double(a$min_mass), TRUE, FALSE, 5)
      pts <- matrix(grid$points, ncol = 2, byrow = TRUE)
      vel <- matrix(grid$velocity, ncol = 2, byrow = TRUE)
      list(arrows = tag(data.frame(x = pts[, 1], y = pts[, 2], vx = vel[, 1], vy = vel[, 2]), g),
           step = mean(c(diff(grid$xs[1:2]), diff(grid$ys[1:2]))))
    } else {
      idx <- if (nrow(g) > a$n_arrows) sort(with_seed(1, sample(nrow(g), a$n_arrows))) else seq_len(nrow(g))
      list(arrows = tag(data.frame(x = g$x[idx], y = g$y[idx], vx = g$vx[idx], vy = g$vy[idx]), g),
           width = diff(range(g$x)))
    }
  })
  res <- Filter(Negate(is.null), res)
  if (!length(res)) return(list())
  head_len <- function(k) ggplot2::unit(k * a$arrow_size, "inches")
  if (type == "stream") {
    # widths are relative to the fastest flow over all panels, so panels compare
    lmax <- max(vapply(res, `[[`, numeric(1), "lmax"))
    lines <- do.call(rbind, lapply(res, `[[`, "lines"))
    arrows <- do.call(rbind, lapply(res, `[[`, "arrows"))
    # unique line ids across groups
    gid <- rep(seq_along(res), vapply(res, function(r) nrow(r$lines), integer(1)))
    lines$line <- paste(gid, lines$line)
    lines$lw <- a$linewidth * 2 * pmin(lines$speed / lmax, 1) / ggplot2::.pt * 1.5
    arrows <- arrows[(arrows$x1 - arrows$x0)^2 + (arrows$y1 - arrows$y0)^2 > 0, , drop = FALSE]
    list(
      # I(): identity widths without adding a linewidth scale to the user's plot
      ggplot2::geom_path(data = lines, ggplot2::aes(x = .data$x, y = .data$y, group = .data$line, linewidth = I(.data$lw)),
                         colour = color, lineend = "round", linejoin = "round",
                         inherit.aes = FALSE),
      ggplot2::geom_segment(data = arrows, ggplot2::aes(x = .data$x0, y = .data$y0, xend = .data$x1, yend = .data$y1),
                            colour = color, linewidth = 0.3, inherit.aes = FALSE,
                            arrow = grid::arrow(type = "closed", angle = 22, length = head_len(0.075)))
    )
  } else {
    ar <- do.call(rbind, lapply(res, `[[`, "arrows"))
    len <- a$arrow_length
    if (is.null(len)) {
      q <- stats::quantile(sqrt(ar$vx^2 + ar$vy^2), 0.95, names = FALSE, na.rm = TRUE)
      ref <- if (type == "grid") 0.95 * min(vapply(res, `[[`, numeric(1), "step"))
             else 0.04 * max(vapply(res, `[[`, numeric(1), "width"))
      len <- ref / q
    }
    ar$xend <- ar$x + ar$vx * len
    ar$yend <- ar$y + ar$vy * len
    list(ggplot2::geom_segment(data = ar, ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
                               colour = color, linewidth = if (type == "grid") 0.3 else 0.25, lineend = "round",
                               inherit.aes = FALSE,
                               arrow = grid::arrow(type = "closed", angle = if (type == "grid") 22 else 25,
                                                   length = head_len(if (type == "grid") 0.045 else 0.035))))
  }
}

# Cells with coordinates and projected velocities (x, y, vx, vy).
velocity_frame <- function(object, v, reduction, dims = 1:2) {
  emb <- get_embedding(object, reduction, v)[, dims, drop = FALSE]
  ve <- v$velocity_embeddings[[reduction]]
  if (is.null(ve)) stop(sprintf("No velocities projected on '%s'. Run `velocity_embedding(object, \"%s\")`.", reduction, reduction), call. = FALSE)
  if (is.null(colnames(emb)) || is.null(colnames(ve))) {
    # unnamed embedding (some SCE reducedDims): projection used dims 1:2
    if (ncol(ve) < max(dims)) stop("Velocity embedding has fewer dimensions than `dims`.", call. = FALSE)
    ve <- ve[, dims, drop = FALSE]
  } else {
    miss <- setdiff(colnames(emb), colnames(ve))
    if (length(miss)) {
      stop(sprintf("Velocities were not projected on %s. Run `velocity_embedding(object, \"%s\", dims = c(%s))`.",
                   paste(miss, collapse = ", "), reduction, paste(dims, collapse = ", ")), call. = FALSE)
    }
    ve <- ve[, colnames(emb), drop = FALSE]
  }
  data.frame(x = emb[, 1], y = emb[, 2], vx = ve[, 1], vy = ve[, 2], row.names = v$cells)
}
