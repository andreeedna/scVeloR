# h5ad input, per-sample analysis and Seurat-style velocity plots.
# Rscript validation/r_split_seurat.R  (after validation/data/pancreas_new.h5ad exists)
suppressPackageStartupMessages({
  library(scVeloR)
  library(Seurat)
  library(ggplot2)
})
out <- "validation/out/figures_split"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
save <- function(p, name, w = 9, h = 5) ggsave(file.path(out, name), p, width = w, height = h, dpi = 120, bg = "white")
tm <- function(label, expr) {
  t0 <- Sys.time()
  r <- force(expr)
  cat(sprintf("%-32s %6.2f s\n", label, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  r
}

# --- a Seurat object built from the h5ad (Cell Ranger-style cell names)
h5 <- "validation/data/pancreas_new.h5ad"
dat <- tm("read_h5ad", read_h5ad(h5))
counts <- dat$spliced + dat$unspliced
colnames(counts) <- paste0(colnames(counts), "-1")
seu <- CreateSeuratObject(counts)
seu$clusters <- dat$cell_attrs$clusters
seu$sample <- as.character(dat$cell_attrs$sample)
umap <- dat$reductions$umap
rownames(umap) <- colnames(seu)
colnames(umap) <- c("umap_1", "umap_2")
seu[["umap"]] <- CreateDimReducObject(umap, key = "umap_", assay = "RNA")
pal <- c("Ductal" = "#8fbc8f", "Ngn3 low EP" = "#f4a460", "Ngn3 high EP" = "#fdbf6f", "Pre-endocrine" = "#ff7f00",
         "Beta" = "#b2df8a", "Alpha" = "#1f78b4", "Delta" = "#6a3d9a", "Epsilon" = "#cab2d6")
Idents(seu) <- "clusters"
seu <- NormalizeData(seu, verbose = FALSE)

# --- joint analysis straight from the h5ad path
joint <- tm("run_velocity joint (h5ad path)", run_velocity(seu, h5, mode = "stochastic"))
save(VelocityDimPlot(joint, group.by = "clusters", cols = pal, label = TRUE, alpha = 0.4) + NoLegend(),
     "joint_dimplot.png", 7, 6)
save(VelocityDimPlot(joint, group.by = "clusters", split.by = "sample", cols = pal), "joint_dimplot_split.png", 12, 5.5)
save(DimPlot(joint, group.by = "clusters", cols = pal) |> add_velocity(joint, type = "grid"), "joint_add_velocity_grid.png", 8, 6)

# --- separate analysis per sample
sep <- tm("run_velocity split dynamical", run_velocity(seu, h5, mode = "dynamical", split_by = "sample"))
print(get_velo(sep))
cat("latent_time NA:", sum(is.na(sep$latent_time)), " range:", range(sep$latent_time, na.rm = TRUE), "\n")
cat("A vs B latent time vs joint-free reference: cor by sample\n")
print(tapply(seq_len(ncol(sep)), sep$sample, function(i) cor(sep$latent_time[i], sep$velocity_pseudotime[i], method = "spearman")))
save(VelocityDimPlot(sep, group.by = "clusters", cols = pal), "split_dimplot.png", 12, 5.5)
save(VelocityFeaturePlot(sep, c("latent_time", "Pcsk2"), velocity.args = list(density = 1.5)), "split_featureplot.png", 10, 9)
save(plot_velocity_stream(sep, color_by = "clusters", palette = pal), "split_stream.png", 12, 5.5)
save(plot_phase(sep, c("Pcsk2", "Cpe", "Top2a"), color_by = "clusters", palette = pal), "split_phase.png", 10, 6)
save(plot_heatmap(sep, sample = "A", group_by = "clusters", palette = pal), "split_heatmap_A.png", 6, 7)
print(head(rank_velocity_genes(sep, "clusters", n_genes = 2), 6))
cat("done\n")
