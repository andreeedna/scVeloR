# End-to-end test of the R package on the pancreas data, Seurat workflow.
# Rscript validation/r_pancreas.R
suppressPackageStartupMessages({
  library(scVeloR)
  library(Seurat)
  library(ggplot2)
})
options(scVeloR.threads = 0)
dd <- "validation/data"
ref <- "validation/out/pancreas"
out <- "validation/out/figures"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
readbin <- function(name, shape) {
  sh <- jsonlite::fromJSON(file.path(ref, "shapes.json"))[[name]]
  x <- readBin(file.path(ref, paste0(name, ".bin")), "double", n = prod(sh))
  if (length(sh) == 2) matrix(x, sh[1], sh[2], byrow = TRUE) else x
}
tm <- function(label, expr) {
  t0 <- Sys.time()
  r <- force(expr)
  cat(sprintf("%-28s %6.2f s\n", label, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  r
}

# --- a Seurat object like a user would have (Cell Ranger-style names)
loom <- tm("read_loom", read_loom(file.path(dd, "pancreas_velocyto.loom")))
cells <- paste0(readLines(file.path(dd, "pancreas_cells.txt")), "-1")
counts <- loom$spliced + loom$unspliced
colnames(counts) <- cells
seu <- CreateSeuratObject(counts, min.cells = 0, min.features = 0)
seu$clusters <- factor(readLines(file.path(dd, "pancreas_clusters.txt")),
                       levels = c("Ductal", "Ngn3 low EP", "Ngn3 high EP", "Pre-endocrine", "Beta", "Alpha", "Delta", "Epsilon"))
umap <- as.matrix(read.csv(file.path(dd, "pancreas_umap.csv"), header = FALSE))
dimnames(umap) <- list(cells, c("umap_1", "umap_2"))
seu[["umap"]] <- CreateDimReducObject(umap, key = "umap_", assay = "RNA")
pal <- c("Ductal" = "#8fbc8f", "Ngn3 low EP" = "#f4a460", "Ngn3 high EP" = "#fdbf6f", "Pre-endocrine" = "#ff7f00",
         "Beta" = "#b2df8a", "Alpha" = "#1f78b4", "Delta" = "#6a3d9a", "Epsilon" = "#cab2d6")

# --- stochastic model, step by step
seu <- tm("add_velocity_data", add_velocity_data(seu, loom))
seu <- tm("filter_normalize", velocity_filter_normalize(seu, min_shared_counts = 20, n_top_genes = 2000))
v <- get_velo(seu)
ref_genes <- readLines(file.path(ref, "genes.txt"))
cat(sprintf("genes: %d, identical to scVelo selection: %s\n", length(v$genes), identical(v$genes, ref_genes)))
ref_s <- readbin("spliced")
cat(sprintf("normalised spliced max |diff| vs scVelo: %.2e\n", max(abs(t(as.matrix(v$spliced)) - ref_s))))

# exact parity: reuse scanpy's neighbour graph
knn_ref <- readbin("knn_idx") + 1
seu_p <- tm("moments (scanpy kNN)", velocity_moments(seu, knn = knn_ref))
seu_p <- tm("velocity stochastic", velocity(seu_p, mode = "stochastic"))
seu_p <- tm("velocity_graph", velocity_graph(seu_p))
seu_p <- tm("velocity_embedding", velocity_embedding(seu_p, reduction = "umap"))
vp <- get_velo(seu_p)
cmp <- function(label, a, b) cat(sprintf("  %-26s cor=%.6f  max|diff|=%.2e\n", label, cor(as.vector(a), as.vector(b)), max(abs(a - b))))
cmp("Ms", t(vp$layers$Ms), readbin("Ms"))
cmp("velocity", t(vp$layers$velocity), readbin("velocity"))
cmp("velocity_umap", vp$velocity_embeddings$umap, readbin("velocity_umap"))
seu_p <- tm("velocity_confidence", velocity_confidence(seu_p))
cmp("velocity_confidence", seu_p$velocity_confidence, readbin("velocity_confidence"))
seu_p <- tm("terminal_states", terminal_states(seu_p))
cat(sprintf("  root_cells cor=%.4f end_points cor=%.4f\n", cor(seu_p$root_cells, readbin("root_cells")),
            cor(seu_p$end_points, readbin("end_points"))))
seu_p <- tm("velocity_pseudotime", velocity_pseudotime(seu_p))
cat(sprintf("  velocity_pseudotime cor=%.4f (spearman %.4f)\n", cor(seu_p$velocity_pseudotime, readbin("velocity_pseudotime"), use = "complete"),
            cor(seu_p$velocity_pseudotime, readbin("velocity_pseudotime"), method = "spearman", use = "complete")))

# --- natural workflow: our own exact kNN on a PCA of the spliced counts
seu <- tm("moments (own PCA+kNN)", velocity_moments(seu))
seu <- tm("velocity stochastic", velocity(seu))
seu <- tm("velocity_graph", velocity_graph(seu))
seu <- tm("velocity_embedding", velocity_embedding(seu, "umap"))
cmp("velocity vs scVelo", t(get_velo(seu)$layers$velocity), readbin("velocity"))
cmp("velocity_umap vs scVelo", get_velo(seu)$velocity_embeddings$umap, readbin("velocity_umap"))

p <- tm("plot stream", plot_velocity_stream(seu, color_by = "clusters", palette = pal, title = "pancreas - stochastic"))
ggsave(file.path(out, "stream_stochastic.png"), p, width = 7, height = 6, dpi = 150, bg = "white")
ggsave(file.path(out, "grid_stochastic.png"), plot_velocity_grid(seu, color_by = "clusters", palette = pal), width = 7, height = 6, dpi = 150, bg = "white")
ggsave(file.path(out, "arrows_stochastic.png"), plot_velocity_embedding(seu, color_by = "clusters", palette = pal), width = 7, height = 6, dpi = 150, bg = "white")
ggsave(file.path(out, "proportions.png"), plot_proportions(seu, "clusters"), width = 7, height = 4, dpi = 150, bg = "white")
ggsave(file.path(out, "phase_stochastic.png"), plot_phase(seu, c("Cpe", "Gnao1", "Ins2", "Adk"), color_by = "clusters", palette = pal, ncol = 4),
       width = 12, height = 3.4, dpi = 150, bg = "white")
seu <- velocity_confidence(seu)
ggsave(file.path(out, "confidence.png"), plot_cells(seu, "velocity_confidence", continuous_palette = "spectral"), width = 6, height = 5, dpi = 150, bg = "white")

# --- dynamical model
seu <- tm("velocity dynamical", velocity(seu, mode = "dynamical"))
seu <- tm("velocity_graph (dyn)", velocity_graph(seu))
seu <- tm("velocity_embedding (dyn)", velocity_embedding(seu, "umap"))
seu <- tm("latent_time", latent_time(seu))
lt_ref <- readbin("latent_time")
cat(sprintf("  latent_time vs scVelo: pearson %.4f spearman %.4f\n", cor(seu$latent_time, lt_ref), cor(seu$latent_time, lt_ref, method = "spearman")))
ggsave(file.path(out, "stream_dynamical.png"), plot_velocity_stream(seu, color_by = "clusters", palette = pal, title = "pancreas - dynamical"),
       width = 7, height = 6, dpi = 150, bg = "white")
ggsave(file.path(out, "latent_time.png"), plot_velocity_stream(seu, color_by = "latent_time", continuous_palette = "gnuplot", point_alpha = 0.8, legend = "right"),
       width = 7, height = 6, dpi = 150, bg = "white")
top <- top_likelihood_genes(seu, 4)
ggsave(file.path(out, "phase_dynamical.png"), plot_phase(seu, top, color_by = "clusters", palette = pal, ncol = 4), width = 12, height = 3.4, dpi = 150, bg = "white")
ggsave(file.path(out, "heatmap.png"), plot_heatmap(seu, top_likelihood_genes(seu, 40), sort_by = "latent_time", group_by = "clusters", palette = pal),
       width = 6, height = 7, dpi = 150, bg = "white")
ggsave(file.path(out, "gene_panel.png"), plot_velocity_gene(seu, top[1], color_by = "clusters"), width = 12, height = 3.8, dpi = 150, bg = "white")
print(head(rank_velocity_genes(seu, "clusters", n_genes = 3), 12))
print(get_velo(seu))
saveRDS(seu, "validation/out/pancreas_seurat.rds")
cat("done\n")
