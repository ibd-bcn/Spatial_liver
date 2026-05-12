# ==============================================================================
# Script: supfig6.R
# Description: Generates panels for Supplementary Figure 6 (Slingshot Trajectories & Heatmap)
# ==============================================================================

suppressPackageStartupMessages({
  library(slingshot)
  library(tradeSeq)
  library(viridis)
  library(ggplot2)
  library(dplyr)
  library(pheatmap)
  library(RColorBrewer)
})

options(bitmapType = "cairo")

# ==============================================================================
# Configuration & Paths
# Update these paths according to your local environment
# ==============================================================================
SCE_SLINGSHOT_PATH <- "/path/to/Post_neigh/sce_slingshot.rds"
SCE_TRADESEQ_PATH  <- "/path/to/Post_neigh/sce_tradeseq.rds"

OUT_DIR <- "figures/outs/"

if (!dir.exists(OUT_DIR)) {
  dir.create(OUT_DIR, recursive = TRUE)
}

# ==============================================================================
# 1. Load Data
# ==============================================================================
sce_slingshot <- readRDS(SCE_SLINGSHOT_PATH)
sce_tradeseq  <- readRDS(SCE_TRADESEQ_PATH)

# ==============================================================================
# Helper Function for Trajectory Plotting
# ==============================================================================
plot_trajectory <- function(sce_obj, traj_num, out_filename) {
  
  plot_df <- data.frame(
    UMAP_1 = reducedDim(sce_obj, 'UMAP')[, 1],
    UMAP_2 = reducedDim(sce_obj, 'UMAP')[, 2],
    pseudotime = slingPseudotime(sce_obj)[, traj_num]
  )
  
  curves <- slingCurves(sce_obj)
  curve_df <- as.data.frame(curves[[traj_num]]$s)
  colnames(curve_df) <- c("UMAP_1", "UMAP_2")
  
  p <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = pseudotime), size = 0.01, alpha = 0.8) +
    geom_path(data = curve_df, aes(x = UMAP_1, y = UMAP_2), linewidth = 1.2, color = "black") +
    scale_color_viridis_c(option = "magma", na.value = "lightgrey") +
    labs(
      title = "Slingshot Pseudotime on UMAP",
      subtitle = paste0("Trajectory ", traj_num, " with Principal Curve"),
      x = "UMAP 1",
      y = "UMAP 2",
      color = "Pseudotime"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
  
  png(filename = out_filename, width = 9.2, height = 6, units = "in", res = 800)
  print(p)
  dev.off()
}

# ==============================================================================
# Supplementary Figure 6A: Slingshot Trajectories (1-4)
# ==============================================================================
# Loop through trajectories 1 to 4 and generate plots
for (traj_idx in 1:4) {
  out_path <- file.path(OUT_DIR, paste0("supfig6A_traj", traj_idx, ".png"))
  plot_trajectory(sce_slingshot, traj_idx, out_path)
}

# ==============================================================================
# Supplementary Figure 6B: TradeSeq Expression Heatmap
# ==============================================================================
selected_clusters <- c(6, 2, 0, 1)
sel <- sce_tradeseq$leiden0.2 %in% selected_clusters

my_genes <- c(
  "SERPINA1", "XBP1", "FASN", "CEACAM1", "SERPINA3", "APOA1", 
  "COL16A1", "CXCL17", "IL7R", "NDRG1", "ICAM3", "INSR", 
  "ITGB1", "IFITM1", "VTN", "BIRC3", "TWIST1", "IL1RAP", 
  "PTEN", "IL1RN", "MET", "LINC01781", "FZD3", "CLEC14A", 
  "SMO", "EFNA1", "SRSF2", "ITGA1", "PSD3", "ZBTB16", 
  "PFN1", "NEAT1", "FN1"
)

# Extract and format count data
heatdata <- assays(sce_tradeseq)$counts[my_genes, sel]
heatclus <- factor(sce_tradeseq$leiden0.2[sel], levels = selected_clusters)

# Calculate cluster-wise average expression
heatdata_avg <- sapply(levels(heatclus), function(cl){
  rowMeans(heatdata[, heatclus == cl, drop = FALSE])
})

# Setup annotations, colors, and breaks
annotation_col <- data.frame(Cluster = factor(levels(heatclus), levels = selected_clusters))
rownames(annotation_col) <- levels(heatclus)

my_colors <- colorRampPalette(c("#053061", "white", "#67001F"))(100)
my_breaks <- seq(-2, 2, length.out = 101)

# Plot and save heatmap
png(filename = file.path(OUT_DIR, "supfig6B_heatmap.png"), width = 5, height = 8, units = "in", res = 600)
pheatmap(
  mat = log1p(heatdata_avg),       
  scale = "row",                   
  cluster_cols = FALSE,            
  cluster_rows = TRUE,             
  treeheight_row = 0,              
  color = my_colors,               
  breaks = my_breaks,              
  annotation_col = annotation_col, 
  show_rownames = TRUE,            
  main = "Average Expression by Cluster: 6 -> 2 -> 0 -> 1",
  border_color = NA
)
dev.off()