# ==============================================================================
# Script: supfig2.R
# Description: Generates panels for Supplementary Figure 2 (Stacked Bars, Volcano)
# ==============================================================================

# Load necessary libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
})

options(bitmapType = "cairo")

# ==============================================================================
# Configuration & Paths
# Update these paths according to your local environment
# ==============================================================================
# Input Path
SEURAT_PATH <- "/path/to/Objects/seurats_annotated.RDS"

# Output Path
OUT_DIR <- "figures/outs/"

# Ensure output directory exists
if (!dir.exists(OUT_DIR)) {
  dir.create(OUT_DIR, recursive = TRUE)
}

# ==============================================================================
# 1. Load Pre-Annotated Seurat Object & Define Palettes
# ==============================================================================
message("Loading annotated Seurat object...")
seu <- readRDS(SEURAT_PATH)
meta <- seu@meta.data

# Define Color Palettes
wide_col <- c(
  "Hepatocytes"     = "#ED9824",
  "Myeloid_cells"   = "#5387C9",
  "B_cell_lineage"  = "#772D8B",
  "T_cell_lineage"  = "#A8C686",
  "Non_parenchymal" = "#F05365"
)

refined_col <- c(
  "Gamma_delta_T_cells"       = "#304796", "NK_cells"                  = "#C2F7F5",
  "NKT_cells"                 = "#8DD1EA", "Regulatory_T_cells"        = "#3571FA",
  "Tem_Trm_cytotoxic_T_cells" = "#4ae9ff", "Effector_helper_T_cells"   = "#0091AB",
  "Naive_T_cells"             = "#BBD6DB", "Memory_B_cells"            = "#E34183",
  "Plasma_cells"              = "#F1B8EF", "Naive_B_cells"             = "#FE64F9",
  "Monocytes"                 = "#a11191", "KC1"                       = "#ff0077",
  "M2_LYVE1"                  = "#ffbf00", "KC2"                       = "#56c9f2",
  "DCs_CD1C"                  = "#02fa8d", "M1"                        = "#9c78fe",
  "Endothelial_cells_2"       = "#FFDC5F", "Fibroblasts"               = "#DB9925",
  "Endothelial_cells_1"       = "#F9F452", "Smooth_muscle_cells"       = "#FF8D08",
  "Endothelial_cells_4"       = "#CCC618", "Endothelial_cells_3"       = "#EADE8D",
  "Hepatocyte_2"              = "#c315f9", "Hepatocyte_1"              = "#E23F36",
  "Hepatocyte_6"              = "#3571FA", "Hepatocyte_3"              = "#DB9925",
  "Hepatocyte_5"              = "#20aa87", "Hepatocyte_4"              = "#CCC618",
  "Cholangiocytes"            = "#E23F36", "Hepatocytes"               = "#96307B",
  "Cycling_B_lineage_cells"   = "coral",   "Other"                     = "#393939"
)

# Set the factor order for etiology consistently
meta$etiology <- factor(meta$etiology, levels = c("HC", "HBV", "HDV RNA+", "HDV RNA-"))
seu@meta.data <- meta

# ==============================================================================
# 2. Reusable Plotting Functions
# ==============================================================================

# Helper function for Stacked Bar Plots
plot_stacked_bars <- function(df, x_var, fill_var, palette, out_file) {
  p <- ggplot(df, aes(x = !!sym(x_var), fill = !!sym(fill_var))) +
    geom_bar(position = "fill", color = "black", linewidth = 0.5) + 
    scale_fill_manual(values = palette) + 
    theme_classic(base_size = 14) +
    labs(
      title = "Distribution of Categories by Etiology",
      x = "Etiology",
      y = "Proportion",
      fill = "Category"
    ) +
    theme(
      axis.line = element_line(linewidth = 1.2, color = "black"),
      axis.ticks = element_line(linewidth = 1.2),
      axis.ticks.length = unit(0.3, "cm"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
      axis.text.y = element_text(size = 12, color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  png(filename = out_file, width = 11, height = 8, units = "in", res = 1200)
  print(p)
  dev.off()
}

# Helper function for Volcano Plots
volcano <- function(anot = "subset", ct, dif_col = "tissue", seu_obj, id1, id2) {
  
  if(ct != "all") {
    object <- subset(seu_obj, subset = !!sym(anot) == ct)
  } else {
    object <- seu_obj
  }
  
  object <- NormalizeData(object)
  object <- ScaleData(object)
  object <- SetIdent(object, value = object@meta.data[[dif_col]])
  
  deg_results <- FindMarkers(object, ident.1 = id1, ident.2 = id2)
  deg_results <- na.omit(deg_results)
  deg_results$genes <- rownames(deg_results)
  
  # Categorize Differential Expression
  deg_results$diffexpressed <- "NO"
  deg_results$diffexpressed[deg_results$avg_log2FC > log2(1.2) & deg_results$p_val < 0.05] <- "UP"
  deg_results$diffexpressed[deg_results$avg_log2FC < -log2(1.2) & deg_results$p_val < 0.05] <- "DOWN"
  deg_results$diffexpressed[deg_results$avg_log2FC > log2(1.2) & deg_results$p_val_adj < 0.05] <- "UPP"
  deg_results$diffexpressed[deg_results$avg_log2FC < -log2(1.2) & deg_results$p_val_adj < 0.05] <- "DOWNN"
  
  deg_results$delabel <- ifelse(deg_results$diffexpressed != "NO", deg_results$genes, NA)
  deg_results$p_val <- ifelse(deg_results$p_val < 1e-300, 1e-300, deg_results$p_val)
  
  diff_colors <- c("UPP" = "#803800", "DOWNN" = "#003F54", "UP" = "#B47846", "DOWN" = "steelblue")
  
  # Only label the significant adjusted p-value genes
  genes_to_label <- deg_results$delabel[deg_results$diffexpressed %in% c("UPP", "DOWNN")]
  deg_results$delabel <- ifelse(deg_results$genes %in% genes_to_label, deg_results$genes, NA)
  
  # Plot
  p <- ggplot(data = deg_results, aes(x = avg_log2FC, y = -log10(p_val), col = diffexpressed)) +
    geom_point(size = 2) +
    scale_color_manual(values = diff_colors) +
    theme_classic(base_size = 22) +
    guides(color = guide_legend(override.aes = list(shape = 1))) +
    theme(
      legend.position = "none",
      plot.title = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.line = element_line(linewidth = 2),
      axis.ticks.length = unit(0.5, "cm"),
      axis.ticks = element_line(linewidth = 2),
      text = element_text(family = "Helvetica")
    ) +
    geom_label_repel(
      aes(label = delabel),
      size = 18 / .pt,
      segment.color = "black",
      segment.size = 1,
      label.padding = unit(0.4, "lines"),
      label.size = 1,
      fontface = "bold",
      na.rm = TRUE
    )
  
  return(list(plot = p, data = deg_results))
}

# ==============================================================================
# Sup Fig 2A, 2B, 2C (Pending specific code additions)
# ==============================================================================

# ==============================================================================
# Sup Fig 2D - Stacked Bar Plots
# ==============================================================================
message("Generating Stacked Bar Plots (Sup Fig 2D)...")

# 1. Wide Categories Distribution
plot_stacked_bars(
  df = meta, 
  x_var = "etiology", 
  fill_var = "wide", 
  palette = wide_col, 
  out_file = file.path(OUT_DIR, "stackplot_supfig3.png")
)

# 2. Hepatocytes Distribution
plot_stacked_bars(
  df = meta[meta$new_anot == "Hepatocytes", ], 
  x_var = "etiology", 
  fill_var = "refined", 
  palette = refined_col, 
  out_file = file.path(OUT_DIR, "stackplot_supfig3_hepatos.png")
)

# 3. Myeloids Distribution
plot_stacked_bars(
  df = meta[meta$wide == "Myeloid_cells", ], 
  x_var = "etiology", 
  fill_var = "refined", 
  palette = refined_col, 
  out_file = file.path(OUT_DIR, "stackplot_supfig3_myeloids.png")
)

# ==============================================================================
# Sup Fig 2E - Volcano Plot
# ==============================================================================
message("Generating Volcano Plot (Sup Fig 2E)...")

volcano_res <- volcano(
  anot = "new_anot", 
  ct = "Hepatocytes", 
  dif_col = "etiology", 
  seu_obj = seu, 
  id1 = "HDV RNA+", 
  id2 = "HBV"
)

png(filename = file.path(OUT_DIR, "volcanogirl.png"), width = 11, height = 8, units = "in", res = 1200)
print(volcano_res$plot)
dev.off()

message("Supplementary Figure 2 execution complete.")