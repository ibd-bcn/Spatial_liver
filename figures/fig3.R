# ==============================================================================
# Script: fig3.R
# Description: Generates panels for Figure 3 (Spatial maps, Volcanos, Profiles)
# ==============================================================================

# Load necessary libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(paletteer)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(ggforce)
  library(png)
  library(grid)
})

options(bitmapType = "cairo")

# ==============================================================================
# Configuration & Paths
# Update these paths according to your local environment
# ==============================================================================
# Input Paths
SEURAT_PATH  <- "/path/to/Objects/seurats_annotated.RDS"
POLY_DIR     <- "/path/to/Polygons/"
ALIGN_DIR    <- "/path/to/Alignment/aligment/"
HDAG_DIR     <- "/path/to/Selection_regions/excels_polygons_HDAg/"

# Output Paths (Relative to the repository root)
OUT_DIR <- "figures/outs/"

# Ensure output directory exists
if (!dir.exists(OUT_DIR)) {
  dir.create(OUT_DIR, recursive = TRUE)
}

# ==============================================================================
# 1. Load Data & Define Palettes
# ==============================================================================
seu <- readRDS(SEURAT_PATH)
meta <- seu@meta.data

# Clean up Type_D NAs
meta$Type_D <- ifelse(is.na(meta$Type_D), "none", meta$Type_D)
seu@meta.data <- meta

# Palettes
wide_pal <- c(
  "Hepatocytes"     = "#ED9824",
  "Myeloid_cells"   = "#5387C9",
  "B_cell_lineage"  = "#772D8B",
  "T_cell_lineage"  = "#A8C686",
  "Non_parenchymal" = "#F05365"
)

final_pal <- c(
  "none" = "#858383",
  "+"    = "#B47846",
  "-"    = "steelblue"
)

# ==============================================================================
# 2. Reusable Functions
# ==============================================================================

# Spatial Plotting Function
plot.spatial <- function(meta, polygon_file = FALSE, x = "x", y = "y", col_fov = "fov",
                         col_cell_names = "cell_names", col_sample = "sample",
                         sample = FALSE, fov = FALSE, color = FALSE, palette = FALSE,
                         ptsize = 1, viridis = "C", scale = FALSE, per = "#000000",
                         aligment = FALSE, alpha = 1, dark = FALSE) {
  
  if (!isFALSE(fov)) meta <- meta[meta[[col_fov]] %in% fov, ]
  if (!isFALSE(sample)) meta <- meta[meta[[col_sample]] == sample, ]
  
  # Point Plot Logic
  if (isFALSE(polygon_file)) {
    if (!isFALSE(color) && isFALSE(scale)) meta[[color]] <- as.factor(meta[[color]])
    
    p <- if(!isFALSE(aligment)) aligment + coord_fixed() else ggplot(data = meta, aes(x = .data[[x]], y = .data[[y]]))
    
    if (!isFALSE(color)) {
      p <- p + geom_point(aes(color = .data[[color]]), size = ptsize, alpha = alpha)
    } else {
      p <- p + geom_point(size = ptsize, alpha = alpha)
    }
    
    # Polygon Plot Logic
  } else if (is.data.frame(polygon_file)) {
    cells <- meta[[col_cell_names]]
    poly <- polygon_file[polygon_file[[col_cell_names]] %in% cells, ]
    
    if (!isFALSE(color)) {
      poly[[color]] <- plyr::mapvalues(x = poly[[col_cell_names]], from = meta[[col_cell_names]], to = meta[[color]], warn_missing = FALSE)
      poly[[color]] <- as.factor(poly[[color]])
    }
    
    p <- if(!isFALSE(aligment)) aligment + coord_fixed() else ggplot(poly, aes(x = .data[[x]], y = .data[[y]]))
    
    if (!isFALSE(color)) {
      p <- p + geom_polygon(aes(group = .data[[col_cell_names]], fill = .data[[color]]), color = per, linewidth = ptsize, alpha = alpha)
    } else {
      p <- p + geom_polygon(aes(group = .data[[col_cell_names]]), color = per, linewidth = ptsize, alpha = alpha)
    }
  }
  
  # Theming
  if (!isFALSE(dark)) p <- p + ggdark::dark_theme_gray(base_family = "Fira Sans Condensed Light", base_size = 20)
  p <- p + theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    axis.title = element_blank(), axis.text = element_blank(),
    axis.ticks = element_blank(), panel.background = element_blank()
  ) + labs(x = "x", y = "y")
  
  if (!isFALSE(palette)) p <- p + scale_fill_manual(values = palette)
  if (!isFALSE(scale)) p <- p + scale_fill_viridis_c(option = viridis)
  
  return(p)
}

# Image Alignment Function
align_image <- function(patient, fov) {
  dig_alignment <- c("Slide_1" = "S1", "Slide_2" = "S2", "Slide_3" = "S3")
  img_path <- file.path(ALIGN_DIR, paste0(dig_alignment[patient], "_F", fov, ".png"))
  
  img <- readPNG(img_path, native = FALSE)
  img_dims <- dim(img)[1:2] 
  
  if (patient == "Slide_3") {
    x_scale <- 0.149; y_scale <- 0.13; x_shift <- 2; y_shift <- 43
  } else if (patient == "Slide_2") {
    x_scale <- 0.102; y_scale <- 0.146; x_shift <- 98; y_shift <- 0
  } else if (patient == "Slide_1") {
    x_scale <- 0.15; y_scale <- 0.13; x_shift <- 8; y_shift <- 50
  }
  
  xmin_adj <- (-x_shift) / x_scale
  xmax_adj <- (img_dims[2] - x_shift) / x_scale
  ymin_adj <- (-y_shift) / y_scale
  ymax_adj <- (img_dims[1] - y_shift) / y_scale
  
  ggplot() + annotation_custom(rasterGrob(img), xmin = xmin_adj, xmax = xmax_adj, ymin = ymin_adj, ymax = ymax_adj)
}

# Volcano Plot Function
volcano <- function(anot = "subset", ct, dif_col = "tissue", seu_obj, id1, id2) {
  
  object <- if(ct != "all") subset(seu_obj, subset = !!sym(anot) == ct) else seu_obj
  object <- NormalizeData(object)
  object <- ScaleData(object)
  object <- SetIdent(object, value = object@meta.data[[dif_col]])
  
  deg_results <- FindMarkers(object, ident.1 = id1, ident.2 = id2)
  deg_results <- na.omit(deg_results)
  deg_results$genes <- rownames(deg_results)
  
  deg_results$diffexpressed <- "NO"
  deg_results$diffexpressed[deg_results$avg_log2FC > log2(1.2) & deg_results$p_val < 0.05] <- "UP"
  deg_results$diffexpressed[deg_results$avg_log2FC < -log2(1.2) & deg_results$p_val < 0.05] <- "DOWN"
  deg_results$diffexpressed[deg_results$avg_log2FC > log2(1.2) & deg_results$p_val_adj < 0.05] <- "UPP"
  deg_results$diffexpressed[deg_results$avg_log2FC < -log2(1.2) & deg_results$p_val_adj < 0.05] <- "DOWNN"
  
  deg_results$delabel <- ifelse(deg_results$diffexpressed != "NO", deg_results$genes, NA)
  deg_results$p_val <- ifelse(deg_results$p_val < 1e-300, 1e-300, deg_results$p_val)
  
  genes_to_label <- deg_results$delabel[deg_results$diffexpressed %in% c("UPP", "DOWNN")]
  deg_results$delabel <- ifelse(deg_results$genes %in% genes_to_label, deg_results$genes, NA)
  
  diff_colors <- c("UPP" = "#803800", "DOWNN" = "#003F54", "UP" = "#B47846", "DOWN" = "steelblue")
  
  p <- ggplot(data = deg_results, aes(x = avg_log2FC, y = -log10(p_val), col = diffexpressed)) +
    geom_point(size = 2) +
    scale_color_manual(values = diff_colors) +
    theme_classic(base_size = 22) +
    guides(color = guide_legend(override.aes = list(shape = 1))) +
    theme(
      legend.position = "none",
      plot.title = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),
      axis.line = element_line(linewidth = 2), axis.ticks.length = unit(0.5, "cm"),
      axis.ticks = element_line(linewidth = 2), text = element_text(family = "Helvetica")
    ) +
    geom_label_repel(aes(label = delabel), size = 18 / .pt, segment.color = "black", segment.size = 2,
                     label.padding = unit(0.4, "lines"), label.size = 1, fontface = "bold", na.rm = TRUE)
  
  return(list(plot = p, data = deg_results))
}


# ==============================================================================
# Figure 3A: Spatial Map S+ and S- 
# ==============================================================================

# Slide 3 FOV 3 (S+)
alig_s3f3 <- align_image(patient = "Slide_3", fov = 3)
p_s3f3 <- plot.spatial(
  meta = meta, polygon_file = read.csv(file.path(POLY_DIR, "Slide_3.csv")),
  x = "x_local_px", y = "y_local_px", col_fov = "fov", fov = 3,
  palette = wide_pal, color = "wide", per = "black", alpha = 0,
  aligment = alig_s3f3, ptsize = 0.1
) 

png(filename = file.path(OUT_DIR, "fig3A_S3F3.png"), width = 10, height = 8, units = "in", res = 600)
print(p_s3f3)
dev.off()

# Slide 3 FOV 1 (S-)
alig_s3f1 <- align_image(patient = "Slide_3", fov = 1)
p_s3f1 <- plot.spatial(
  meta = meta, polygon_file = read.csv(file.path(POLY_DIR, "Slide_3.csv")),
  x = "x_local_px", y = "y_local_px", col_fov = "fov", fov = 1,
  palette = wide_pal, color = "wide", per = "black", alpha = 0,
  aligment = alig_s3f1, ptsize = 0.1
) 

png(filename = file.path(OUT_DIR, "fig3A_S3F1.png"), width = 10, height = 8, units = "in", res = 600)
print(p_s3f1)
dev.off()


# ==============================================================================
# Figure 3B: Volcano Plot Type S
# ==============================================================================

png(filename = file.path(OUT_DIR, "fig3B_volc_typeS.png"), width = 10, height = 7, units = "in", res = 600)
print(volcano(anot = "new_anot", ct = "Hepatocytes", dif_col = "Type_S", seu_obj = seu, id1 = "+", id2 = "-")[[1]])
dev.off()


# ==============================================================================
# Figure 3C: Spatial Map D+ vs D- (Slide 2 FOV 8)
# ==============================================================================

alig_s2f8 <- align_image(patient = "Slide_2", fov = 8)

# B-Spline Plot
pols_hdag <- read_csv(file.path(HDAG_DIR, "polygons_Slide_2_fov_8.csv"), show_col_types = FALSE)
pols_hdag$group <- ifelse(pols_hdag$polygon_id == 4, "-", "+")

p_spline <- alig_s2f8 + coord_fixed() +
  geom_bspline_closed(data = pols_hdag, aes(x = x, y = y, group = polygon_id, fill = group), color = "black", linewidth = 1, alpha = 0.15) +
  scale_fill_manual(values = final_pal) + theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.title = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), legend.position = "none")

png(filename = file.path(OUT_DIR, "fig3C_S2F8_alig.png"), width = 8, height = 8, units = "in", res = 800)
print(p_spline)
dev.off()

# Standard Spatial Plot
cell_pols_s2 <- read.csv(file.path(POLY_DIR, "Slide_2.csv"))
cells_pols_s2 <- cell_pols_s2[cell_pols_s2$fov == 8,]
meta_cut_s2 <- meta[meta$tissue == "Slide_2" & meta$fov == 8,]

type_d_vals <- meta_cut_s2$Type_D
names(type_d_vals) <- meta_cut_s2$cell_names
cells_pols_s2$type_d <- type_d_vals[cells_pols_s2$cell_names]

p_s2f8 <- plot.spatial(
  meta = meta, polygon_file = cell_pols_s2, x = "x_local_px", y = "y_local_px", 
  col_fov = "fov", fov = 8, palette = final_pal, color = "Type_D", 
  per = "black", alpha = 1, aligment = alig_s2f8, ptsize = 0.1
) 

png(filename = file.path(OUT_DIR, "fig3C_S2F8.png"), width = 10, height = 8, units = "in", res = 800)
print(p_s2f8)
dev.off()


# ==============================================================================
# Figure 3D: Volcano Plot Type D
# ==============================================================================

seu@meta.data$wide <- ifelse(seu@meta.data$subset == "Im_hepatocytes", "Hepatocytes", seu@meta.data$subset)

png(filename = file.path(OUT_DIR, "fig3D_volc_typeD.png"), width = 10, height = 7, units = "in", res = 800)
volcano_res <- volcano(anot = "new_anot", ct = "Hepatocytes", dif_col = "Type_D", seu_obj = seu, id1 = "+", id2 = "-")
my_plot <- volcano_res[[1]]
target_genes <- c("IFITM3", "IRF3", "CXCL9", "B2M", "VTN", "APOE", "CD74", 
                  "IGHG1", "STAT1", "CLU", "RARRES2", "IL32", "SOD2", "CRP", 
                  "NPR2", "ACKR4")
my_plot$data$delabel <- ifelse(my_plot$data$genes %in% target_genes, my_plot$data$genes, NA)
print(my_plot)
dev.off()

