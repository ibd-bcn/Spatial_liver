# ==============================================================================
# Script: fig7.R
# Description: Generates panels for Figure 7 (Spatial Maps, UMAPs, Trajectories)
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(slingshot)
})

options(bitmapType = "cairo")

# ==============================================================================
# Configuration & Paths
# Update these paths according to your local environment
# ==============================================================================
SEURAT_PATH    <- "/path/to/Objects/seurats_annotated.RDS"
POLY_DIR       <- "/path/to/Polygons/"
NEIGH_CSV_PATH <- "/path/to/Neighborhood/seu_neigh.csv"
UMAP_CSV_PATH  <- "/path/to/Neighborhood/umap.csv"
SCE_PATH       <- "/path/to/Post_neigh/sce_slingshot.rds"

OUT_DIR <- "figures/outs/"

if (!dir.exists(OUT_DIR)) {
  dir.create(OUT_DIR, recursive = TRUE)
}

# ==============================================================================
# 1. Load Data & Define Palettes
# ==============================================================================
seu <- readRDS(SEURAT_PATH)
meta <- seu@meta.data

seu_neigh <- read_csv(NEIGH_CSV_PATH, show_col_types = FALSE)
umap_csv  <- read_csv(UMAP_CSV_PATH, show_col_types = FALSE)
colnames(umap_csv)[1] <- "cell_names"

# Map 'wide' annotations
meta$wide <- ifelse(meta$subset == "Im_hepatocytes", "Hepatocytes", meta$subset)

wide_pal <- c(
  "Hepatocytes"     = "#ED9824",
  "Myeloid_cells"   = "#5387C9",
  "B_cell_lineage"  = "#772D8B",
  "T_cell_lineage"  = "#A8C686",
  "Non_parenchymal" = "#F05365"
)

neigh_cols <- c(
  "0" = "#CE3D32FF", "1" = "#F24B91",   "2" = "#BA6338FF", "3" = "#F0E685FF",
  "4" = "#6BD76BFF", "5" = "#749B58FF", "6" = "#5DB1DDFF", "7" = "#5050FFFF",
  "8" = "#466983FF", "9" = "#802268FF", "10"= "#F88A4C",   "11"= "#D5A5DA"
)

# ==============================================================================
# 2. Helper Plotting Functions
# ==============================================================================
plot_polygon_fov <- function(slide_id, fov_id, meta_df, color_col, palette, out_file) {
  poly_file <- file.path(POLY_DIR, paste0(slide_id, ".csv"))
  pols <- read_csv(poly_file, show_col_types = FALSE) %>% filter(fov == fov_id)
  
  meta_cut <- meta_df %>% filter(tissue == slide_id, fov == fov_id)
  pols <- pols %>% filter(cell_names %in% meta_cut$cell_names)
  
  pols[[color_col]] <- plyr::mapvalues(pols$cell_names, from = meta_cut$cell_names, to = meta_cut[[color_col]], warn_missing = FALSE)
  pols[[color_col]] <- as.factor(pols[[color_col]])
  
  p <- ggplot(pols, aes(x = x_global_px, y = y_global_px)) +
    geom_polygon(aes(group = cell_names, fill = .data[[color_col]]), color = "#000000") +
    ggdark::dark_theme_gray(base_family = "Fira Sans Condensed Light", base_size = 20) +
    theme(
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      axis.title = element_blank(), axis.text = element_blank(),
      axis.ticks = element_blank(), panel.background = element_blank()
    ) + 
    scale_fill_manual(values = palette) + guides(fill = "none")
  
  png(filename = out_file, width = 8, height = 8, units = "in", res = 1200)
  print(p)
  dev.off()
}

plot_point_fov <- function(slide_id, fov_id, meta_df, color_col, palette, out_file) {
  meta_cut <- meta_df %>% filter(tissue == slide_id, fov == fov_id)
  meta_cut[[color_col]] <- as.factor(meta_cut[[color_col]])
  
  p <- ggplot(meta_cut, aes(x = CenterX_global_px, y = CenterY_global_px, color = .data[[color_col]])) +
    geom_point(size = 3) +
    ggdark::dark_theme_gray(base_family = "Fira Sans Condensed Light", base_size = 20) +
    theme(
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      axis.title = element_blank(), axis.text = element_blank(),
      axis.ticks = element_blank(), panel.background = element_blank()
    ) + 
    scale_color_manual(values = palette) + guides(color = "none")
  
  png(filename = out_file, width = 8, height = 8, units = "in", res = 1200)
  print(p)
  dev.off()
}

# ==============================================================================
# Figure 7A: Spatial Maps (Wide & Neighborhood)
# ==============================================================================
# Wide Annotations
plot_polygon_fov("Slide_3", 9, meta, "wide", wide_pal, file.path(OUT_DIR, "fig7A_wide_S3F9.png"))
plot_polygon_fov("Slide_2", 7, meta, "wide", wide_pal, file.path(OUT_DIR, "fig7A_wide_S2F7.png"))
plot_polygon_fov("Slide_1", 25, meta, "wide", wide_pal, file.path(OUT_DIR, "fig7A_wide_S1F25.png"))
plot_point_fov("Slide_healthy", 4, meta, "wide", wide_pal, file.path(OUT_DIR, "fig7A_wide_SHCF4.png"))

# Neighborhood Annotations
seu_neigh$leiden_res0.2 <- as.factor(seu_neigh$leiden_res0.2)
plot_polygon_fov("Slide_1", 25, seu_neigh, "leiden_res0.2", neigh_cols, file.path(OUT_DIR, "fig7A_neigh_S1F25.png"))
plot_polygon_fov("Slide_3", 9, seu_neigh, "leiden_res0.2", neigh_cols, file.path(OUT_DIR, "fig7A_neigh_S3F9.png"))
plot_polygon_fov("Slide_2", 7, seu_neigh, "leiden_res0.2", neigh_cols, file.path(OUT_DIR, "fig7A_neigh_S2F7.png"))
plot_point_fov("Slide_healthy", 4, seu_neigh, "leiden_res0.2", neigh_cols, file.path(OUT_DIR, "fig7A_neigh_SHCF4.png"))

# ==============================================================================
# Figure 7B: Overall Proportions Barplot
# ==============================================================================
df_bar <- seu_neigh %>%
  mutate(
    broad_group = case_when(
      leiden_res0.2 %in% c(11, 6, 2) ~ "healthy/recovery",
      leiden_res0.2 %in% c(0, 3, 4, 7, 8, 10, 5) ~ "immune_active",
      leiden_res0.2 %in% c(1, 9) ~ "anti-inflamatory",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(broad_group)) %>%
  group_by(etiology, broad_group) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  group_by(etiology) %>%
  mutate(total_cells = sum(cell_count), proportion = cell_count / total_cells) %>%
  ungroup()

df_bar$etiology <- factor(df_bar$etiology, levels = c("HC", "HBV", "HDV RNA+", "HDV RNA-"))
df_bar$broad_group <- factor(df_bar$broad_group, levels = c("healthy/recovery", "immune_active", "anti-inflamatory"))

group_colors <- c("healthy/recovery" = "#4daf4a", "immune_active" = "#e41a1c", "anti-inflamatory" = "#377eb8")

p_7b <- ggplot(df_bar, aes(x = etiology, y = proportion, fill = broad_group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, color = "black") +
  scale_fill_manual(values = group_colors) + theme_classic() + 
  labs(x = "Etiology", y = "Overall Proportion of Cells", fill = "Cluster Group") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), text = element_text(size = 25),
    axis.line = element_line(linewidth = 1.5), axis.ticks = element_line(linewidth = 1.2),
    axis.ticks.length = unit(0.3, "cm")
  )

png(filename = file.path(OUT_DIR, "fig7B_barplot_overall_proportions.png"), width = 9.2, height = 6, units = "in", res = 1200)
print(p_7b)
dev.off()

# ==============================================================================
# Figure 7C: UMAP Split by Etiology
# ==============================================================================
umap_df <- umap_csv %>%
  left_join(seu_neigh[, c("cell_names", "leiden_res0.2")], by = "cell_names") %>%
  left_join(meta[, c("cell_names", "etiology")], by = "cell_names") %>%
  mutate(etiology = factor(etiology, levels = c("HC", "HBV", "HDV RNA+", "HDV RNA-")))

x_limits <- range(umap_df$UMAP1, na.rm = TRUE)
y_limits <- range(umap_df$UMAP2, na.rm = TRUE)

for (eti in unique(umap_df$etiology)) {
  subset_umap <- filter(umap_df, etiology == eti)
  p_7c <- ggplot(subset_umap, aes(x = UMAP1, y = UMAP2, color = as.factor(leiden_res0.2))) +
    geom_point(size = 0.5) + scale_color_manual(values = neigh_cols) + 
    coord_cartesian(xlim = x_limits, ylim = y_limits) + theme_void() + theme(legend.position = "none")
  
  png(filename = file.path(OUT_DIR, paste0("fig7C_UMAP_etiology_", eti, ".png")), width = 8, height = 8, units = "in", res = 1200)
  print(p_7c)
  dev.off()
}

# ==============================================================================
# Figure 7D: Slingshot Trajectories
# ==============================================================================
sce_slingshot <- readRDS(SCE_PATH)

for (traj in c(5, 6)) {
  plot_df <- data.frame(
    UMAP_1 = reducedDim(sce_slingshot, 'UMAP')[, 1],
    UMAP_2 = reducedDim(sce_slingshot, 'UMAP')[, 2],
    pseudotime = slingPseudotime(sce_slingshot)[, traj]
  )
  
  curves <- slingCurves(sce_slingshot)
  curve_df <- as.data.frame(curves[[traj]]$s)
  colnames(curve_df) <- c("UMAP_1", "UMAP_2")
  
  p_7d <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = pseudotime), size = 0.01, alpha = 0.8) +
    geom_path(data = curve_df, aes(x = UMAP_1, y = UMAP_2), linewidth = 1.2, color = "black") +
    scale_color_viridis_c(option = "magma", na.value = "lightgrey") + theme_minimal(base_size = 14) +
    labs(title = "Slingshot Pseudotime on UMAP", subtitle = paste0("Trajectory ", traj), color = "Pseudotime") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5),
      panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()
    )
  
  png(filename = file.path(OUT_DIR, paste0("fig7D_Slingshot_Traj", traj, ".png")), width = 9.2, height = 6, units = "in", res = 800)
  print(p_7d)
  dev.off()
}

# ==============================================================================
# Figure 7E: Module Scores (Seurat Expression)
# ==============================================================================
umap_ordered <- umap_df[match(colnames(seu), umap_df$cell_names), ]
umap_matrix <- as.matrix(umap_ordered[, c("UMAP1", "UMAP2")])
rownames(umap_matrix) <- colnames(seu)
colnames(umap_matrix) <- c("UMAP_1", "UMAP_2") 

seu[["umap"]] <- CreateDimReducObject(embeddings = umap_matrix, key = "UMAP_", assay = DefaultAssay(seu))
seu <- NormalizeData(seu)

signatures <- list(
  Healthy_parenchyma = c("APOA1", "APOC1", "APOE", "TTR", "HPGDS", "FABP4", "CYP1B1", "INSR", "IGF1", "IGF2", "GLUL", "SLC2A1", "HSD17B2"),
  Immune_activation  = c("B2M", "TAP1", "TAP2", "HLA-DRA", "HLA-DPA1", "CIITA", "STAT1", "IFIT3", "IFIT1", "ISG15", "MX1", "OAS1"),
  Antigen_presenting = c("CSF1R", "CLEC4A", "CD14", "CD1C", "TLR2", "IL15RA", "NOD2", "IDO1", "HLA-DPA1", "CXCL17", "CCL20"),
  Tissue_remodelling = c("ANGPT2", "COL16A1", "COL15A1", "TWIST1", "MMP7", "TIMP1", "TAGLN", "NR1H2", "NR2F2", "CXCL17"),
  Immune_tolerance   = c("FOXP3", "IL10", "TGFB1", "CTLA4", "IL2RA", "ICOS", "CD274", "HAVCR2", "GATA3", "LGALS9", "STAT5B", "AHR", "TGFBR2", "IL1RN", "C1QA"),
  Immune_exhaustion  = c("PDCD1", "LAG3", "TIGIT", "TOX", "EOMES", "CXCL13", "CD38", "HIF1A", "BCL2", "GZMK", "IL7R", "STAT3", "NR3C1", "IRF4"),
  Immune_suppression = c("ARG1", "IDO1", "IL10RA", "IL10RB", "TGFB2", "TGFBR1", "MERTK", "CSF1R", "CD163", "MRC1", "VCAN", "CSTB", "S100A9", "S100A8", "LGALS3")
)

for (sig_name in names(signatures)) {
  seu <- AddModuleScore(object = seu, features = list(signatures[[sig_name]]), name = sig_name, ctrl = 25)
  target_feature <- paste0(sig_name, "1")
  
  p_7e <- FeaturePlot(seu, features = target_feature, pt.size = 0.05, order = TRUE, raster = FALSE, min.cutoff = "q80", max.cutoff = "q99") + 
    scale_colour_gradient(low = "grey90", high = "#FF0033") + 
    theme(aspect.ratio = 1, axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) +
    ggtitle(gsub("_", " ", sig_name)) 
  
  png(filename = file.path(OUT_DIR, paste0("fig7E_UMAP_", sig_name, ".png")), width = 9.2, height = 6, units = "in", res = 1200)
  print(p_7e)
  dev.off()
}