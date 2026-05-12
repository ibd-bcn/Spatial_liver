# ==============================================================================
# Script: 5.pseudotime.R
# Description: Slingshot trajectory inference, DimReduc, and tradeSeq analysis
# ==============================================================================

# ==============================================================================
# 0. Environment & Dependencies
# ==============================================================================
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

pacman::p_load(
  # Data wrangling
  dplyr, tidyr, purrr, stringr, readr, readxl, jsonlite, plyr,
  # Seurat & SingleCellExperiment
  Seurat, SeuratObject, SingleCellExperiment,
  # Trajectory & DiffExp
  slingshot, tradeSeq,
  # Spatial & plotting
  ggplot2, viridis, ggrepel, patchwork, paletteer, ggdark, ComplexHeatmap,
  circlize, visNetwork, DT, grDevices, RColorBrewer,
  # Enrichment
  clusterProfiler, msigdbr, org.Hs.eg.db,
  # Misc
  Matrix, BiocParallel, progress
)
options(stringsAsFactors = FALSE)
set.seed(20250714)

# ==============================================================================
# Configuration & Paths
# Update these paths according to your local environment
# ==============================================================================
SEURAT_IN_PATH <- "path/to/Post_analysis/Selecting_regions/seu.RDS"
UMAP_CSV_PATH  <- "path/to/Post_analysis/Neighborhood/Run/30/umap.csv"
NEIGH_CSV_PATH <- "path/to/Post_analysis/Neighborhood/Run/30/seu_neigh.csv"

OUT_DIR        <- "path/to/Post_analysis/Post_neigh/"
OUT_PLOT_DIR   <- file.path(OUT_DIR, "Plots")

# Ensure output directories exist
invisible(lapply(c(OUT_DIR, OUT_PLOT_DIR), dir.create, recursive = TRUE, showWarnings = FALSE))

# Global Analysis Parameters
TARGET_TRAJECTORY <- 5
N_CORES           <- 8 # Used for tradeSeq fitGAM parallelization

# ==============================================================================
# 1. Annotation Dictionaries
# ==============================================================================
subset_type <- c(
  "hepatocytes"       = "Hepatocytes",
  "myeloids"          = "Myeloid_cells",
  "plasmas"           = "B_cell_lineage",
  "tcells"            = "T_cell_lineage",
  "non-parenquimal"   = "Non_parenchymal",
  "IM-hepatocytes"    = "Im_hepatocytes"
)

grouped_anot <- c(
  "Non_parenchymal"   = "Hepatocytes_Non_Parenchymal",
  "Hepatocytes"       = "Hepatocytes_Non_Parenchymal",
  "Myeloid_cells"     = "Immune_cells",
  "B_cell_lineage"    = "Immune_cells",
  "T_cell_lineage"    = "Immune_cells",
  "Im_hepatocytes"    = "Hepatocytes_Non_Parenchymal"
)

cell_types <- c(
  "Gamma-delta T cells"       = "Gamma_delta_T_cells",
  "NK cells"                  = "NK_cells",
  "NKT cells"                 = "NKT_cells",
  "Regulatory T cells"        = "Regulatory_T_cells",
  "Tem/Trm cytotoxic T cells" = "Tem_Trm_cytotoxic_T_cells",
  "Effector helper T cells"   = "Effector_helper_T_cells",
  "Naive T cells"             = "Naive_T_cells",
  "Memory B cells"            = "Memory_B_cells",
  "Plasma cells"              = "Plasma_cells",
  "Naive B cells"             = "Naive_B_cells",
  "Monocytes"                 = "Monocytes",
  "KC1"                       = "KC1",
  "M2-LYVE1"                  = "M2_LYVE1",
  "KC2"                       = "KC2",
  "DCs CD1C"                  = "DCs_CD1C",
  "M1"                        = "M1",
  "Endothelial cells 2"       = "Endothelial_cells_2",
  "Fibroblasts"               = "Fibroblasts",
  "Endothelial cells 1"       = "Endothelial_cells_1",
  "Smooth muscle cells"       = "Smooth_muscle_cells",
  "Endothelial cells 4"       = "Endothelial_cells_4",
  "Endothelial cells 3"       = "Endothelial_cells_3",
  "Hepatocyte 2"              = "Hepatocyte_2",
  "Hepatocyte 1"              = "Hepatocyte_1",
  "Hepatocyte 6"              = "Hepatocyte_6",
  "Hepatocyte 3"              = "Hepatocyte_3",
  "Hepatocyte 5"              = "Hepatocyte_5",
  "Hepatocyte 4"              = "Hepatocyte_4",
  "Cholangiocytes"            = "Cholangiocytes",
  "Cycling B lineage cells"   = "Cycling_B_lineage_cells"
)

new_anot <- c(
  "Gamma_delta_T_cells"       = "T_cells",
  "NK_cells"                  = "T_cells",
  "NKT_cells"                 = "T_cells",
  "Regulatory_T_cells"        = "T_cells",
  "Tem_Trm_cytotoxic_T_cells" = "T_cells",
  "Effector_helper_T_cells"   = "T_cells",
  "Naive_T_cells"             = "T_cells",
  "Memory_B_cells"            = "B_cells",
  "Plasma_cells"              = "B_cells",
  "Naive_B_cells"             = "B_cells",
  "Monocytes"                 = "Myeloid_cells",
  "KC1"                       = "Kupffer_cells",
  "M2_LYVE1"                  = "Myeloid_cells",
  "KC2"                       = "Kupffer_cells",
  "DCs_CD1C"                  = "Myeloid_cells",
  "M1"                        = "Myeloid_cells",
  "Endothelial_cells_2"       = "Non_parenchymal",
  "Fibroblasts"               = "Non_parenchymal",
  "Endothelial_cells_1"       = "Non_parenchymal",
  "Smooth_muscle_cells"       = "Non_parenchymal",
  "Endothelial_cells_4"       = "Non_parenchymal",
  "Endothelial_cells_3"       = "Non_parenchymal",
  "Hepatocyte_2"              = "Hepatocytes",
  "Hepatocyte_1"              = "Hepatocytes",
  "Hepatocyte_6"              = "Hepatocytes",
  "Hepatocyte_3"              = "Hepatocytes",
  "Hepatocyte_5"              = "Hepatocytes",
  "Hepatocyte_4"              = "Hepatocytes",
  "Cholangiocytes"            = "Non_parenchymal",
  "Cycling_B_lineage_cells"   = "B_cells"
)

# ==============================================================================
# 2. Load Objects & Join Layers
# ==============================================================================
stopifnot("Seurat path not found!" = file.exists(SEURAT_IN_PATH))
seu <- readRDS(SEURAT_IN_PATH) 

# Ensure custom Spatial_Package function works, otherwise standard JoinLayers
if(exists("JoinLayers", where="package:Seurat")) {
  seu <- JoinLayers(seu)
} else {
  seu <- Spatial_Package::JoinLayers(seu)
}

# ==============================================================================
# 3. Metadata Enrichment & Recoding
# ==============================================================================
meta <- seu@meta.data

meta <- meta %>%
  mutate(
    # 3.1 Mark specific slide FOVs
    red = (tissue == "Slide_2" & fov %in% c(21, 22, 23, 24)),
    orange = case_when(
      tissue == "Slide_1" & fov == 25                            ~ TRUE,
      tissue == "Slide_2" & fov %in% c(2, 3, 4, 6, 7, 10)        ~ TRUE,
      tissue == "Slide_3" & fov %in% c(6, 12, 14, 23)            ~ TRUE,
      TRUE                                                       ~ FALSE
    ),
    
    # 3.2 Assign Antigen Levels
    antigen_level = case_when(
      tissue == "Slide_1" & fov %in% c(1:8, 25)                  ~ "neg",
      tissue == "Slide_3" & fov %in% c(1, 4, 11, 16, 24, 25)     ~ "neg",
      tissue == "Slide_2" & fov %in% c(13, 14, 15, 17, 19, 20)   ~ "S_low",
      tissue == "Slide_3" & fov == 13                            ~ "S_low",
      tissue == "Slide_2" & fov %in% c(16, 18)                   ~ "S_high",
      tissue == "Slide_3" & fov %in% c(2, 3, 5, 6, 12, 14, 23)   ~ "S_high",
      tissue == "Slide_3" & fov %in% c(7, 8, 9, 10, 15, 17:22)   ~ "S+/-",
      tissue == "Slide_2" & fov %in% c(1, 5, 6, 7, 9, 11, 12)    ~ "S&D_low",
      tissue == "Slide_2" & fov %in% c(2, 3, 4, 8, 10)           ~ "S&D",
      TRUE                                                       ~ NA_character_
    ),
    
    # [Insert Etiology & Patient_Type definitions here]
    
    # 3.3 Vectorised Recoding
    subset       = recode(subset, !!!subset_type),
    refined      = recode(refined, !!!cell_types),
    new_anot     = recode(refined, !!!new_anot),
    grouped_anot = recode(subset,  !!!grouped_anot)
  )

# ==============================================================================
# 4. Integrate UMAP & Leiden Labels
# ==============================================================================
umap_df   <- read_csv(UMAP_CSV_PATH, show_col_types = FALSE)
seu_neigh <- read_csv(NEIGH_CSV_PATH, show_col_types = FALSE)
cell_names <- rownames(meta)

# Use standard matching (Assuming the first column contains cell names)
meta$umap1     <- umap_df$UMAP1[match(cell_names, umap_df[[1]])]
meta$umap2     <- umap_df$UMAP2[match(cell_names, umap_df[[1]])]
meta$leiden0.2 <- seu_neigh$leiden_res0.2[match(cell_names, seu_neigh$cell_names)]

seu@meta.data <- meta

# ==============================================================================
# 5. DimReduc & Normalization
# ==============================================================================
umap_coords <- as.matrix(meta[, c("umap1", "umap2")])
rownames(umap_coords) <- rownames(meta)

seu[["umap"]] <- CreateDimReducObject(
  embeddings = umap_coords,
  key        = "UMAP_",
  assay      = DefaultAssay(seu)
)

seu <- NormalizeData(seu)
seu <- ScaleData(seu)

# ==============================================================================
# 6. Slingshot Trajectory Inference
# ==============================================================================
sce <- as.SingleCellExperiment(seu)

sce <- slingshot(
  sce,
  clusterLabels = "leiden0.2",
  reducedDim    = "UMAP"
)

saveRDS(sce, file = file.path(OUT_DIR, "sce_slingshot.rds"))
message("Slingshot completed -> Output saved as 'sce_slingshot.rds'")

# ==============================================================================
# 7. Visualization & Plotting
# ==============================================================================
# Save Base R Plots
pdf(file.path(OUT_PLOT_DIR, "Slingshot_Base_Plots.pdf"), width=8, height=6)

# 7.1 Pseudotime Coloring
colors <- colorRampPalette(brewer.pal(11, 'Spectral')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]
plot(reducedDims(sce)$UMAP, col = plotcol, pch=16, asp = 1, main="Slingshot Pseudotime")
lines(SlingshotDataSet(sce), lwd=2, col='black')

# 7.2 Cluster Coloring
plot(reducedDims(sce)$UMAP, col = brewer.pal(9, 'Set1')[as.factor(sce$leiden0.2)], pch=16, asp = 1, main="Slingshot by Leiden Clusters")
lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black')

dev.off()

# 7.3 Ggplot Trajectory Output
plot_df <- data.frame(
  UMAP_1     = reducedDim(sce, 'UMAP')[, 1],
  UMAP_2     = reducedDim(sce, 'UMAP')[, 2],
  pseudotime = slingPseudotime(sce)[, TARGET_TRAJECTORY]
)

curves <- slingCurves(sce)
curve1_df <- as.data.frame(curves[[TARGET_TRAJECTORY]]$s)
colnames(curve1_df) <- c("UMAP_1", "UMAP_2")

umap_plot_with_line <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = pseudotime), size = 0.01, alpha = 0.8) +
  geom_path(data = curve1_df, aes(x = UMAP_1, y = UMAP_2), linewidth = 1.2, color = "black") +
  scale_color_viridis(option = "magma", na.value = "lightgrey") +
  labs(
    title = "Slingshot Pseudotime on UMAP",
    subtitle = paste0("Trajectory ", TARGET_TRAJECTORY, " with Principal Curve"),
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

ggsave(filename = file.path(OUT_PLOT_DIR, paste0("Slingshot_Trajectory_", TARGET_TRAJECTORY, ".png")), 
       plot = umap_plot_with_line, width = 8, height = 6, dpi = 300)

# ==============================================================================
# 8. TradeSeq Dynamic Expression Analysis
# ==============================================================================
message("Running tradeSeq fitGAM (This may take a while depending on cell count...)")

# Use BiocParallel to drastically speed up GAM fitting
bp_param <- MulticoreParam(workers = N_CORES, progressbar = TRUE)

sce <- fitGAM(sce, BPPARAM = bp_param)
ATres <- associationTest(sce, lineages = TRUE)

saveRDS(ATres, file = file.path(OUT_DIR, "ATres.RDS"))
saveRDS(sce, file = file.path(OUT_DIR, "sce_tradeseq.rds"))

message("tradeSeq completed -> Output saved as 'sce_tradeseq.rds' & 'ATres.RDS'")