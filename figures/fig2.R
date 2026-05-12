# ==============================================================================
# Script: fig2.R
# Description: Generates panels for Figure 2 (Spatial Polygons, Volcano, Pathways)
# ==============================================================================

# Load necessary libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(paletteer)
  library(plyr)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(circlize)
  library(clusterProfiler)
  library(msigdbr)
  library(org.Hs.eg.db)
})

options(bitmapType = "cairo")

# ==============================================================================
# Configuration & Paths
# Update these paths according to your local environment
# ==============================================================================
# Input Paths
SEURAT_PATH <- "/path/to/Objects/seurats_annotated.RDS"
POLY_DIR    <- "/path/to/Polygons/"

# Output Paths (Relative to the repository root)
OUT_DIR <- "figures/outs/"

# Ensure output directory exists
if (!dir.exists(OUT_DIR)) {
  dir.create(OUT_DIR, recursive = TRUE)
}

# ==============================================================================
# 1. Load Pre-Annotated Seurat Object & Define Palettes
# ==============================================================================
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

# ==============================================================================
# FIG 2A & 2B (Pending code from Angela)
# ==============================================================================

# ==============================================================================
# FIG 2C - Spatial Polygons
# ==============================================================================

# Helper function to prevent repetitive plotting code
plot_spatial_fov <- function(slide_csv, target_slide, target_fov, meta_df, color_col, palette, out_file) {
  
  # Load and filter polygons
  pols <- read_csv(slide_csv, show_col_types = FALSE)
  pols <- pols[pols$fov == target_fov, ]
  
  # Filter metadata
  meta_cut <- meta_df[meta_df$tissue == target_slide & meta_df$fov == target_fov, ]
  pols <- pols[pols$cell_names %in% meta_cut$cell_names, ]
  
  # Map metadata to polygons
  pols[[color_col]] <- plyr::mapvalues(
    x = pols$cell_names, 
    from = meta_cut$cell_names, 
    to = meta_cut[[color_col]]
  )
  pols[[color_col]] <- as.factor(pols[[color_col]])
  
  # Generate Plot
  p <- ggplot(pols, aes(x = .data[["x_global_px"]], y = .data[["y_global_px"]])) +
    geom_polygon(aes(group = .data[["cell_names"]], fill = .data[[color_col]]), color = "#000000") +
    ggdark::dark_theme_gray(base_family = "Fira Sans Condensed Light", base_size = 20) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.background = element_blank()
    ) + 
    scale_fill_manual(values = palette) + 
    guides(fill = "none")
  
  # Save Plot
  png(filename = out_file, width = 8, height = 8, units = "in", res = 1200)
  print(p)
  dev.off()
}

# --- Plot Execution ---

# 1. Wide Annotations
plot_spatial_fov(file.path(POLY_DIR, "Slide_2.csv"), "Slide_2", 2, meta, "wide", wide_col, file.path(OUT_DIR, "cslide2_fov2_wide.png"))
plot_spatial_fov(file.path(POLY_DIR, "Slide_1.csv"), "Slide_1", 1, meta, "wide", wide_col, file.path(OUT_DIR, "cslide1_fov1_wide.png"))
plot_spatial_fov(file.path(POLY_DIR, "Slide_2.csv"), "Slide_2", 16, meta, "wide", wide_col, file.path(OUT_DIR, "cslide2_fov16_wide.png"))

# 2. Refined Annotations (Hepatocytes Only)
hep_types <- c("Hepatocyte_1", "Hepatocyte_2", "Hepatocyte_3", "Hepatocyte_4", "Hepatocyte_5", "Hepatocyte_6")
meta$refined_hep <- ifelse(meta$refined %in% hep_types, meta$refined, "Other")

plot_spatial_fov(file.path(POLY_DIR, "Slide_2.csv"), "Slide_2", 2, meta, "refined_hep", refined_col, file.path(OUT_DIR, "cslide2_fov2_refined.png"))
plot_spatial_fov(file.path(POLY_DIR, "Slide_2.csv"), "Slide_2", 16, meta, "refined_hep", refined_col, file.path(OUT_DIR, "cslide2_fov16_refined.png"))

# 3. Refined Annotations (Myeloids Only)
myeloid_types <- c("M2_LYVE1", "KC1", "KC2", "Monocytes", "DCs_CD1C", "M1")
meta$refined_myel <- ifelse(meta$refined %in% myeloid_types, meta$refined, "Other")

plot_spatial_fov(file.path(POLY_DIR, "Slide_1.csv"), "Slide_1", 1, meta, "refined_myel", refined_col, file.path(OUT_DIR, "cslide1_fov1_refined_myel.png"))
plot_spatial_fov(file.path(POLY_DIR, "Slide_2.csv"), "Slide_2", 2, meta, "refined_myel", refined_col, file.path(OUT_DIR, "cslide2_fov2_refined_myel.png"))

# ==============================================================================
# FIG 2D (Pending code from Angela)
# ==============================================================================

# ==============================================================================
# FIG 2E - Volcano & Enrichment Analysis
# ==============================================================================

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
  deg_results$diffexpressed[deg_results$avg_log2FC > log2(1.2) & deg_results$p_val < 0.05] <- "p.val<0.05 & FC>1.2"
  deg_results$diffexpressed[deg_results$avg_log2FC < -log2(1.2) & deg_results$p_val < 0.05] <- "p.val<0.05 & FC<0.83"
  deg_results$diffexpressed[deg_results$avg_log2FC > log2(1.2) & deg_results$p_val_adj < 0.05] <- "p.adj<0.05 & FC>1.2"
  deg_results$diffexpressed[deg_results$avg_log2FC < -log2(1.2) & deg_results$p_val_adj < 0.05] <- "p.adj<0.05 & FC<0.83"
  
  deg_results$delabel <- ifelse(deg_results$diffexpressed != "NO", deg_results$genes, NA)
  deg_results$p_val <- ifelse(deg_results$p_val < 1e-300, 1e-300, deg_results$p_val)
  
  # Plot
  p <- ggplot(data = deg_results, aes(x = avg_log2FC, y = -log10(p_val), col = diffexpressed, label = delabel)) + 
    geom_point() + 
    theme_bw(base_size = 18) +
    geom_text_repel() +
    scale_color_manual(values = c("p.val<0.05 & FC<0.83" = "green", "p.adj<0.05 & FC<0.83" = "darkgreen", 
                                  "p.val<0.05 & FC>1.2" = "red", "p.adj<0.05 & FC>1.2" = "darkred")) +
    geom_vline(xintercept = c(-log2(1.2), log2(1.2)), col = "black", linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), col = "black", linetype = "dashed") +
    ggtitle(paste0(ct, ": ", id1, " vs ", id2))
  
  return(list(plot = p, data = deg_results))
}

# --- Execute Volcano ---
hc_vs_hbv_res <- volcano(anot = "new_anot", ct = "Hepatocytes", seu_obj = seu, dif_col = "etiology", id1 = "HDV RNA+", id2 = "HBV")
hc_vs_hbv <- hc_vs_hbv_res$data

# --- Pathway Enrichment ---
msigdb_hallmark <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
gene_sets <- msigdb_hallmark %>% dplyr::select(gs_name, entrez_gene)

# Process Upregulated
up_genes <- hc_vs_hbv$genes[hc_vs_hbv$diffexpressed %in% c("p.adj<0.05 & FC>1.2", "p.val<0.05 & FC>1.2")]
gene_list_up <- bitr(up_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

up_enrichment <- as.data.frame(enricher(gene = gene_list_up$ENTREZID, TERM2GENE = gene_sets))
up_enrichment <- up_enrichment %>% 
  mutate(GeneRatio = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
         s1 = "Upregulated")

# Process Downregulated
down_genes <- hc_vs_hbv$genes[hc_vs_hbv$diffexpressed %in% c("p.adj<0.05 & FC<0.83", "p.val<0.05 & FC<0.83")]
gene_list_down <- bitr(down_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

down_enrichment <- as.data.frame(enricher(gene = gene_list_down$ENTREZID, TERM2GENE = gene_sets))
down_enrichment <- down_enrichment %>% 
  mutate(GeneRatio = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
         s1 = "Downregulated")

# Target Pathways
pathways <- c(
  "GOBP_VIRAL_PROCESS",
  "GOBP_RESPONSE_TO_VIRUS",
  "GOBP_INTERFERON_MEDIATED_SIGNALING_PATHWAY",
  "GOBP_RESPONSE_TO_TYPE_I_INTERFERON",
  "GOBP_RESPONSE_TO_TYPE_II_INTERFERON",
  "GOBP_RESPONSE_TO_TYPE_III_INTERFERON",
  "GOBP_INTERLEUKIN_6_MEDIATED_SIGNALING_PATHWAY",
  "GOBP_INTERLEUKIN_6_PRODUCTION",
  "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION",
  "GOBP_T_CELL_ACTIVATION",
  "GOBP_T_CELL_PROLIFERATION",
  "GOBP_CYTOKINE_PRODUCTION",
  "GOBP_TRANSFORMING_GROWTH_FACTOR_BETA2_PRODUCTION",
  "GOBP_MACROPHAGE_ACTIVATION",
  "GOBP_ERK1_AND_ERK2_CASCADE",
  "GOBP_MAPK_CASCADE",
  "GOBP_REGULATION_OF_CYTOPLASMIC_PATTERN_RECOGNITION_RECEPTOR_SIGNALING_PATHWAY",
  "GOBP_CELL_CYCLE_G2_M_PHASE_TRANSITION",
  "GOBP_PHOSPHORYLATION",
  "GOBP_RESPONSE_TO_TUMOR_NECROSIS_FACTOR",
  "GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION",
  "GOBP_RESPONSE_TO_OXYGEN_LEVELS",
  "GOBP_RESPONSE_TO_OXIDATIVE_STRESS",
  "GOBP_CELL_MATRIX_ADHESION",
  "GOBP_CELL_CELL_ADHESION",
  "GOBP_TISSUE_REMODELING",
  "GOBP_TISSUE_REGENERATION",
  "GOBP_LIPID_METABOLIC_PROCESS",
  "GOBP_FATTY_ACID_BETA_OXIDATION",
  "GOBP_REGULATION_OF_GLUCONEOGENESIS",
  "GOBP_XENOBIOTIC_METABOLIC_PROCESS"
)

# Combine and Filter Pathways
final_enrichment <- bind_rows(up_enrichment, down_enrichment) %>%
  filter(ID %in% pathways) %>%
  arrange(GeneRatio)

# --- Plot Enrichment ---
p_path <- ggplot(final_enrichment, aes(x = s1, y = Description, size = GeneRatio, color = s1)) +
  geom_point() +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +
  scale_size_continuous(name = "Gene ratio", range = c(3, 10), guide = guide_legend(override.aes = list(color = "black", fill = "black"))) +
  theme_bw(base_size = 15) +
  labs(x = "", y = "Pathway description", color = "Pathway") +
  theme(
    axis.line = element_line(linewidth = 1.2, colour = "black"),
    axis.ticks = element_line(linewidth = 1.2, colour = "black"),
    axis.ticks.length = unit(6, "pt"),
    panel.background = element_rect(fill = NA, colour = NA),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2)
  )

png(filename = file.path(OUT_DIR, "pathway_enrichment.png"), width = 11, height = 8, units = "in", res = 1200)
print(p_path)
dev.off()

write.csv(final_enrichment, file.path(OUT_DIR, "pathway_cosmx.csv"), row.names = FALSE)