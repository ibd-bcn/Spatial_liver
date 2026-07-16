# ==============================================================================
# Script: supfig2.R
# Description: Generates panels for Supplementary Figure 2 (UMAPs, Barplots, Volcano, GO)
# ==============================================================================

# Load necessary libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(grid)
})

options(bitmapType = "cairo")

# ==============================================================================
# Configuration & Paths
# Update these paths according to your local environment
# ==============================================================================
# Spatial Object Path
SEURAT_PATH <- "/path/to/Objects/seurats_annotated.RDS"

# scRNAseq Object Paths
ALL_SEURAT_PATH         <- "/path/to/seurat_objects/todas.rds"
MYELOIDS_SEURAT_PATH    <- "/path/to/seurat_objects/myeloids.rds"
HEPATOCYTES_SEURAT_PATH <- "/path/to/seurat_objects/hepatocytes.rds"
PLASMAS_SEURAT_PATH     <- "/path/to/seurat_objects/plasmas.rds"
TCELLS_SEURAT_PATH      <- "/path/to/seurat_objects/tcells.rds"
PARENQUIMAL_SEURAT_PATH <- "/path/to/seurat_objects/parenquimal.rds"

# Output Path
OUT_DIR <- "figures/outs/"

# Ensure output directory exists
if (!dir.exists(OUT_DIR)) {
  dir.create(OUT_DIR, recursive = TRUE)
}

# ==============================================================================
# 1. Load Pre-Annotated Seurat Object & Define Palettes
# ==============================================================================
message("Loading annotated Seurat objects and defining palettes...")
seu <- readRDS(SEURAT_PATH)
meta <- seu@meta.data

# Set the factor order for etiology consistently
meta$etiology <- factor(meta$etiology, levels = c("HC", "HBV", "HDV RNA+", "HDV RNA-"))
seu@meta.data <- meta

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

# Broad annotation colours.
wide <- c(
  hepatocytes     = refined_col[["Hepatocytes"]],
  myeloids        = refined_col[["KC2"]],
  plasmas         = refined_col[["Plasma_cells"]],
  tcells          = refined_col[["Tem_Trm_cytotoxic_T_cells"]],
  non_parenquimal = refined_col[["Fibroblasts"]]
)

# Refined annotation colours mapping matching labels in annotation metadata
refined_col_annot <- c(
  "Memory B cells"          = refined_col[["Memory_B_cells"]],
  "pDC"                     = "#7A7A7A",
  "Plasma cells"            = refined_col[["Plasma_cells"]],
  "Naive B cells"           = refined_col[["Naive_B_cells"]],
  "Cycling B lineage cells" = refined_col[["Cycling_B_lineage_cells"]],
  "NKT cells"               = refined_col[["NKT_cells"]],
  "Tem cytotoxic T cells"   = refined_col[["Tem_Trm_cytotoxic_T_cells"]],
  "Rb high"                 = "#6D4C41",
  "Trm cytotoxic T cells"   = refined_col[["Tem_Trm_cytotoxic_T_cells"]],
  "Effector helper T cells" = refined_col[["Effector_helper_T_cells"]],
  "Naive T cells"           = refined_col[["Naive_T_cells"]],
  "Regulatory T cells"      = refined_col[["Regulatory_T_cells"]],
  "Cycling T cells"         = "#5E60CE",
  "Gamma-delta T cells"     = refined_col[["Gamma_delta_T_cells"]],
  "NK cells"                = refined_col[["NK_cells"]],
  "Cholangiocytes"          = refined_col[["Cholangiocytes"]],
  "Hepatocyte 1"            = refined_col[["Hepatocyte_1"]],
  "Hepatocyte 2"            = refined_col[["Hepatocyte_2"]],
  "Hepatocyte 5"            = refined_col[["Hepatocyte_5"]],
  "Hepatocyte 6"            = refined_col[["Hepatocyte_6"]],
  "Hepatocyte 3"            = refined_col[["Hepatocyte_3"]],
  "Hepatocyte 4"            = refined_col[["Hepatocyte_4"]],
  "KC1"                     = refined_col[["KC1"]],
  "M2-LYVE1"                = refined_col[["M2_LYVE1"]],
  "DCs CD1C"                = refined_col[["DCs_CD1C"]],
  "Mast cells"              = "#B56576",
  "Monocytes"               = refined_col[["Monocytes"]],
  "Neutrophils"             = "#8C564B",
  "M1"                      = refined_col[["M1"]],
  "KC2"                     = refined_col[["KC2"]],
  "Fibroblasts"             = refined_col[["Fibroblasts"]],
  "Endothelial cells 1"     = refined_col[["Endothelial_cells_1"]],
  "Endothelial cells 2"     = refined_col[["Endothelial_cells_2"]],
  "Endothelial cells 4"     = refined_col[["Endothelial_cells_4"]],
  "Smooth muscle cells"     = refined_col[["Smooth_muscle_cells"]],
  "Endothelial cells 3"     = refined_col[["Endothelial_cells_3"]],
  "Schwann cells"           = "#2A9D8F"
)

# ==============================================================================
# 2. Reusable Plotting Functions
# ==============================================================================

# Helper for UMAP Generation
plot_umap <- function(seurat_obj, group_var, palette, out_file) {
  p <- DimPlot(seurat_obj, reduction = "umap", group.by = group_var, cols = palette, label = FALSE, repel = TRUE) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "right",
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
  png(filename = out_file, width = 8, height = 8, units = "in", res = 1200)
  print(p)
  dev.off()
}

# Helper for Composition Barplots
plot_composition_barplot <- function(seurat_obj, group_var, palette, out_file) {
  meta_subset <- seurat_obj@meta.data[, c("pathology", group_var, "sample")]
  meta_subset$pathology <- factor(meta_subset$pathology, levels = c("HC", "HBV", "HDV"))
  
  meta_aggregated <- meta_subset %>%
    dplyr::group_by(pathology, !!sym(group_var)) %>%
    dplyr::summarise(total_count = n(), .groups = "drop")
  
  p <- ggplot(meta_aggregated, aes(x = pathology, y = total_count, fill = !!sym(group_var))) +
    geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.5) +
    theme_classic() +
    scale_fill_manual(values = palette) +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(x = "Pathology", y = "Proportion", fill = tools::toTitleCase(group_var)) +
    ggtitle("") +
    RotatedAxis()
  
  png(filename = out_file, width = 8, height = 8, units = "in", res = 800)
  print(p)
  dev.off()
}

# Helper for Stacked Bar Plots (Spatial)
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

# Helper for Volcano Plots (Spatial)
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
  
  deg_results$diffexpressed <- "NO"
  deg_results$diffexpressed[deg_results$avg_log2FC > log2(1.2) & deg_results$p_val < 0.05] <- "UP"
  deg_results$diffexpressed[deg_results$avg_log2FC < -log2(1.2) & deg_results$p_val < 0.05] <- "DOWN"
  deg_results$diffexpressed[deg_results$avg_log2FC > log2(1.2) & deg_results$p_val_adj < 0.05] <- "UPP"
  deg_results$diffexpressed[deg_results$avg_log2FC < -log2(1.2) & deg_results$p_val_adj < 0.05] <- "DOWNN"
  
  deg_results$delabel <- ifelse(deg_results$diffexpressed != "NO", deg_results$genes, NA)
  deg_results$p_val <- ifelse(deg_results$p_val < 1e-300, 1e-300, deg_results$p_val)
  
  diff_colors <- c("UPP" = "#803800", "DOWNN" = "#003F54", "UP" = "#B47846", "DOWN" = "steelblue")
  
  genes_to_label <- deg_results$delabel[deg_results$diffexpressed %in% c("UPP", "DOWNN")]
  deg_results$delabel <- ifelse(deg_results$genes %in% genes_to_label, deg_results$genes, NA)
  
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
    geom_label_repel(aes(label = delabel), size = 18 / .pt, segment.color = "black", segment.size = 1,
                     label.padding = unit(0.4, "lines"), label.size = 1, fontface = "bold", na.rm = TRUE)
  
  return(list(plot = p, data = deg_results))
}


# ==============================================================================
# Sup Fig 2A - UMAPS scRNAseq
# ==============================================================================
message("Loading single cell objects for UMAPs...")
todas       <- readRDS(ALL_SEURAT_PATH)
myeloids    <- readRDS(MYELOIDS_SEURAT_PATH)
hepatocytes <- readRDS(HEPATOCYTES_SEURAT_PATH)
plasmas     <- readRDS(PLASMAS_SEURAT_PATH)
tcells      <- readRDS(TCELLS_SEURAT_PATH)
parenquimal <- readRDS(PARENQUIMAL_SEURAT_PATH)

message("Generating UMAPs...")
plot_umap(todas, "subset", wide, file.path(OUT_DIR, "supfig2a_todas_umap.png"))
plot_umap(myeloids, "annotation", refined_col_annot, file.path(OUT_DIR, "supfig2a_myeloids_umap.png"))
plot_umap(hepatocytes, "annotation", refined_col_annot, file.path(OUT_DIR, "supfig2a_hepatocytes_umap.png"))
plot_umap(plasmas, "annotation", refined_col_annot, file.path(OUT_DIR, "supfig2a_plasmas_umap.png"))
plot_umap(tcells, "annotation", refined_col_annot, file.path(OUT_DIR, "supfig2a_tcells_umap.png"))
plot_umap(parenquimal, "annotation", refined_col_annot, file.path(OUT_DIR, "supfig2a_parenquimal_umap.png"))


# ==============================================================================
# Sup Fig 2B - Barplots scRNAseq
# ==============================================================================
message("Generating composition barplots...")
plot_composition_barplot(hepatocytes, "subset", wide, file.path(OUT_DIR, "supfig2b_hepatocytes_subset_barplot.png"))
plot_composition_barplot(hepatocytes, "annotation", refined_col_annot, file.path(OUT_DIR, "supfig2b_hepatocytes_annot_barplot.png"))
plot_composition_barplot(myeloids, "annotation", refined_col_annot, file.path(OUT_DIR, "supfig2b_myeloids_annot_barplot.png"))
plot_composition_barplot(todas, "subset", wide, file.path(OUT_DIR, "supfig2b_todas_subset_barplot.png"))


# ==============================================================================
# Sup Fig 2C - Stacked Bar Plots (Spatial)
# ==============================================================================
message("Generating Stacked Bar Plots (Spatial)...")
plot_stacked_bars(df = meta, x_var = "etiology", fill_var = "wide", palette = wide_col, 
                  out_file = file.path(OUT_DIR, "supfig2c_stackplot_wide.png"))
plot_stacked_bars(df = meta[meta$new_anot == "Hepatocytes", ], x_var = "etiology", fill_var = "refined", palette = refined_col, 
                  out_file = file.path(OUT_DIR, "supfig2c_stackplot_hepatos.png"))
plot_stacked_bars(df = meta[meta$wide == "Myeloid_cells", ], x_var = "etiology", fill_var = "refined", palette = refined_col, 
                  out_file = file.path(OUT_DIR, "supfig2c_stackplot_myeloids.png"))


# ==============================================================================
# Sup Fig 2D - Volcano plots
# ==============================================================================
message("Generating Volcano Plots...")

# Part 1: scRNAseq Volcano (Hepatocytes HDV vs HBV)
Idents(hepatocytes) <- "pathology"
deg_HDV_vs_HBV <- FindMarkers(hepatocytes, ident.1 = "HDV", ident.2 = "HBV", logfc.threshold = 0.01)

deg_results <- deg_HDV_vs_HBV
deg_results$gene <- rownames(deg_results)
deg_results$p_val <- ifelse(deg_results$p_val < 1e-300, 1e-300, deg_results$p_val)

deg_results$sign <- "0"
deg_results$sign[deg_results$avg_log2FC > log2(1.2) & deg_results$p_val < 0.05] <- "UP"
deg_results$sign[deg_results$avg_log2FC < -log2(1.2) & deg_results$p_val < 0.05] <- "DW"
deg_results$sign[deg_results$avg_log2FC > log2(1.2) & deg_results$p_val_adj < 0.05] <- "UPP"
deg_results$sign[deg_results$avg_log2FC < -log2(1.2) & deg_results$p_val_adj < 0.05] <- "DWW"
deg_results$sign <- factor(deg_results$sign, levels = c("UPP", "UP", "0", "DW", "DWW"))

genes_of_interest <- c("IL32", "SOD1", "GAPDH", "IFITM2", "IFITM3", "APOO", "HLA-F", "CD74", "SERPINA1", "NOTEN1", "PTEN", "NEAT1", "RORA")
label_data <- deg_results[deg_results$gene %in% genes_of_interest, ]
data_up <- label_data[label_data$avg_log2FC >= 0, ]
data_dw <- label_data[label_data$avg_log2FC < 0, ]

colors_volcano <- c("UPP" = "#803800", "DWW" = "#003F54", "UP" = "#B47846", "DW" = "steelblue")

p_vol_sc <- ggplot(deg_results, aes(x = avg_log2FC, y = -log10(p_val), col = sign)) +
  geom_point(size = 1) +
  scale_color_manual(values = colors_volcano) +
  theme_classic() +
  theme(text = element_text(family = "Helvetica", size = 18), legend.position = "none", plot.title = element_text(face = "bold"),
        axis.line = element_line(linewidth = 0.5), axis.ticks.length = unit(0.1, "cm")) +
  geom_vline(xintercept = c(-log2(1.2), log2(1.2)), col = "black", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = "dashed") +
  ggtitle("Hepatocytes: HDV vs HBV") + xlab("avg_log2FC") + ylab("-log10(p_val)") +
  geom_point(data = data_up, shape = 21, color = "black", fill = "#911704", size = 3, stroke = 0.4) +
  geom_point(data = data_dw, shape = 21, color = "black", fill = "#376D38", size = 3, stroke = 0.4) +
  geom_label_repel(data = data_dw, aes(label = gene), size = 9 / .pt, fontface = "bold", color = "black", fill = "white",
                   segment.color = "black", box.padding = 0.35, point.padding = 0.2, force = 2, max.overlaps = Inf, min.segment.length = 0) +
  geom_label_repel(data = data_up, aes(label = gene), size = 9 / .pt, fontface = "bold", color = "black", fill = "white",
                   segment.color = "black", box.padding = 0.35, point.padding = 0.2, force = 2, max.overlaps = Inf, min.segment.length = 0)

png(filename = file.path(OUT_DIR, "supfig2d_volcano_hepatocytes.png"), width = 10, height = 10, units = "in", res = 1200)
print(p_vol_sc)
dev.off()

# Part 2: Spatial Transcriptomics Volcano (Hepatocytes HDV vs HBV)
volcano_res <- volcano(anot = "new_anot", ct = "Hepatocytes", dif_col = "etiology", seu_obj = seu, id1 = "HDV RNA+", id2 = "HBV")
png(filename = file.path(OUT_DIR, "supfig2d_volcano_spatial.png"), width = 11, height = 8, units = "in", res = 1200)
print(volcano_res$plot)
dev.off()


# ==============================================================================
# Sup Fig 2E - Hepatocyte pathways scRNAseq (GO Enrichment)
# ==============================================================================
message("Generating GO enrichment dot plot...")

deg_HDV_vs_HC <- FindMarkers(hepatocytes, ident.1 = "HDV", ident.2 = "HC", logfc.threshold = 0.01, pos.only = TRUE)
deg_HBV_vs_HC <- FindMarkers(hepatocytes, ident.1 = "HBV", ident.2 = "HC", logfc.threshold = 0.01, pos.only = TRUE)

deg_HDV_sig <- deg_HDV_vs_HC[deg_HDV_vs_HC$p_val_adj < 0.05, ]
deg_HBV_sig <- deg_HBV_vs_HC[deg_HBV_vs_HC$p_val_adj < 0.05, ]

entrez_HDV <- clusterProfiler::bitr(rownames(deg_HDV_sig), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
entrez_HBV <- clusterProfiler::bitr(rownames(deg_HBV_sig), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

go_ids <- c("GO:0016032", "GO:0009615", "GO:0019058", "GO:0140888", "GO:0034340", "GO:0034341", "GO:0034342", "GO:0032635", 
            "GO:0070102", "GO:0019882", "GO:0042110", "GO:0001816", "GO:0042116", "GO:0070371", "GO:0044839", "GO:0034612", 
            "GO:0001837", "GO:0070482", "GO:0006979", "GO:0007160", "GO:0098609", "GO:0048771", "GO:0042246", "GO:0006111", "GO:0006805")

ego_HDV <- as.data.frame(enrichGO(gene = entrez_HDV$ENTREZID, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", readable = TRUE))
ego_HBV <- as.data.frame(enrichGO(gene = entrez_HBV$ENTREZID, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", readable = TRUE))

ego_HDV_df <- ego_HDV[ego_HDV$ID %in% go_ids, ]
ego_HBV_df <- ego_HBV[ego_HBV$ID %in% go_ids, ]

ego_HDV_df$Comparison <- "HDV vs HC"
ego_HBV_df$Comparison <- "HBV vs HC"

combined <- rbind(ego_HDV_df, ego_HBV_df)
combined$EnrichmentScore <- -log10(combined$p.adjust)

p_enrich <- ggplot(combined, aes(x = Comparison, y = Description, color = EnrichmentScore, size = Count)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() +
  labs(color = "-log10(adj p)", size = "Gene Count") +
  theme(axis.text.y = element_text(size = 10))

png(filename = file.path(OUT_DIR, "supfig2e_enrichment.png"), width = 10, height = 12, units = "in", res = 1200)
print(p_enrich)
dev.off()

message("Supplementary Figure 2 generation complete.")