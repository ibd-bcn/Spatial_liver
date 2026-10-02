# ==============================================================================
# Script: fig1.R
# Description: Generates panels for Figure 1 (Spatial Polygons, Volcano, Pathways)
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
  library(DOSE)
  library(grid)
  library(limma)
  library(stringr)
  library(readxl)
  library(ggvenn)    # Added for Fig 1D
  library(pheatmap)  # Added for Fig 1E
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
  "Monocytes"                 = "#a11191", "i-KC"                       = "#ff0077",
  "h-Mac"                  = "#ffbf00", "h-KC"                       = "#56c9f2",
  "DCs_CD1C"                  = "#02fa8d", "i-Mac"                        = "#9c78fe",
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
# FIG 1A - Spatial Polygons
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
plot_spatial_fov(file.path(POLY_DIR, "Slide_1.csv"), "Slide_1", 1, meta, "refined_hep", refined_col, file.path(OUT_DIR, "cslide1_fov1_refined.png"))
plot_spatial_fov(file.path(POLY_DIR, "Slide_2.csv"), "Slide_2", 16, meta, "refined_hep", refined_col, file.path(OUT_DIR, "cslide2_fov16_refined.png"))

# 3. Refined Annotations (Myeloids Only)
myeloid_types <- c("h-Mac", "i-KC", "h-KC", "Monocytes", "DCs_CD1C", "i-Mac")
meta$refined_myel <- ifelse(meta$refined %in% myeloid_types, meta$refined, "Other")

plot_spatial_fov(file.path(POLY_DIR, "Slide_2.csv"), "Slide_2", 2, meta, "refined_myel", refined_col, file.path(OUT_DIR, "cslide2_fov2_refined_myel.png"))
plot_spatial_fov(file.path(POLY_DIR, "Slide_1.csv"), "Slide_1", 1, meta, "refined_myel", refined_col, file.path(OUT_DIR, "cslide1_fov1_refined_myel.png"))
plot_spatial_fov(file.path(POLY_DIR, "Slide_2.csv"), "Slide_2", 16, meta, "refined_myel", refined_col, file.path(OUT_DIR, "cslide2_fov16_refined_myel.png"))

# ==============================================================================
# FIG 1B - Volcano & Enrichment Analysis
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

png(filename = file.path(OUT_DIR, "fig1_pathway_enrichment.png"), width = 11, height = 8, units = "in", res = 1200)
print(p_path)
dev.off()

write.csv(final_enrichment, file.path(OUT_DIR, "fig1_pathway_cosmx.csv"), row.names = FALSE)

# ==============================================================================
# FIG 1C - scRNAseq HDV-active vs HBV 
# ==============================================================================

hepatocytes <- readRDS("/path/to/hepatocytes_harmony.RDS")

# Differential expression analysis: HDV vs HBV hepatocytes
hep <- hepatocytes
Idents(hep) <- "pathology"

deg_HDV_vs_HBV <- FindMarkers(
  hep,
  ident.1 = "HDV",
  ident.2 = "HBV",
  logfc.threshold = 0.01
)

deg_HDV_sig <- deg_HDV_vs_HBV[
  deg_HDV_vs_HBV$p_val_adj < 0.05 &
    deg_HDV_vs_HBV$avg_log2FC > 1.2,
]

deg_HBV_sig <- deg_HDV_vs_HBV[
  deg_HDV_vs_HBV$p_val_adj < 0.05 &
    deg_HDV_vs_HBV$avg_log2FC < -1.2,
]

# Convert gene symbols to Entrez IDs 
genes_HDV <- rownames(deg_HDV_sig)
genes_HBV <- rownames(deg_HBV_sig)

entrez_HDV <- clusterProfiler::bitr(
  genes_HDV,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

entrez_HBV <- clusterProfiler::bitr(
  genes_HBV,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# Selected GO Biological Process terms 
go_ids <- c(
  "GO:0016032",  # viral process
  "GO:0019058",  # viral life cycle / viral infectious cycle
  "GO:0019079",  # viral genome replication
  "GO:0034612",  # response to TNF
  "GO:0034341",  # response to type II interferon
  "GO:0071346",  # cellular response to type II interferon
  "GO:0006979",  # response to oxidative stress
  "GO:0006935",  # chemotaxis
  "GO:0060326",  # cell chemotaxis
  "GO:0006638",  # triglyceride metabolic process
  "GO:0019216",  # regulation of lipid metabolic process
  "GO:0140467",  # integrated stress response signaling
  "GO:0006457",  # protein folding
  "GO:0061077",  # chaperone-mediated protein folding
  "GO:0022604"   # regulation of cell morphogenesis
)

# GO Biological Process enrichment 
ego_HDV <- enrichGO(
  gene = entrez_HDV$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  readable = TRUE
)

ego_HBV <- enrichGO(
  gene = entrez_HBV$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  readable = TRUE
)

# Format enrichment results 
ego_HDV_df <- as.data.frame(ego_HDV)
ego_HBV_df <- as.data.frame(ego_HBV)

ego_HDV_df <- ego_HDV_df[ego_HDV_df$ID %in% go_ids, ]
ego_HBV_df <- ego_HBV_df[ego_HBV_df$ID %in% go_ids, ]

ego_HDV_df$Comparison <- "Upregulated"
ego_HBV_df$Comparison <- "Downregulated"

combined <- rbind(ego_HDV_df, ego_HBV_df)
combined$EnrichmentScore <- -log10(combined$p.adjust)

# Plot pathway enrichment 
final_plot <- ggplot(
  combined,
  aes(x = Comparison, y = Description, color = Comparison, size = Count)
) +
  geom_point() +
  scale_size(range = c(4, 10)) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +
  theme_bw() +
  labs(
    title = "", x = "", y = "Pathway description",
    size = "Gene count", color = "Comparison"
  ) +
  theme(
    axis.text = element_text(size = 15),
    legend.text = element_text(size = 15),
    legend.title = element_text(face = "bold", size = 15),
    axis.line = element_line(linewidth = 1.2, colour = "black"),
    axis.ticks = element_line(linewidth = 1.2, colour = "black"),
    axis.ticks.length = unit(6, "pt"),
    panel.background = element_rect(fill = NA, colour = NA),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2)
  )

# Export figure
png(filename = file.path(OUT_DIR, "fig1c.png"), width = 11, height = 8, units = "in", res = 1200)
print(final_plot)
dev.off()


# ==============================================================================
# FIG 1D - nCounter - Venn diagram 
# ==============================================================================

ALL_normalized_data_with_clinical_variables_Nanostring <- read_excel("/path/to/file")
df <- ALL_normalized_data_with_clinical_variables_Nanostring

# Format expression matrix to the only data we need 
df <- df[1:39, ]
df <- as.data.frame(df)
colnames(df) <- df[1, ]
df <- df[2:39, ]
rownames(df) <- df[, 1]
df <- df[, 2:783]
df <- t(df)
df <- df[10:782, ]

expr_matrix <- df
mode(expr_matrix) <- "numeric"
sample_names <- colnames(expr_matrix)

# Create sample metadata
condition <- ifelse(
  str_starts(sample_names, "D"), "delta",
  ifelse(
    str_starts(sample_names, "C"), "control",
    ifelse(
      str_starts(sample_names, "N"), "negative",
      ifelse(str_starts(sample_names, "B"), "vb", NA)
    )
  )
)

metadata <- data.frame(
  sample = sample_names,
  condition = factor(condition)
)
rownames(metadata) <- metadata$sample
expr_matrix <- expr_matrix[, metadata$sample]

# Differential expression: HDV RNA+ vs Control
metadata2 <- metadata[metadata$condition %in% c("delta", "control"), ]
metadata2$condition <- droplevels(metadata2$condition)
expr_matrix2 <- expr_matrix[, metadata2$sample]

design <- model.matrix(~0 + condition, data = metadata2)
colnames(design) <- levels(metadata2$condition)
fit <- lmFit(expr_matrix2, design)
contrast_matrix <- makeContrasts(delta - control, levels = design)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)
deg_table <- topTable(fit2, number = Inf, adjust.method = "BY", sort.by = "P")
deg_sig <- deg_table[deg_table$adj.P.Val < 0.05 & abs(deg_table$logFC) > 1, ]

# Differential expression: HDV RNA− vs Control
metadata2 <- metadata[metadata$condition %in% c("negative", "control"), ]
metadata2$condition <- droplevels(metadata2$condition)
expr_matrix2 <- expr_matrix[, metadata2$sample]

design <- model.matrix(~0 + condition, data = metadata2)
colnames(design) <- levels(metadata2$condition)
fit <- lmFit(expr_matrix2, design)
contrast_matrix <- makeContrasts(negative - control, levels = design)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)
deg_table <- topTable(fit2, number = Inf, adjust.method = "BY", sort.by = "P")
deg_signeg <- deg_table[deg_table$adj.P.Val < 0.05 & abs(deg_table$logFC) > 1, ]

# Differential expression: HBV vs Control
metadata2 <- metadata[metadata$condition %in% c("vb", "control"), ]
metadata2$condition <- droplevels(metadata2$condition)
expr_matrix2 <- expr_matrix[, metadata2$sample]

design <- model.matrix(~0 + condition, data = metadata2)
colnames(design) <- levels(metadata2$condition)
fit <- lmFit(expr_matrix2, design)
contrast_matrix <- makeContrasts(vb - control, levels = design)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)
deg_table <- topTable(fit2, number = Inf, adjust.method = "BY", sort.by = "P")
deg_sigvb <- deg_table[deg_table$adj.P.Val < 0.05 & abs(deg_table$logFC) > 1, ]

# Venn diagram
up_delta <- rownames(deg_sig[deg_sig$adj.P.Val < 0.05 & deg_sig$logFC > 1.2, ])
up_neg <- rownames(deg_signeg[deg_signeg$adj.P.Val < 0.05 & deg_signeg$logFC > 1.2, ])
up_vb <- rownames(deg_sigvb[deg_sigvb$adj.P.Val < 0.05 & deg_sigvb$logFC > 1.2, ])

venn_list <- list(
  "HDV RNA +" = up_delta,
  "HDV RNA -" = up_neg,
  "HBV"       = up_vb
)

p <- ggvenn(
  venn_list,
  fill_color = c("#fed766", "#D866FE", "#2ab7ca"),
  stroke_size = 0.5,
  set_name_size = 5,
  text_size = 4,
  show_percentage = FALSE
)

# Export figure
png(filename = file.path(OUT_DIR, "fig1d.png"), width = 4, height = 4, units = "in", res = 1200)
p
dev.off()


# ==============================================================================
# FIG 1E - nCounter Heatmap 
# ==============================================================================

# Pathway file
pathway_file <- file.path("/path/to/file")

# Read selected pathways 
lines <- readLines(pathway_file, warn = FALSE)
lines <- lines[nzchar(trimws(lines))]

parts <- strsplit(lines, "\\s+")
pathways_raw <- do.call(rbind, lapply(parts, function(x) c(x[1], x[2], x[3])))
pathways_raw <- as.data.frame(pathways_raw, stringsAsFactors = FALSE)
colnames(pathways_raw) <- c("Pathway", "GO", "Module")

pathways_raw$Pathway_clean <- pathways_raw$Pathway %>%
  gsub("^GOBP_", "", .) %>%
  gsub("_", " ", .)

go_ids <- unique(pathways_raw$GO)
go_ids <- go_ids[grepl("^GO:", go_ids)]

# Get GO:BP genes from MSigDB 
m_bp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
m_bp_sub <- m_bp %>% filter(gs_exact_source %in% go_ids)

go2genes <- split(m_bp_sub$gene_symbol, m_bp_sub$gs_exact_source)
pathway2genes <- setNames(
  lapply(pathways_raw$GO, function(go) unique(go2genes[[go]])),
  pathways_raw$Pathway_clean
)

# Clean pathways and match to NanoString panel
pathway2genes$`REGULATION OF CYTOPLASMIC PATTERN RECOGNITION RECEPTOR SIGNALING PATHWAY` <- NULL

rownames(expr_matrix) <- gsub("-mRNA", "", rownames(expr_matrix))
genes_measured <- rownames(expr_matrix)

pathway2genes <- lapply(pathway2genes, intersect, y = genes_measured)
pathway_sizes <- sapply(pathway2genes, length)
pathway2genes <- pathway2genes[pathway_sizes >= 3]

# Calculate pathway scores 
expr_z <- t(scale(t(expr_matrix)))
expr_z[is.na(expr_z)] <- 0

pathway_scores <- vapply(
  names(pathway2genes),
  function(pw) {
    colMeans(expr_z[pathway2genes[[pw]], , drop = FALSE])
  },
  FUN.VALUE = numeric(ncol(expr_z))
)

pathway_scores <- t(pathway_scores)
pathway_z <- t(scale(t(pathway_scores)))
pathway_z[is.na(pathway_z)] <- 0

# Heatmap annotation 
metadata$condition2 <- recode(
  metadata$condition,
  control  = "HC",
  vb       = "HBV",
  delta    = "HDV RNA+",
  negative = "HDV RNA-"
)

metadata$condition2 <- factor(
  metadata$condition2,
  levels = c("HC", "HBV", "HDV RNA+", "HDV RNA-")
)

sample_order <- rownames(metadata)[order(metadata$condition2)]

annotation_col <- data.frame(Condition = metadata[sample_order, "condition2"])
rownames(annotation_col) <- sample_order

ann_colors <- list(
  Condition = c(
    "HC"       = "#fe4a49",
    "HBV"      = "#2ab7ca",
    "HDV RNA+" = "#fed766",
    "HDV RNA-" = "#D866FE"
  )
)

heat_cols <- colorRampPalette(c("#2166ac", "white", "#b2182b"))(100)

# Export heatmap 
png(filename = file.path(OUT_DIR, "fig1e.png"), width = 8.5, height = 3.5, units = "in", res = 800)

pheatmap(
  pathway_z[, sample_order, drop = FALSE],
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = heat_cols,
  border_color = NA,
  fontsize_row = 9,
  fontsize_col = 7
)

dev.off()