# ==============================================================================
# Script: fig6.R
# Description: Generates panels for Figure 6 (Volcano, Pathways, SCOTIA, Co-occurrence)
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(clusterProfiler)
  library(msigdbr)
  library(org.Hs.eg.db)
})

options(bitmapType = "cairo")

# ==============================================================================
# Configuration & Paths
# Update these paths according to your local environment
# ==============================================================================
SEURAT_PATH  <- "/path/to/Objects/seurats_annotated.RDS"
SCOTIA_PATH  <- "/path/to/SCOTIA/Results/all_int.csv"
COOCCUR_PATH <- "/path/to/Celltype_enrichment/enrichment_files/all.csv"

OUT_DIR <- "figures/outs/"

if (!dir.exists(OUT_DIR)) {
  dir.create(OUT_DIR, recursive = TRUE)
}

# ==============================================================================
# 1. Load Data, Palettes & Target Pathways
# ==============================================================================
seu <- readRDS(SEURAT_PATH)
meta <- seu@meta.data

hepatocytes <- c("Hepatocyte 1", "Hepatocyte 2", "Hepatocyte 3", "Hepatocyte 5", "Hepatocyte 6")

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

target_pathways <- c(
  "GOBP_VIRAL_PROCESS", "GOBP_RESPONSE_TO_VIRUS", "GOBP_INTERFERON_MEDIATED_SIGNALING_PATHWAY",
  "GOBP_RESPONSE_TO_TYPE_I_INTERFERON", "GOBP_RESPONSE_TO_TYPE_II_INTERFERON",
  "GOBP_RESPONSE_TO_TYPE_III_INTERFERON", "GOBP_INTERLEUKIN_6_MEDIATED_SIGNALING_PATHWAY",
  "GOBP_INTERLEUKIN_6_PRODUCTION", "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION",
  "GOBP_T_CELL_ACTIVATION", "GOBP_T_CELL_PROLIFERATION", "GOBP_CYTOKINE_PRODUCTION",
  "GOBP_TRANSFORMING_GROWTH_FACTOR_BETA2_PRODUCTION", "GOBP_MACROPHAGE_ACTIVATION",
  "GOBP_ERK1_AND_ERK2_CASCADE", "GOBP_MAPK_CASCADE", 
  "GOBP_REGULATION_OF_CYTOPLASMIC_PATTERN_RECOGNITION_RECEPTOR_SIGNALING_PATHWAY",
  "GOBP_CELL_CYCLE_G2_M_PHASE_TRANSITION", "GOBP_PHOSPHORYLATION", 
  "GOBP_RESPONSE_TO_TUMOR_NECROSIS_FACTOR", "GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION",
  "GOBP_RESPONSE_TO_OXYGEN_LEVELS", "GOBP_RESPONSE_TO_OXIDATIVE_STRESS",
  "GOBP_CELL_MATRIX_ADHESION", "GOBP_CELL_CELL_ADHESION", "GOBP_TISSUE_REMODELING",
  "GOBP_TISSUE_REGENERATION", "GOBP_LIPID_METABOLIC_PROCESS", "GOBP_FATTY_ACID_BETA_OXIDATION",
  "GOBP_REGULATION_OF_GLUCONEOGENESIS", "GOBP_XENOBIOTIC_METABOLIC_PROCESS"
)

# ==============================================================================
# 2. Helper Functions
# ==============================================================================
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
    geom_point(size = 2) + scale_color_manual(values = diff_colors) + theme_classic(base_size = 22) +
    guides(color = guide_legend(override.aes = list(shape = 1))) +
    theme(
      legend.position = "none", plot.title = element_blank(), axis.title.x = element_blank(), 
      axis.title.y = element_blank(), axis.line = element_line(linewidth = 2), 
      axis.ticks.length = unit(0.5, "cm"), axis.ticks = element_line(linewidth = 2), 
      text = element_text(family = "Helvetica")
    ) +
    geom_label_repel(aes(label = delabel), size = 18 / .pt, segment.color = "black", segment.size = 2,
                     label.padding = unit(0.4, "lines"), label.size = 1, fontface = "bold", na.rm = TRUE)
  
  return(list(plot = p, data = deg_results))
}

calc_interaction_freq <- function(int_df, target_cells, group_col, ref_level, condition_name) {
  total_counts <- int_df %>% group_by(!!sym(group_col)) %>% summarise(all_int = n(), .groups = "drop")
  
  sub_int <- int_df %>% filter(refined_receptor %in% target_cells | refined_source %in% target_cells)
  sub_counts <- sub_int %>% group_by(!!sym(group_col)) %>% summarise(target_int = n(), .groups = "drop")
  
  res <- total_counts %>% 
    left_join(sub_counts, by = group_col) %>%
    mutate(target_int = replace_na(target_int, 0), perc = (target_int / all_int) * 100)
  
  ref_val <- res$perc[res[[group_col]] == ref_level]
  if(length(ref_val) == 0 || ref_val == 0) stop(paste("Reference level", ref_level, "not found or is zero."))
  
  res$norm <- res$perc / ref_val
  res$Condition <- condition_name
  return(res)
}

# ==============================================================================
# Figure 6A: Volcano Plot (HDV RNA- vs HDV RNA+)
# ==============================================================================
res_6a <- volcano(anot = "new_anot", ct = "Hepatocytes", seu_obj = seu, dif_col = "etiology", id1 = "HDV RNA-", id2 = "HDV RNA+")

png(filename = file.path(OUT_DIR, "figure6_A_volcano.png"), width = 10, height = 7, units = "in", res = 800)
print(res_6a$plot)
dev.off()

# ==============================================================================
# Figure 6B: Pathway Enrichment (HDV RNA- vs HDV RNA+)
# ==============================================================================
hdv_deg <- res_6a$data

msigdb_hallmark <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
gene_sets <- msigdb_hallmark %>% dplyr::select(gs_name, entrez_gene)

# Upregulated
up_genes <- hdv_deg$genes[hdv_deg$diffexpressed %in% c("p.adj<0.05 & FC>1.2", "p.val<0.05 & FC>1.2")]
gene_list_up <- bitr(up_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

up_enrich <- as.data.frame(enricher(gene = gene_list_up$ENTREZID, TERM2GENE = gene_sets))
up_enrich <- up_enrich %>% 
  mutate(GeneRatio = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])), s1 = "Upregulated")

# Downregulated
down_genes <- hdv_deg$genes[hdv_deg$diffexpressed %in% c("p.adj<0.05 & FC<0.83", "p.val<0.05 & FC<0.83")]
gene_list_down <- bitr(down_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

down_enrich <- as.data.frame(enricher(gene = gene_list_down$ENTREZID, TERM2GENE = gene_sets))
down_enrich <- down_enrich %>% 
  mutate(GeneRatio = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])), s1 = "Downregulated")

# Combine and Filter
final_enrichment <- bind_rows(up_enrich, down_enrich) %>%
  filter(ID %in% target_pathways) %>%
  arrange(GeneRatio)

p_6b <- ggplot(final_enrichment, aes(x = s1, y = Description, size = GeneRatio, color = s1)) +
  geom_point() + scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +
  scale_size_continuous(name = "Gene ratio", range = c(3, 10), guide = guide_legend(override.aes = list(color = "black", fill = "black"))) +
  theme_bw(base_size = 15) + labs(x = "", y = "Pathway description", color = "Pathway") +
  theme(
    axis.line = element_line(linewidth = 1.2, colour = "black"), axis.ticks = element_line(linewidth = 1.2, colour = "black"),
    axis.ticks.length = unit(6, "pt"), panel.background = element_rect(fill = NA, colour = NA),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2), plot.title = element_blank()
  )

png(filename = file.path(OUT_DIR, "figure6_B_pathway.png"), width = 11, height = 8, units = "in", res = 1200)
print(p_6b)
dev.off()

# ==============================================================================
# Figure 6C: SCOTIA Interactions (HDV RNA+ vs HDV RNA-)
# ==============================================================================
all_int <- read_csv(SCOTIA_PATH, show_col_types = FALSE)

logy_map <- meta$etiology
names(logy_map) <- meta$cell_names
all_int$etiology <- logy_map[all_int$id_source]

cut_int_6c <- all_int %>% filter(etiology %in% c("HDV RNA+", "HDV RNA-"), likelihood > 0.5)

df_hep_6c <- calc_interaction_freq(cut_int_6c, hepatocytes, "etiology", "HDV RNA-", "Hepatocytes All")
df_kc1_6c <- calc_interaction_freq(cut_int_6c, c(hepatocytes, "KC1"), "etiology", "HDV RNA-", "Hep + KC1")
df_kc2_6c <- calc_interaction_freq(cut_int_6c, c(hepatocytes, "KC2"), "etiology", "HDV RNA-", "Hep + KC2")

all_data_6c <- bind_rows(df_hep_6c, df_kc1_6c, df_kc2_6c)
all_data_6c$etiology <- factor(all_data_6c$etiology, levels = c("HDV RNA+", "HDV RNA-"))
all_data_6c$Condition <- factor(all_data_6c$Condition, levels = c("Hepatocytes All", "Hep + KC1", "Hep + KC2"))

p_6c <- ggplot(all_data_6c, aes(x = Condition, y = norm, fill = etiology)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black", width = 0.7) +
  scale_fill_manual(values = c("HDV RNA-" = "#D866FE", "HDV RNA+" = "#fed766")) + 
  theme_minimal() +
  theme(
    panel.grid = element_blank(), axis.line = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1), axis.ticks.length = unit(0.2, "cm"),
    axis.text.x = element_text(size = 12, face = "bold", color = "black"), axis.text.y = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 12, face = "bold", color = "black"), axis.title.x = element_blank(), 
    legend.position = "top"
  ) + labs(y = "Normalized Frequency (Ref: HDV RNA-)", fill = "Etiology")

png(filename = file.path(OUT_DIR, "figure6_C_scotia.png"), width = 8, height = 6, units = "in", res = 400)
print(p_6c)
dev.off()

# ==============================================================================
# Figure 6D: Co-occurrence Enrichment
# ==============================================================================
all_cooccur <- read_csv(COOCCUR_PATH, show_col_types = FALSE)

all_kc1_6d <- all_cooccur %>% 
  filter(from == "Hepatocyte", to %in% c("KC1", "KC2"), etiology %in% c("HDV RNA+", "HDV RNA-"))

all_kc1_6d$interval_numeric <- factor(all_kc1_6d$bin, levels = sort(unique(all_kc1_6d$bin)))
all_kc1_6d$etiology <- factor(all_kc1_6d$etiology, levels = c("HDV RNA+", "HDV RNA-"))

p_6d <- ggplot(all_kc1_6d, aes(x = bin, y = enrichment, group = to, color = to)) +
  geom_smooth(alpha = 0.05, linewidth = 1.5, method = "loess", formula = y ~ x) + 
  geom_hline(yintercept = 0, color = "black", linewidth = 1, linetype = "dashed") +
  facet_wrap(~ etiology) + theme_linedraw() + scale_color_manual(values = refined_col) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10), labels = NULL) +
  theme(
    panel.grid = element_blank(), panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
    axis.line = element_blank(), legend.position = "none",
    axis.text.x = element_blank(), axis.text.y = element_blank(),
    axis.ticks.x = element_line(linewidth = 1.5, color = "black"),
    axis.ticks.y = element_line(linewidth = 1.5, color = "black"),
    axis.ticks.length = unit(6, "pt"), text = element_text(size = 20, face = "bold")
  ) + labs(x = NULL, y = NULL)

png(filename = file.path(OUT_DIR, "figure6_D_coocurrence.png"), width = 9.2, height = 6, units = "in", res = 1200)
print(p_6d)
dev.off()