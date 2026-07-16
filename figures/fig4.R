# ==============================================================================
# Script: fig4.R
# Description: Generates panels for Figure 4 (SCOTIA Interactions & Co-occurrence)
# ==============================================================================

# Load necessary libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(readxl)
  library(tidyr)
  library(stringr)
  library(circlize)
  library(liana)
  library(reticulate)
})

options(bitmapType = "cairo")

# ==============================================================================
# Configuration & Paths
# Update these paths according to your local environment
# ==============================================================================
# Input Paths
SEURAT_PATH  <- "/path/to/Objects/seurats_annotated.RDS"
SCOTIA_PATH  <- "/path/to/SCOTIA/Results/all_int.csv"
COOCCUR_PATH <- "/path/to/Celltype_enrichment/enrichment_files/all.csv"
LIANA_DIR    <- "/path/to/liana"
TENSOR_FILE  <- "/path/to/tensor.pkl"

# Output Paths (Relative to the repository root)
OUT_DIR <- "figures/outs/"

# Ensure output directory exists
if (!dir.exists(OUT_DIR)) {
  dir.create(OUT_DIR, recursive = TRUE)
}

# ==============================================================================
# 1. Load Data & Define Palettes
# ==============================================================================
message("Loading annotated Seurat object metadata...")
seu <- readRDS(SEURAT_PATH)
meta <- seu@meta.data

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

etiology_colors <- c("HC" = "#fe4a49", "HBV" = "#2ab7ca", "HDV RNA+" = "#fed766")
health_order <- c("HC", "HBV", "HDV")

# ==============================================================================
# Figure 4A: Chord plot hepatocytes per health 
# ==============================================================================

# Read LIANA results 
liana_res_hdv <- read_excel(file.path(LIANA_DIR, "liana_res_hdv.xlsx"))
liana_res_hbv <- read_excel(file.path(LIANA_DIR, "liana_res_hbv.xlsx"))
liana_res_hc  <- read_excel(file.path(LIANA_DIR, "liana_res_hc.xlsx"))

liana_res_hdv$disease <- "HDV"
liana_res_hbv$disease <- "HBV"
liana_res_hc$disease  <- "HC"

liana_all <- bind_rows(liana_res_hc, liana_res_hbv, liana_res_hdv)

# Column names 
source_col   <- "source"
target_col   <- "target"
ligand_col   <- "ligand_complex"
receptor_col <- "receptor_complex"
pval_col     <- "cellphone_pvals"
score_col    <- "magnitude_rank"

required_cols <- c(source_col, target_col, ligand_col, receptor_col, pval_col, score_col, "disease")
missing_cols <- setdiff(required_cols, colnames(liana_all))

if (length(missing_cols) > 0) {
  stop("Missing columns in liana_all: ", paste(missing_cols, collapse = ", "))
}

# Collapse hepatocyte subtypes 
hepatocyte_labels <- c("Hepatocyte 1", "Hepatocyte 2", "Hepatocyte 3", 
                       "Hepatocyte 4", "Hepatocyte 5", "Hepatocyte 6", "Hepatocytes")

collapse_hepatocytes <- function(x) {
  ifelse(x %in% hepatocyte_labels, "Hepatocytes", x)
}

liana2 <- liana_all %>%
  mutate(
    source2  = collapse_hepatocytes(.data[[source_col]]),
    target2  = collapse_hepatocytes(.data[[target_col]]),
    ligand   = .data[[ligand_col]],
    receptor = .data[[receptor_col]],
    pval     = as.numeric(.data[[pval_col]]),
    score    = as.numeric(.data[[score_col]])
  )

# Keep only hepatocyte outgoing significant interactions 
hep_sig <- liana2 %>%
  filter(source2 == "Hepatocytes" & !is.na(pval) & pval < 0.05) %>%
  mutate(
    lr_pair = paste(ligand, receptor, sep = " - "),
    interaction_type = ifelse(target2 == "Hepatocytes", "Hepatocytes -> Hepatocytes", "Hepatocytes -> partner"),
    partner = target2
  )

top_n_lr <- 100
top_lr <- hep_sig %>%
  group_by(disease) %>%
  slice_min(order_by = score, n = top_n_lr, with_ties = FALSE) %>%
  ungroup()

cat("\nTop outgoing hepatocyte LR interactions per disease:\n")
print(table(top_lr$disease))

all_ann_colors <- refined_col 

get_plot_colors <- function(labels) {
  missing <- setdiff(labels, names(all_ann_colors))
  if (length(missing) > 0) {
    extra_cols <- rep("#BDBDBD", length(missing))
    names(extra_cols) <- missing
    c(all_ann_colors, extra_cols)
  } else {
    all_ann_colors
  }
}

# Summarize for chord plot
chord_df <- top_lr %>% count(disease, source2, target2, name = "weight")

# Chord plot function 
plot_hep_chord <- function(df_sub, disease_name, color_map) {
  if (nrow(df_sub) == 0) {
    plot.new()
    title(main = paste0(disease_name, " (no interactions)"))
    return(invisible(NULL))
  }
  
  sectors <- unique(c(df_sub$source2, df_sub$target2))
  sectors <- c("Hepatocytes", setdiff(sectors, "Hepatocytes"))
  plot_cols <- get_plot_colors(sectors)
  
  circos.clear()
  circos.par(start.degree = 90, gap.degree = 6, track.margin = c(0.01, 0.01), cell.padding = c(0, 0, 0, 0))
  
  chordDiagram(
    x = df_sub[, c("source2", "target2", "weight")],
    order = sectors, grid.col = plot_cols[sectors], transparency = 0.25,
    directional = 1, direction.type = c("arrows"), link.arr.type = "big.arrow",
    annotationTrack = "grid", preAllocateTracks = list(track.height = 0.12)
  )
  
  circos.trackPlotRegion(
    track.index = 1, panel.fun = function(x, y) {
      sector_name <- get.cell.meta.data("sector.index")
      xlim <- get.cell.meta.data("xlim")
      ylim <- get.cell.meta.data("ylim")
      circos.text(x = mean(xlim), y = ylim[1] + 0.1, labels = sector_name,
                  facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.8)
    }, bg.border = NA
  )
  
  title(main = disease_name, cex.main = 1.25)
}

# Export chord plot 
png(filename = file.path(OUT_DIR, "fig4a.png"), width = 18, height = 6, units = "in", res = 1200)
par(mfrow = c(1, 3), mar = c(1, 1, 4, 1), oma = c(0, 0, 3, 0))

for (d in c("HC", "HBV", "HDV")) {
  df_sub <- chord_df %>% filter(disease == d)
  plot_hep_chord(df_sub, d, all_ann_colors)
}

mtext("Hepatocyte as sender: top 100 interactions across etiologies", outer = TRUE, side = 3, line = 0, cex = 1.5, font = 2)
circos.clear()
dev.off()

# ==============================================================================
# Figure 4B: FACTOR 5 (cell2cell Tensor Analysis)
# ==============================================================================

message("Loading cell2cell tensor and plotting Factor 5 sender-receiver loadings...")

py$TENSOR_FILE <- TENSOR_FILE
py$OUT_FILE <- file.path(OUT_DIR, "fig4B_cell2cell_factor5_loadings.pdf")

py_run_string(r"(
import pickle
import matplotlib.pyplot as plt
import cell2cell as c2c

with open(TENSOR_FILE, 'rb') as f:
    tensor = pickle.load(f)

# Select factor to plot
selected_factor = 'Factor 5'

# Get loadings
loading_product = c2c.analysis.tensor_downstream.get_joint_loadings(
    tensor.factors,
    dim1='Sender Cells',
    dim2='Receiver Cells',
    factor=selected_factor
)

# Renaming cell type
loading_product = loading_product.rename(
    index={'Smooth muscle cells': 'Hepatic stellate cells'},
    columns={'Smooth muscle cells': 'Hepatic stellate cells'}
)

# Plotting Clustermap
lprod_cm = c2c.plotting.loading_clustermap(
    loadings=loading_product.T,   # Remove .T if you want the opposite orientation
    use_zscore=False,
    figsize=(8, 8),
    filename=OUT_FILE,
    cbar_label='Loading Product'
)
plt.close('all')
print(f'Saved Factor 5 sender-receiver loading clustermap to: {OUT_FILE}')
)")

# Right panel (Context Loadings)
py_run_string(r"(
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# -----------------------
# Inputs
# -----------------------
group_order = ['HC', 'HBV', 'HDV']
factor_human = 5
factor_idx = factor_human - 1  

palette = {
    "HC": "#fe4a49",
    "HBV": "#2ab7ca",
    "HDV": "#fed766"
}

# -----------------------
# Context loadings
# -----------------------
context_loadings = tensor.factors['Contexts']

try:
    import torch
    if isinstance(context_loadings, torch.Tensor):
        context_loadings = context_loadings.detach().cpu().numpy()
except Exception:
    pass

context_loadings = np.asarray(context_loadings)
contexts = list(metadict.keys())

assert len(contexts) == context_loadings.shape[0], (
    f"metadict has {len(contexts)} contexts but loadings has {context_loadings.shape[0]}. "
    "Provide the correct ordered context list."
)

df = pd.DataFrame({
    "Context": contexts,
    "Group": [metadict[c] for c in contexts],
    "Loading": context_loadings[:, factor_idx]
})

df = df[df["Group"].isin(group_order)].copy()
df["Group"] = pd.Categorical(df["Group"], categories=group_order, ordered=True)

# Plot
fig, ax = plt.subplots(figsize=(4, 3))

sns.boxplot(
    data=df, x="Group", y="Loading", order=group_order, hue="Group",
    palette=palette, dodge=False, showfliers=False, ax=ax, legend=False
)

sns.stripplot(
    data=df, x="Group", y="Loading", order=group_order,
    ax=ax, color="black", size=3, jitter=0.15, alpha=0.8
)

ax.set_title(f"Context factor {factor_human}")
ax.set_xlabel("")
ax.set_ylabel("Context loading")

ymin, ymax = ax.get_ylim()
y_range = ymax - ymin
ax.set_ylim(ymin, ymax + 0.3 * y_range)
ax.set_ylim(ymin, ymax + 0.1)

plt.tight_layout()

# Save
out = "context_fact5_1200dpi_colors.png"
fig.savefig(out, dpi=1200, bbox_inches="tight")
plt.close(fig)
print(f"Saved: {out}")
)")

# ==============================================================================
# Figure 4C: Hepatocyte-KC interaction proportion
# ==============================================================================
health_list <- list(HC = liana_res_hc, HBV = liana_res_hbv, HDV = liana_res_hdv)
kc_order <- c("KC2", "KC1")

is_hep <- function(x) { grepl("^Hepatocyte", x) }

get_hep_kc_counts <- function(df, condition_name) {
  df %>%
    filter((is_hep(source) & target %in% c("KC1", "KC2")) | (source %in% c("KC1", "KC2") & is_hep(target))) %>%
    mutate(kc = ifelse(source %in% c("KC1", "KC2"), source, target)) %>%
    count(kc, name = "n_interactions") %>%
    mutate(condition = condition_name)
}

count_df <- bind_rows(lapply(names(health_list), function(cond) { get_hep_kc_counts(health_list[[cond]], cond) })) %>%
  complete(condition = health_order, kc = c("KC1", "KC2"), fill = list(n_interactions = 0)) %>%
  mutate(condition = factor(condition, levels = health_order), kc = factor(kc, levels = kc_order)) %>%
  arrange(condition, kc)

total_interactions_df <- bind_rows(lapply(names(health_list), function(cond) {
  tibble(condition = cond, total_interactions = nrow(health_list[[cond]]))
})) %>% mutate(condition = factor(condition, levels = health_order))

# Divide Hepatocyte-KC interactions by total interactions
count_df <- count_df %>%
  left_join(total_interactions_df, by = "condition") %>%
  mutate(pct_total_interactions = n_interactions / total_interactions)

# Normalize each KC subtype to HC
hc_ref <- count_df %>% filter(condition == "HC") %>% dplyr::select(kc, pct_total_HC = pct_total_interactions)

plot_df <- count_df %>%
  left_join(hc_ref, by = "kc") %>%
  mutate(pct_vs_HC = ifelse(pct_total_HC > 0, (pct_total_interactions / pct_total_HC) * 100, NA_real_)) %>%
  arrange(condition, kc)

p <- ggplot(plot_df, aes(x = kc, y = pct_vs_HC, fill = condition)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65, color = "#333333", linewidth = 0.6) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "grey45", linewidth = 0.6) +
  scale_fill_manual(values = etiology_colors) +
  labs(title = "Hepatocyte ↔ KC subtype interactions normalized to HC", x = "KC subtype", y = "Interaction proportion (% of HC)", fill = NULL) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold")
  )

png(filename = file.path(OUT_DIR, "fig4c.png"), width = 10, height = 8, units = "in", res = 1200)
print(p)
dev.off()

# ==============================================================================
# Figure 4E: Hepatocyte - KC interactions
# ==============================================================================
pval_threshold <- 0.05
plot_width  <- 2.8
plot_height <- 4.2
dot_size_range <- c(0.18, 1.4)
y_text_size <- 3.7
strip_text_size <- 5.5
title_size <- 7
legend_text_size <- 4.3
legend_title_size <- 5.2

# Selected interactions
hep_to_kc1_pairs <- c("CCN1 → ITGB2", "CD99 → PILRA", "SAA1 → FPR1", "SAA1 → TLR2", "SAA1 → CD36", "HP → TLR4", 
                      "HP → ITGAM", "HP → CD163", "HP → ITGB2", "SERPINA1 → LRP1", "ALB → B2M-FCGRT", "VTN → CD47", 
                      "FN1 → CD44", "APOB → ITGB2", "KNG1 → ITGB2", "PLG → ITGB2")

hep_to_kc2_pairs <- c("ARF6 → SMAP1", "BMP1 → BMPR1A", "ADAM17 → IL6R", "SAA1 → TLR2", "SAA1 → FPR1", "SAA1 → CD36", 
                      "SAA1 → SCARB1", "HP → TLR4", "HP → ITGAM", "HP → CD163", "ALB → B2M-FCGRT", "APOA2 → LRP1", 
                      "APOA1 → ABCA1", "APOA1 → LRP1", "APOC3 → LRP1")

# Helpers
clean_lr_label <- function(x) { stringr::str_squish(stringr::str_replace_all(x, "_", "-")) }

add_liana_plot_columns <- function(df, disease_name) {
  df %>% mutate(
    disease = disease_name, original_source = source, original_target = target,
    ligand.complex = clean_lr_label(ligand_complex), receptor.complex = clean_lr_label(receptor_complex),
    lr_pair_clean = paste(ligand.complex, receptor.complex, sep = " → "),
    significant = !is.na(cellphone_pvals) & cellphone_pvals < pval_threshold
  )
}

collapse_lr_pairs <- function(df) {
  df %>% group_by(disease, lr_pair_clean) %>% arrange(desc(lr_means), cellphone_pvals, .by_group = TRUE) %>% slice_head(n = 1) %>% ungroup()
}

prepare_hep_to_kc_dotplot <- function(health_list, kc_name, selected_pairs) {
  selected_pairs_clean <- clean_lr_label(selected_pairs)
  panel_label <- paste0("Hepatocytes_to_", kc_name)
  
  liana_df <- bind_rows(lapply(names(health_list), function(disease_name) {
    health_list[[disease_name]] %>% filter(grepl("^Hepatocyte", source), target == kc_name) %>% add_liana_plot_columns(disease_name)
  }))
  
  plot_df <- liana_df %>%
    filter(lr_pair_clean %in% selected_pairs_clean) %>%
    collapse_lr_pairs() %>%
    mutate(source = factor(disease, levels = health_order), target = factor(panel_label), lr_pair_clean = factor(lr_pair_clean, levels = selected_pairs_clean)) %>%
    arrange(lr_pair_clean, disease) %>%
    as.data.frame()
  
  list(plot_res = plot_df, panel_label = panel_label)
}

plot_liana_interactions <- function(dotplot_obj, plot_title) {
  panel_target <- unique(as.character(dotplot_obj$plot_res$target))
  
  liana_dotplot(
    liana_res = dotplot_obj$plot_res, source_groups = health_order, target_groups = panel_target, ntop = NULL,
    magnitude = "lr_means", specificity = "cellphone_pvals", invert_specificity = TRUE, invert_magnitude = FALSE,
    colour.label = "LR mean\nexpression", size.label = "CPDB\np-value", show_complex = TRUE, size_range = dot_size_range
  ) + ggtitle(plot_title) +
    theme(
      plot.title = element_text(size = title_size, face = "bold", hjust = 0.5), strip.text.x = element_text(size = strip_text_size, face = "plain"),
      axis.text.y = element_text(size = y_text_size), axis.title.y = element_text(size = 5.2, face = "bold"),
      axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(),
      legend.title = element_text(size = legend_title_size, face = "bold"), legend.text = element_text(size = legend_text_size),
      legend.key.size = unit(0.20, "cm"), plot.margin = margin(2, 2, 2, 2)
    )
}

# Generate plots
hep_to_kc1_dot <- prepare_hep_to_kc_dotplot(health_list, "KC1", hep_to_kc1_pairs)
hep_to_kc2_dot <- prepare_hep_to_kc_dotplot(health_list, "KC2", hep_to_kc2_pairs)

p_hep_to_kc1 <- plot_liana_interactions(hep_to_kc1_dot, "Hepatocytes → KC1")
p_hep_to_kc2 <- plot_liana_interactions(hep_to_kc2_dot, "Hepatocytes → KC2")

# Export plots
ggsave(filename = file.path(OUT_DIR, "fig4e_hepatocytes_to_kc1.png"), plot = p_hep_to_kc1, width = plot_width, height = plot_height, dpi = 1200)
ggsave(filename = file.path(OUT_DIR, "fig4e_hepatocytes_to_kc1.pdf"), plot = p_hep_to_kc1, width = plot_width, height = plot_height)
ggsave(filename = file.path(OUT_DIR, "fig4e_hepatocytes_to_kc2.png"), plot = p_hep_to_kc2, width = plot_width, height = plot_height, dpi = 1200)
ggsave(filename = file.path(OUT_DIR, "fig4e_hepatocytes_to_kc2.pdf"), plot = p_hep_to_kc2, width = plot_width, height = plot_height)