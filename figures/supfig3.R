# ==============================================================================
# Script: supfig3.R
# Description: Generates panels for Supplementary Figure 3 (LIANA, cell2cell Tensor)
# ==============================================================================

# Load necessary libraries
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(circlize)
  library(reticulate)
  library(liana)      # Added for Panel 3E
})

options(bitmapType = "cairo")

# ==============================================================================
# Configuration & Paths
# Update these paths according to your local environment
# ==============================================================================
# Input Paths
LIANA_DIR     <- "/path/to/liana"
TENSOR_FILE   <- "/path/to/cell2cell/tensor.pkl"
METADICT_FILE <- "/path/to/metadict.pkl"

# Output Path
OUT_DIR <- "figures/outs/"

# Ensure output directory exists
if (!dir.exists(OUT_DIR)) {
  dir.create(OUT_DIR, recursive = TRUE)
}

# ==============================================================================
# 1. Global Data Loading & Palettes
# ==============================================================================
message("Loading LIANA results globally for panels 3A and 3E...")

liana_res_hdv <- read_excel(file.path(LIANA_DIR, "liana_res_hdv.xlsx"))
liana_res_hbv <- read_excel(file.path(LIANA_DIR, "liana_res_hbv.xlsx"))
liana_res_hc  <- read_excel(file.path(LIANA_DIR, "liana_res_hc.xlsx"))

liana_res_hdv$disease <- "HDV"
liana_res_hbv$disease <- "HBV"
liana_res_hc$disease  <- "HC"

# Combined for chord plots (3A)
liana_all <- bind_rows(liana_res_hc, liana_res_hbv, liana_res_hdv)

# List for dotplots (3E)
health_list <- list(
  HC  = liana_res_hc,
  HBV = liana_res_hbv,
  HDV = liana_res_hdv
)

health_order <- c("HC", "HBV", "HDV")

# Define Global Colors
all_ann_colors <- c(
  "Hepatocytes" = "#E69F00", "Cholangiocytes" = "#4DBBD5",
  "Monocytes"   = "#1F78B4", "Neutrophils" = "#A6CEE3",
  "i-KC"         = "#33A02C", "h-KC"         = "#B2DF8A",
  "i-Mac"          = "#FB9A99", "M2-LYVE1"    = "#E31A1C",
  "DCs CD1C"    = "#FF7F00", "Mast cells"  = "#6A3D9A",
  "Memory B cells"          = "#FFD1A1", "pDC"                     = "#FFB066",
  "Plasma cells"            = "#FF8C42", "Naive B cells"           = "#E76F00",
  "Cycling B lineage cells" = "#B55200", "Memory_B_cells"          = "#E34183",
  "Plasma_cells"            = "#F1B8EF", "Naive_B_cells"           = "#FE64F9",
  "NKT cells"               = "#5F9DF7", "Tem cytotoxic T cells"   = "#00B4D8",
  "Rb high"                 = "#9B5DE5", "Trm cytotoxic T cells"   = "#00F5D4",
  "Effector helper T cells" = "#118AB2", "Naive T cells"           = "#ADE8F4",
  "Regulatory T cells"      = "#4361EE", "Cycling T cells"         = "#F72585",
  "Gamma-delta T cells"     = "#3A0CA3", "NK cells"                = "#4CC9F0",
  "Hepatic stellate cells"  = "#8E7DBE", "Fibroblasts"             = "#8E7DBE",
  "Endothelial cells 1"     = "#4CAF50", "Endothelial cells 2"     = "#7CB342",
  "Endothelial cells 3"     = "#009688", "Endothelial cells 4"     = "#C0CA33",
  "Smooth muscle cells"     = "#D81B60", "Schwann cells"           = "#795548"
)

# ==============================================================================
# FIG Sup 3A - Liana chord plot (Incoming to Hepatocytes)
# ==============================================================================
message("Generating Sup Fig 3A (Chord Plots)...")

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

hepatocyte_labels <- c("Hepatocyte 1", "Hepatocyte 2", "Hepatocyte 3", 
                       "Hepatocyte 4", "Hepatocyte 5", "Hepatocyte 6", "Hepatocytes")

collapse_hepatocytes <- function(x) { ifelse(x %in% hepatocyte_labels, "Hepatocytes", x) }

liana2 <- liana_all %>%
  mutate(
    source2  = collapse_hepatocytes(.data[[source_col]]),
    target2  = collapse_hepatocytes(.data[[target_col]]),
    ligand   = .data[[ligand_col]],
    receptor = .data[[receptor_col]],
    pval     = as.numeric(.data[[pval_col]]),
    score    = as.numeric(.data[[score_col]])
  )

# Filter for significant incoming interactions to hepatocytes
hep_sig <- liana2 %>%
  filter(target2 == "Hepatocytes", !is.na(pval), pval < 0.05) %>%
  mutate(
    lr_pair = paste(ligand, receptor, sep = " - "),
    interaction_type = ifelse(source2 == "Hepatocytes", "Hepatocytes -> Hepatocytes", "partner -> Hepatocytes"),
    partner = source2
  )

# Keep top 100 per disease
top_n_lr <- 100
top_lr <- hep_sig %>%
  group_by(disease) %>%
  slice_min(order_by = score, n = top_n_lr, with_ties = FALSE) %>%
  ungroup()

chord_df <- top_lr %>% count(disease, source2, target2, name = "weight")

# Chord Helpers
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
    x = df_sub[, c("source2", "target2", "weight")], order = sectors, grid.col = plot_cols[sectors], 
    transparency = 0.25, directional = 1, direction.type = c("arrows"), link.arr.type = "big.arrow", 
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

png(filename = file.path(OUT_DIR, "supfig3a_hepatocyte_incoming_chord_100.png"), width = 18, height = 6, units = "in", res = 1200)
par(mfrow = c(1, 3), mar = c(1, 1, 3, 1))

for (d in c("HC", "HBV", "HDV")) {
  plot_hep_chord(chord_df %>% filter(disease == d), d, all_ann_colors)
}
circos.clear()
dev.off()


# ==============================================================================
# FIG Sup 3B - Heatmap all factors (cell2cell)
# ==============================================================================
message("Generating Sup Fig 3B (cell2cell Ligand-Receptor Pair Loadings)...")

py$TENSOR_FILE <- TENSOR_FILE
py$OUT_FILE <- file.path(OUT_DIR, "supfig3b_ligand_receptor_pair_loadings_heatmap.pdf")

py_run_string(r"(
import pickle
import cell2cell as c2c
import matplotlib.pyplot as plt

with open(TENSOR_FILE, 'rb') as f:
    tensor = pickle.load(f)

c2c.plotting.loading_clustermap(
    loadings=tensor.factors['Ligand-Receptor Pairs'],
    loading_threshold=0.1,
    use_zscore=False,
    figsize=(28, 8),
    filename=OUT_FILE,
    row_cluster=False,
    tick_fontsize=12,
    dendrogram_ratio=0.15
)
plt.close('all')
print(f'Saved ligand-receptor loading heatmap to: {OUT_FILE}')
)")

# ==============================================================================
# FIG Sup 3C - Heatmap per factor (cell2cell)
# ==============================================================================
message("Generating Sup Fig 3C (cell2cell Sender-Receiver Joint Loadings)...")

py$TENSOR_FILE <- TENSOR_FILE
py$OUT_DIR <- OUT_DIR

py_run_string(r"(
import os
import pickle
import matplotlib.pyplot as plt
import cell2cell as c2c

with open(TENSOR_FILE, 'rb') as f:
    tensor = pickle.load(f)

selected_factors = ['Factor 1', 'Factor 2', 'Factor 3', 'Factor 4']

for selected_factor in selected_factors:
    loading_product = c2c.analysis.tensor_downstream.get_joint_loadings(
        tensor.factors,
        dim1='Sender Cells',
        dim2='Receiver Cells',
        factor=selected_factor
    )

    out_file = os.path.join(OUT_DIR, 'supfig3c_' + selected_factor.replace(' ', '_').lower() + '_sender_receiver_loadings.pdf')

    c2c.plotting.loading_clustermap(
        loadings=loading_product.T,
        use_zscore=False,
        figsize=(8, 8),
        filename=out_file,
        cbar_label='Loading Product'
    )
    plt.close('all')
    print(f'Saved: {out_file}')
)")

# ==============================================================================
# FIG Sup 3D - Boxplot Context Loadings (cell2cell)
# ==============================================================================
message("Generating Sup Fig 3D (cell2cell Context Loadings)...")

py$TENSOR_FILE <- TENSOR_FILE
py$METADICT_FILE <- METADICT_FILE
py$OUT_DIR <- OUT_DIR

py_run_string(r"(
import os
import pickle
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

with open(TENSOR_FILE, 'rb') as f:
    tensor = pickle.load(f)

with open(METADICT_FILE, 'rb') as f:
    metadict = pickle.load(f)

group_order = ['HC', 'HBV', 'HDV']
factors_to_plot = [1, 2, 3, 4]
palette = {'HC': '#fe4a49', 'HBV': '#2ab7ca', 'HDV': '#fed766'}

context_loadings = tensor.factors['Contexts']

try:
    import torch
    if isinstance(context_loadings, torch.Tensor):
        context_loadings = context_loadings.detach().cpu().numpy()
except Exception:
    pass

context_loadings = np.asarray(context_loadings)
contexts = list(metadict.keys())

for factor_human in factors_to_plot:
    factor_idx = factor_human - 1  
    
    df = pd.DataFrame({
        'Context': contexts,
        'Group': [metadict[c] for c in contexts],
        'Loading': context_loadings[:, factor_idx]
    })
    
    df = df[df['Group'].isin(group_order)].copy()
    df['Group'] = pd.Categorical(df['Group'], categories=group_order, ordered=True)
    
    fig, ax = plt.subplots(figsize=(4, 3))
    
    sns.boxplot(
        data=df, x='Group', y='Loading', order=group_order, hue='Group',
        palette=palette, dodge=False, showfliers=False, ax=ax, legend=False
    )
    
    sns.stripplot(
        data=df, x='Group', y='Loading', order=group_order,
        ax=ax, color='black', size=3, jitter=0.15, alpha=0.8
    )
    
    ax.set_title(f'Context factor {factor_human}')
    ax.set_xlabel('')
    ax.set_ylabel('Context loading')
    plt.tight_layout()
    
    out = os.path.join(OUT_DIR, f'supfig3d_context_factor{factor_human}_colors.png')
    fig.savefig(out, dpi=1200, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved: {out}')
)")

# ==============================================================================
# FIG Sup 3E - Dotplot KC -> Hepatocytes (LIANA)
# ==============================================================================
message("Generating Sup Fig 3E (KC -> Hepatocyte Dotplots)...")

pval_threshold <- 0.05
top_n_per_condition <- 4
plot_width  <- 2.8
plot_height <- 4.2
dot_size_range <- c(0.18, 1.4)
y_text_size        <- 3.7
strip_text_size    <- 5.5
title_size         <- 7
legend_text_size   <- 4.3
legend_title_size  <- 5.2

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

prepare_kc_to_hep_dotplot <- function(health_list, kc_name, top_n = 4) {
  panel_label <- paste0(kc_name, "_to_Hepatocytes")
  
  liana_df <- bind_rows(lapply(names(health_list), function(disease_name) {
    health_list[[disease_name]] %>% filter(source == kc_name, grepl("^Hepatocyte", target)) %>% add_liana_plot_columns(disease_name)
  }))
  
  liana_df <- collapse_lr_pairs(liana_df)
  
  top_pairs <- liana_df %>% filter(significant) %>% group_by(disease) %>%
    arrange(desc(lr_means), cellphone_pvals, .by_group = TRUE) %>%
    mutate(rank_within_condition = row_number()) %>% slice_head(n = top_n) %>% ungroup() %>%
    arrange(factor(disease, levels = health_order), rank_within_condition) %>% pull(lr_pair_clean) %>% unique()
  
  if (length(top_pairs) == 0) stop(paste0("No significant ", kc_name, " -> Hepatocyte interactions found."))
  
  plot_df <- liana_df %>% filter(lr_pair_clean %in% top_pairs) %>%
    mutate(source = factor(disease, levels = health_order), target = factor(panel_label), lr_pair_clean = factor(lr_pair_clean, levels = top_pairs)) %>%
    arrange(lr_pair_clean, disease) %>% as.data.frame()
  
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

# Generate and Export
kc1_to_hep_dot <- prepare_kc_to_hep_dotplot(health_list, "i-KC", top_n_per_condition)
kc2_to_hep_dot <- prepare_kc_to_hep_dotplot(health_list, "h-KC", top_n_per_condition)

p_kc1_to_hep <- plot_liana_interactions(kc1_to_hep_dot, "i-KC → Hepatocytes")
p_kc2_to_hep <- plot_liana_interactions(kc2_to_hep_dot, "h-KC → Hepatocytes")

ggsave(filename = file.path(OUT_DIR, "supfig3e_kc1_to_hepatocytes.png"), plot = p_kc1_to_hep, width = plot_width, height = plot_height, dpi = 1200)
ggsave(filename = file.path(OUT_DIR, "supfig3e_kc2_to_hepatocytes.png"), plot = p_kc2_to_hep, width = plot_width, height = plot_height, dpi = 1200)

message("Supplementary Figure 3 generation complete.")