# ==============================================================================
# Script: supfig6.R
# Description: Generates panels for Supplementary Figure 6 (Slingshot Trajectories & Heatmap)
# ==============================================================================

suppressPackageStartupMessages({
  library(slingshot)
  library(tradeSeq)
  library(viridis)
  library(ggplot2)
  library(dplyr)
  library(pheatmap)
  library(RColorBrewer)
})

options(bitmapType = "cairo")

# ==============================================================================
# Configuration & Paths
# Update these paths according to your local environment
# ==============================================================================
SCE_SLINGSHOT_PATH <- "/path/to/Post_neigh/sce_slingshot.rds"
SCE_TRADESEQ_PATH  <- "/path/to/Post_neigh/sce_tradeseq.rds"

OUT_DIR <- "figures/outs/"

if (!dir.exists(OUT_DIR)) {
  dir.create(OUT_DIR, recursive = TRUE)
}

# ==============================================================================
# 1. Load Data
# ==============================================================================
sce_slingshot <- readRDS(SCE_SLINGSHOT_PATH)
sce_tradeseq  <- readRDS(SCE_TRADESEQ_PATH)

# ==============================================================================
# Helper Function for Trajectory Plotting
# ==============================================================================
plot_trajectory <- function(sce_obj, traj_num, out_filename) {
  
  plot_df <- data.frame(
    UMAP_1 = reducedDim(sce_obj, 'UMAP')[, 1],
    UMAP_2 = reducedDim(sce_obj, 'UMAP')[, 2],
    pseudotime = slingPseudotime(sce_obj)[, traj_num]
  )
  
  curves <- slingCurves(sce_obj)
  curve_df <- as.data.frame(curves[[traj_num]]$s)
  colnames(curve_df) <- c("UMAP_1", "UMAP_2")
  
  p <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = pseudotime), size = 0.01, alpha = 0.8) +
    geom_path(data = curve_df, aes(x = UMAP_1, y = UMAP_2), linewidth = 1.2, color = "black") +
    scale_color_viridis_c(option = "magma", na.value = "lightgrey") +
    labs(
      title = "Slingshot Pseudotime on UMAP",
      subtitle = paste0("Trajectory ", traj_num, " with Principal Curve"),
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
  
  png(filename = out_filename, width = 9.2, height = 6, units = "in", res = 800)
  print(p)
  dev.off()
}

# ==============================================================================
# Supplementary Figure 6A: Slingshot Trajectories (1-4)
# ==============================================================================
# Loop through trajectories 1 to 4 and generate plots
for (traj_idx in 1:4) {
  out_path <- file.path(OUT_DIR, paste0("supfig6A_traj", traj_idx, ".png"))
  plot_trajectory(sce_slingshot, traj_idx, out_path)
}

# ==============================================================================
# Supplementary Figure 6B: TradeSeq Expression Heatmap
# ==============================================================================
selected_clusters <- c(6, 2, 0, 1)
sel <- sce_tradeseq$leiden0.2 %in% selected_clusters

my_genes <- c(
  "SERPINA1", "XBP1", "FASN", "CEACAM1", "SERPINA3", "APOA1", 
  "COL16A1", "CXCL17", "IL7R", "NDRG1", "ICAM3", "INSR", 
  "ITGB1", "IFITM1", "VTN", "BIRC3", "TWIST1", "IL1RAP", 
  "PTEN", "IL1RN", "MET", "LINC01781", "FZD3", "CLEC14A", 
  "SMO", "EFNA1", "SRSF2", "ITGA1", "PSD3", "ZBTB16", 
  "PFN1", "NEAT1", "FN1"
)

# Extract and format count data
heatdata <- assays(sce_tradeseq)$counts[my_genes, sel]
heatclus <- factor(sce_tradeseq$leiden0.2[sel], levels = selected_clusters)

# Calculate cluster-wise average expression
heatdata_avg <- sapply(levels(heatclus), function(cl){
  rowMeans(heatdata[, heatclus == cl, drop = FALSE])
})

# Setup annotations, colors, and breaks
annotation_col <- data.frame(Cluster = factor(levels(heatclus), levels = selected_clusters))
rownames(annotation_col) <- levels(heatclus)

my_colors <- colorRampPalette(c("#053061", "white", "#67001F"))(100)
my_breaks <- seq(-2, 2, length.out = 101)

# Plot and save heatmap
png(filename = file.path(OUT_DIR, "supfig6B_heatmap.png"), width = 5, height = 8, units = "in", res = 600)
pheatmap(
  mat = log1p(heatdata_avg),       
  scale = "row",                   
  cluster_cols = FALSE,            
  cluster_rows = TRUE,             
  treeheight_row = 0,              
  color = my_colors,               
  breaks = my_breaks,              
  annotation_col = annotation_col, 
  show_rownames = TRUE,            
  main = "Average Expression by Cluster: 6 -> 2 -> 0 -> 1",
  border_color = NA
)
dev.off()

# ==============================================================================
# Supplementary Figure 6C: nCounter boxplots
# ==============================================================================

expression_path <- "/path/to/ALL_normalized_data with clinical variables_Nanostring.xlsx"

output_dir <- "/path/to/ncounter/"
boxplot_dir <- file.path(output_dir, "boxplots")

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  boxplot_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==============================================================================
# 2. Load and prepare NanoString nCounter data
# ==============================================================================

df <- read_excel(expression_path)

# Keep expression data
df <- df[1:39, ]
df <- as.data.frame(df)

# First row contains column names
colnames(df) <- df[1, ]

# Remove header row
df <- df[2:39, ]

# First column contains sample IDs
rownames(df) <- df[, 1]

# Keep expression columns
df <- df[, 2:783]

# Transpose:
# rows = genes
# columns = samples
df <- t(df)

# Remove non-gene / technical rows
df <- df[10:782, , drop = FALSE]

# Convert expression values to numeric
expr_matrix <- df
mode(expr_matrix) <- "numeric"
# Remove NanoString target suffix from gene names
rownames(expr_matrix) <- str_remove(
  rownames(expr_matrix),
  "-mRNA$"
)

# Ensure gene names are unique
rownames(expr_matrix) <- make.unique(
  rownames(expr_matrix)
)

# ==============================================================================
# 3. Create sample metadata
# ==============================================================================

sample_names <- colnames(expr_matrix)

condition <- ifelse(
  str_starts(sample_names, "D"), "delta",
  ifelse(
    str_starts(sample_names, "C"), "control",
    ifelse(
      str_starts(sample_names, "N"), "negative",
      ifelse(
        str_starts(sample_names, "B"), "vb",
        NA
      )
    )
  )
)

metadata <- data.frame(
  sample = sample_names,
  condition = condition
)

rownames(metadata) <- metadata$sample


# Rename conditions for figures
metadata$condition2 <- recode(
  metadata$condition,
  control  = "HC",
  vb       = "HBV",
  delta    = "HDV RNA+",
  negative = "HDV RNA-"
)

metadata$condition2 <- factor(
  metadata$condition2,
  levels = c(
    "HC",
    "HBV",
    "HDV RNA+",
    "HDV RNA-"
  )
)

# Ensure expression matrix and metadata have the same sample order
metadata <- metadata[colnames(expr_matrix), , drop = FALSE]


# 4 colors
condition_palette <- c(
  "HC"       = "#fe4a49",
  "HBV"      = "#2ab7ca",
  "HDV RNA+" = "#fed766",
  "HDV RNA-" = "#D866FE"
)


# ==============================================================================
# 5. Boxplots for selected genes
# ==============================================================================

plot_gene_boxplot <- function(gene) {
  
  # Match gene names ignoring upper/lower case
  gene_match <- rownames(expr_matrix)[
    toupper(rownames(expr_matrix)) == toupper(gene)
  ]
  
  if (length(gene_match) == 0) {
    stop(
      "Gene not found in expression matrix: ",
      gene
    )
  }
  
  gene_match <- gene_match[1]
  
  df_plot <- data.frame(
    expr = as.numeric(expr_matrix[gene_match, ]),
    Condition = metadata$condition2,
    Sample = rownames(metadata)
  )
  
  df_plot$Condition <- factor(
    df_plot$Condition,
    levels = c(
      "HC",
      "HBV",
      "HDV RNA+",
      "HDV RNA-"
    )
  )
  
  # Remove missing values if present
  df_plot <- df_plot %>%
    filter(
      !is.na(expr),
      !is.na(Condition)
    )
  
  
  # --------------------------------------------------------------------------
  # Pairwise Wilcoxon tests: all groups vs all groups
  # Benjamini-Hochberg correction across the 6 comparisons
  # --------------------------------------------------------------------------
  
  stat_test <- compare_means(
    expr ~ Condition,
    data = df_plot,
    method = "wilcox.test",
    p.adjust.method = "BH"
  )
  
  # Significance labels based on adjusted p-values
  stat_test$p.adj.signif <- cut(
    stat_test$p.adj,
    breaks = c(
      -Inf,
      0.0001,
      0.001,
      0.01,
      0.05,
      Inf
    ),
    labels = c(
      "****",
      "***",
      "**",
      "*",
      "ns"
    )
  )
  
  
  # --------------------------------------------------------------------------
  # Position significance brackets
  # --------------------------------------------------------------------------
  
  y_min <- min(df_plot$expr, na.rm = TRUE)
  y_max <- max(df_plot$expr, na.rm = TRUE)
  y_range <- y_max - y_min
  
  # Avoid problems if expression range is zero
  if (y_range == 0) {
    y_range <- abs(y_max)
    
    if (y_range == 0) {
      y_range <- 1
    }
  }
  
  stat_test$y.position <- y_max +
    seq(
      0.10,
      0.60,
      length.out = nrow(stat_test)
    ) * y_range
  
  
  # --------------------------------------------------------------------------
  # Plot
  # --------------------------------------------------------------------------
  
  p <- ggplot(
    df_plot,
    aes(
      x = Condition,
      y = expr,
      fill = Condition
    )
  ) +
    geom_boxplot(
      width = 0.55,
      alpha = 0.9,
      outlier.shape = NA,
      color = "black",
      linewidth = 1.2
    ) +
    scale_fill_manual(
      values = condition_palette
    ) +
    scale_y_continuous(
      expand = expansion(
        mult = c(0.05, 0.45)
      )
    ) +
    stat_pvalue_manual(
      stat_test,
      label = "p.adj.signif",
      xmin = "group1",
      xmax = "group2",
      y.position = "y.position",
      tip.length = 0.01,
      bracket.size = 0.8,
      size = 6,
      hide.ns = FALSE
    ) +
    labs(
      title = toupper(gene),
      x = "",
      y = "Normalized expression"
    ) +
    theme_classic(
      base_size = 28,
      base_family = "Arial"
    ) +
    theme(
      panel.grid = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      plot.title = element_text(
        size = 14,
        face = "bold",
        hjust = 0.5
      ),
      
      axis.title = element_text(
        size = 12
      ),
      
      axis.text = element_text(
        size = 10
      ),
      
      axis.text.x = element_text(
        angle = 0,
        hjust = 0.5
      ),
      
      axis.line = element_line(
        linewidth = 1.5
      ),
      
      axis.ticks = element_line(
        linewidth = 1.5
      ),
      
      axis.ticks.length = unit(
        10,
        "pt"
      ),
      
      legend.position = "none"
    )
  
  return(p)
}


# ==============================================================================
# 6. Generate boxplots
# ==============================================================================

genes_boxplot <- c(
  "XBP1",
  "ICAM3",
  "IL7R",
  "SERPINA1"
)

boxplot_list <- lapply(
  genes_boxplot,
  plot_gene_boxplot
)

names(boxplot_list) <- genes_boxplot


# ==============================================================================
# 7. Save boxplots
# ==============================================================================

save_boxplot <- function(plot_object, gene_name) {
  
  output_file <- file.path(
    boxplot_dir,
    paste0(
      tolower(gene_name),
      "_box.png"
    )
  )
  
  png(
    filename = output_file,
    width = 6,
    height = 7,
    units = "in",
    res = 1200
  )
  
  print(plot_object)
  
  dev.off()
  
  message(
    "Saved: ",
    normalizePath(
      output_file,
      mustWork = FALSE
    )
  )
}

purrr::iwalk(
  boxplot_list,
  ~ save_boxplot(.x, .y)
)
