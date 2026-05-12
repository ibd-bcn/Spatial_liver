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

# ==============================================================================
# Figure 4A, 4B, 4C, 4D (Pending specific code additions)
# ==============================================================================

# ==============================================================================
# Figure 4E: SCOTIA Cell-to-Cell Interactions
# ==============================================================================
message("Generating Figure 4E (SCOTIA Interactions)...")

all_int <- read_csv(SCOTIA_PATH, show_col_types = FALSE)

# Map etiology from Seurat metadata to SCOTIA results
logy <- meta$etiology
names(logy) <- meta$cell_names
all_int$etiology <- logy[all_int$id_source]

hepatocytes <- c("Hepatocyte 1", "Hepatocyte 2", "Hepatocyte 3", "Hepatocyte 5", "Hepatocyte 6")

# Filter for relevant etiologies and high likelihood interactions
cut_int <- all_int %>% 
  filter(etiology %in% c("HDV RNA+", "HBV", "HC"), likelihood > 0.5)

# Calculate totals per etiology
int_table_all <- as.data.frame(table(cut_int$etiology))
names(int_table_all) <- c("etiology", "all_int")

# --- 1. Hepatocytes (All) ---
hep_all <- cut_int %>% 
  filter(refined_receptor %in% hepatocytes | refined_source %in% hepatocytes)
int_hep_all <- as.data.frame(table(hep_all$etiology))

df_hep <- int_table_all
df_hep$hep_all <- int_hep_all$Freq
df_hep$perc <- (df_hep$hep_all / df_hep$all_int) * 100
df_hep$norm <- df_hep$perc / df_hep$perc[df_hep$etiology == "HC"] # Safely normalize to HC
df_hep$Condition <- "Hepatocytes All"

# --- 2. Hepatocytes with KC2 ---
hep_kc2 <- cut_int %>% 
  filter(refined_receptor %in% c(hepatocytes, "KC2") & refined_source %in% c(hepatocytes, "KC2"))
int_hep_kc2 <- as.data.frame(table(hep_kc2$etiology))

df_kc2 <- int_table_all
df_kc2$hep_all <- int_hep_kc2$Freq
df_kc2$perc <- (df_kc2$hep_all / df_kc2$all_int) * 100
df_kc2$norm <- df_kc2$perc / df_kc2$perc[df_kc2$etiology == "HC"] # Safely normalize to HC
df_kc2$Condition <- "Hep + KC2"

# Combine data for plotting
all_data <- bind_rows(df_hep, df_kc2)

# Ensure factor ordering
all_data$etiology <- factor(all_data$etiology, levels = c("HC", "HBV", "HDV RNA+"))
all_data$Condition <- factor(all_data$Condition, levels = c("Hepatocytes All", "Hep + KC2"))

# Plot 4E
p_4e <- ggplot(all_data, aes(x = Condition, y = norm, fill = etiology)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black", width = 0.7) +
  scale_fill_manual(values = etiology_colors) + 
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1),
    axis.ticks.length = unit(0.2, "cm"),
    axis.text.x = element_text(size = 12, face = "bold", color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 12, face = "bold", color = "black"),
    axis.title.x = element_blank(), 
    legend.position = "top"
  ) +
  labs(y = "Normalized Frequency (Ref: HC)", fill = "Etiology")

png(filename = file.path(OUT_DIR, "figure4_E.png"), width = 8, height = 6, units = "in", res = 400)
print(p_4e)
dev.off()

# ==============================================================================
# Figure 4F: Cell Co-occurrence Enrichment (Hepatocytes & KC)
# ==============================================================================
message("Generating Figure 4F (Co-occurrence Enrichment)...")

all_cooccur <- read_csv(COOCCUR_PATH, show_col_types = FALSE)

# Filter for Hepatocyte interactions with KC1 and KC2
all_kc1 <- all_cooccur %>% 
  filter(from == "Hepatocyte", to %in% c("KC1", "KC2"), etiology != "HDV RNA-")

# Set factor levels for correct plotting order
all_kc1$interval_numeric <- factor(all_kc1$bin, levels = sort(unique(all_kc1$bin)))
all_kc1$etiology <- factor(all_kc1$etiology, levels = c("HC", "HBV", "HDV RNA+"))

# Plot 4F
p_4f <- ggplot(all_kc1, aes(x = bin, y = enrichment, group = to, color = to)) +
  geom_smooth(alpha = 0.05, linewidth = 1.5, method = "loess", formula = y ~ x) + 
  geom_hline(yintercept = 0, color = "black", linewidth = 1, linetype = "dashed") +
  facet_wrap(~ etiology) + 
  scale_color_manual(values = refined_col) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10), labels = NULL) +
  theme_linedraw() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
    axis.line = element_blank(),
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_line(linewidth = 1.5, color = "black"),
    axis.ticks.y = element_line(linewidth = 1.5, color = "black"),
    axis.ticks.length = unit(6, "pt"),
    text = element_text(size = 20, face = "bold")
  ) +
  labs(x = NULL, y = NULL)

png(filename = file.path(OUT_DIR, "figure4_F_coocurrence_hep_kc1.png"), width = 14, height = 6, units = "in", res = 1200)
print(p_4f)
dev.off()

message("Figure 4 generation complete.")