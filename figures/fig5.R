# ==============================================================================
# Script: fig5.R
# Description: Generates panels for Figure 5 (Explant interaction frequencies)
# Note: Figures 5B and 5D are generated via 7.localcomposition.py
# ==============================================================================

# Load necessary libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
})

options(bitmapType = "cairo")

# ==============================================================================
# Configuration & Paths
# Update these paths according to your local environment
# ==============================================================================
# Input Paths
SEURAT_PATH <- "/path/to/Objects/seurats_annotated.RDS"
SCOTIA_PATH <- "/path/to/SCOTIA/Results/all_int.csv"
COOCCUR_PATH <- "/path/to/Celltype_enrichment/enrichment_files/all.csv"

# Output Paths (Relative to the repository root)
OUT_DIR <- "figures/outs/"

# Ensure output directory exists
if (!dir.exists(OUT_DIR)) {
  dir.create(OUT_DIR, recursive = TRUE)
}

# Define Target Cells
hepatocytes <- c("Hepatocyte 1", "Hepatocyte 2", "Hepatocyte 3", "Hepatocyte 5", "Hepatocyte 6")

# Define Palettes (Required for independent execution)
etiology_colors <- c("HC" = "#fe4a49", "HBV" = "#2ab7ca", "HDV RNA+" = "#fed766")

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
# 1. Load Data
# ==============================================================================
seu <- readRDS(SEURAT_PATH)
all_int <- read_csv(SCOTIA_PATH, show_col_types = FALSE)

# ==============================================================================
# 2. Reusable Helper Function
# ==============================================================================
# This function safely calculates frequencies and normalizes against a control
calc_interaction_freq <- function(int_df, target_cells, group_col, ref_level, condition_name) {
  
  # Total interactions per group
  total_counts <- int_df %>% 
    group_by(!!sym(group_col)) %>% 
    summarise(all_int = n(), .groups = "drop")
  
  # Subset for specific receptor/source interactions
  sub_int <- int_df %>% 
    filter(refined_receptor %in% target_cells | refined_source %in% target_cells)
  
  sub_counts <- sub_int %>% 
    group_by(!!sym(group_col)) %>% 
    summarise(target_int = n(), .groups = "drop")
  
  # Merge, calculate percentages, and normalize against the reference level
  res <- total_counts %>% 
    left_join(sub_counts, by = group_col) %>%
    mutate(
      target_int = replace_na(target_int, 0),
      perc = (target_int / all_int) * 100
    )
  
  # Safe Normalization
  ref_val <- res$perc[res[[group_col]] == ref_level]
  if(length(ref_val) == 0 || ref_val == 0) stop(paste("Reference level", ref_level, "not found or is zero."))
  
  res$norm <- res$perc / ref_val
  res$Condition <- condition_name
  
  return(res)
}

# ==============================================================================
# Figure 5A: SCOTIA Cell-to-Cell Interactions
# ==============================================================================
message("Generating Figure 5A (SCOTIA Interactions)...")

# Map etiology from Seurat metadata to SCOTIA results
logy <- seu@meta.data$etiology
names(logy) <- seu@meta.data$cell_names
all_int$etiology <- logy[all_int$id_source]

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

# --- 2. Hepatocytes with h-KC ---
hep_kc2 <- cut_int %>% 
  filter(refined_receptor %in% c(hepatocytes, "h-KC") & refined_source %in% c(hepatocytes, "h-KC"))
int_hep_kc2 <- as.data.frame(table(hep_kc2$etiology))

df_kc2 <- int_table_all
df_kc2$hep_all <- int_hep_kc2$Freq
df_kc2$perc <- (df_kc2$hep_all / df_kc2$all_int) * 100
df_kc2$norm <- df_kc2$perc / df_kc2$perc[df_kc2$etiology == "HC"] # Safely normalize to HC
df_kc2$Condition <- "Hep + h-KC"

# Combine data for plotting
all_data <- bind_rows(df_hep, df_kc2)

# Ensure factor ordering
all_data$etiology <- factor(all_data$etiology, levels = c("HC", "HBV", "HDV RNA+"))
all_data$Condition <- factor(all_data$Condition, levels = c("Hepatocytes All", "Hep + h-KC"))

# Plot 5A
p_5a <- ggplot(all_data, aes(x = Condition, y = norm, fill = etiology)) +
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

png(filename = file.path(OUT_DIR, "fig5a.png"), width = 8, height = 6, units = "in", res = 400)
print(p_5a)
dev.off()

# ==============================================================================
# Figure 5B
# Note: Panel generated via 7.localcomposition.py
# ==============================================================================

# ==============================================================================
# Figure 5C: Cell Co-occurrence Enrichment (Hepatocytes & KC)
# ==============================================================================
message("Generating Figure 5C (Co-occurrence Enrichment)...")

all_cooccur <- read_csv(COOCCUR_PATH, show_col_types = FALSE)

# Filter for Hepatocyte interactions with i-KC and h-KC
all_kc1 <- all_cooccur %>% 
  filter(from == "Hepatocyte", to %in% c("i-KC", "h-KC"), etiology != "HDV RNA-")

# Set factor levels for correct plotting order
all_kc1$interval_numeric <- factor(all_kc1$bin, levels = sort(unique(all_kc1$bin)))
all_kc1$etiology <- factor(all_kc1$etiology, levels = c("HC", "HBV", "HDV RNA+"))

# Plot 5C
p_5c <- ggplot(all_kc1, aes(x = bin, y = enrichment, group = to, color = to)) +
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

png(filename = file.path(OUT_DIR, "fig5c.png"), width = 14, height = 6, units = "in", res = 1200)
print(p_5c)
dev.off()

# ==============================================================================
# Figure 5D
# Note: Panel generated via 7.localcomposition.py
# ==============================================================================

# ==============================================================================
# Figure 5E: Explant Interactions (Slide 3) -> Ag General
# ==============================================================================
message("Generating Figure 5E (Explant Interactions Ag General)...")

# Filter for Slide 3 specific FOVs and map Antigen Status
slide3_fovs <- c(1, 4, 11, 16, 24, 25, 13, 2, 3, 5, 6, 12, 14, 23)

cut_int_5e <- all_int %>%
  filter(tissue == "Slide_3", fov %in% slide3_fovs, likelihood > 0.5)

# Map Ag_general from Seurat to SCOTIA results
ag_map <- seu@meta.data$Ag_general
names(ag_map) <- seu@meta.data$cell_names
cut_int_5e$ag_general <- ag_map[cut_int_5e$id_source]

# Drop NAs to prevent math errors
cut_int_5e <- cut_int_5e %>% filter(!is.na(ag_general))

# Calculate Interaction Frequencies
df_hep_5e <- calc_interaction_freq(cut_int_5e, hepatocytes, "ag_general", "Neg", "Hepatocytes All")
df_kc1_5e <- calc_interaction_freq(cut_int_5e, c(hepatocytes, "i-KC"), "ag_general", "Neg", "Hep + i-KC")
df_kc2_5e <- calc_interaction_freq(cut_int_5e, c(hepatocytes, "h-KC"), "ag_general", "Neg", "Hep + h-KC")

all_data_5e <- bind_rows(df_hep_5e, df_kc1_5e, df_kc2_5e)

# Factor Ordering
all_data_5e$ag_general <- factor(all_data_5e$ag_general, levels = c("Neg", "S+"))
all_data_5e$Condition  <- factor(all_data_5e$Condition, levels = c("Hepatocytes All", "Hep + i-KC", "Hep + h-KC"))

# Plot 5E
p_5e <- ggplot(all_data_5e, aes(x = Condition, y = norm, fill = ag_general)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black", width = 0.7) +
  scale_fill_manual(values = c("Neg" = "#02876f", "S+" = "#f07801")) + 
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
  labs(y = "Normalized Frequency (Ref: Neg)", fill = "Antigen Status")

png(filename = file.path(OUT_DIR, "fig5e.png"), width = 10, height = 6, units = "in", res = 400)
print(p_5e)
dev.off()


# ==============================================================================
# Figure 5F: Explant Interactions (Slide 2) -> Type D
# ==============================================================================
message("Generating Figure 5F (Explant Interactions Type D)...")

# Filter for Slide 2 specific FOVs and map Type D Status
slide2_fovs <- c(1:12)

cut_int_5f <- all_int %>%
  filter(tissue == "Slide_2", fov %in% slide2_fovs, likelihood > 0.5)

# Map Type_D from Seurat to SCOTIA results
typed_map <- seu@meta.data$Type_D
names(typed_map) <- seu@meta.data$cell_names
cut_int_5f$Type_D <- typed_map[cut_int_5f$id_source]

# Drop NAs to prevent math errors
cut_int_5f <- cut_int_5f %>% filter(!is.na(Type_D) & Type_D != "none")

# Calculate Interaction Frequencies
df_hep_5f <- calc_interaction_freq(cut_int_5f, hepatocytes, "Type_D", "-", "Hepatocytes All")
df_kc1_5f <- calc_interaction_freq(cut_int_5f, c(hepatocytes, "i-KC"), "Type_D", "-", "Hep + i-KC")
df_kc2_5f <- calc_interaction_freq(cut_int_5f, c(hepatocytes, "h-KC"), "Type_D", "-", "Hep + h-KC")

all_data_5f <- bind_rows(df_hep_5f, df_kc1_5f, df_kc2_5f)

# Factor Ordering
all_data_5f$Type_D    <- factor(all_data_5f$Type_D, levels = c("-", "+"))
all_data_5f$Condition <- factor(all_data_5f$Condition, levels = c("Hepatocytes All", "Hep + i-KC", "Hep + h-KC"))

# Plot 5F
p_5f <- ggplot(all_data_5f, aes(x = Condition, y = norm, fill = Type_D)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black", width = 0.7) +
  scale_fill_manual(values = c("-" = "#f9b80d", "+" = "#c51f05")) + 
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
  labs(y = "Normalized Frequency (Ref: -)", fill = "Type D Status")

png(filename = file.path(OUT_DIR, "fig5f.png"), width = 10, height = 6, units = "in", res = 400)
print(p_5f)
dev.off()

message("Figure 5 generation complete.")