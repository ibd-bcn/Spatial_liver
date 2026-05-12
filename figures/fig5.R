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

# Output Paths (Relative to the repository root)
OUT_DIR <- "figures/outs/"

# Ensure output directory exists
if (!dir.exists(OUT_DIR)) {
  dir.create(OUT_DIR, recursive = TRUE)
}

# Define Target Cells
hepatocytes <- c("Hepatocyte 1", "Hepatocyte 2", "Hepatocyte 3", "Hepatocyte 5", "Hepatocyte 6")

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
# Figure 5A: Explant Interactions (Slide 3) -> Ag General
# ==============================================================================

# Filter for Slide 3 specific FOVs and map Antigen Status
slide3_fovs <- c(1, 4, 11, 16, 24, 25, 13, 2, 3, 5, 6, 12, 14, 23)

cut_int_5a <- all_int %>%
  filter(tissue == "Slide_3", fov %in% slide3_fovs, likelihood > 0.5)

# Map Ag_general from Seurat to SCOTIA results
ag_map <- seu$Ag_general
names(ag_map) <- seu$cell_names
cut_int_5a$ag_general <- ag_map[cut_int_5a$id_source]

# Drop NAs to prevent math errors
cut_int_5a <- cut_int_5a %>% filter(!is.na(ag_general))

# Calculate Interaction Frequencies
df_hep_5a <- calc_interaction_freq(cut_int_5a, hepatocytes, "ag_general", "Neg", "Hepatocytes All")
df_kc1_5a <- calc_interaction_freq(cut_int_5a, c(hepatocytes, "KC1"), "ag_general", "Neg", "Hep + KC1")
df_kc2_5a <- calc_interaction_freq(cut_int_5a, c(hepatocytes, "KC2"), "ag_general", "Neg", "Hep + KC2")

all_data_5a <- bind_rows(df_hep_5a, df_kc1_5a, df_kc2_5a)

# Factor Ordering
all_data_5a$ag_general <- factor(all_data_5a$ag_general, levels = c("Neg", "S+"))
all_data_5a$Condition  <- factor(all_data_5a$Condition, levels = c("Hepatocytes All", "Hep + KC1", "Hep + KC2"))

# Plot 5A
p_5a <- ggplot(all_data_5a, aes(x = Condition, y = norm, fill = ag_general)) +
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

png(filename = file.path(OUT_DIR, "figure5_A.png"), width = 10, height = 6, units = "in", res = 400)
print(p_5a)
dev.off()

# ==============================================================================
# Figure 5B
# Note: Can be found in 7.localcomposition.py
# ==============================================================================

# ==============================================================================
# Figure 5C: Explant Interactions (Slide 2) -> Type D
# ==============================================================================

# Filter for Slide 2 specific FOVs and map Type D Status
slide2_fovs <- c(1:12)

cut_int_5c <- all_int %>%
  filter(tissue == "Slide_2", fov %in% slide2_fovs, likelihood > 0.5)

# Map Type_D from Seurat to SCOTIA results
typed_map <- seu$Type_D
names(typed_map) <- seu$cell_names
cut_int_5c$Type_D <- typed_map[cut_int_5c$id_source]

# Drop NAs to prevent math errors
cut_int_5c <- cut_int_5c %>% filter(!is.na(Type_D) & Type_D != "none")

# Calculate Interaction Frequencies
df_hep_5c <- calc_interaction_freq(cut_int_5c, hepatocytes, "Type_D", "-", "Hepatocytes All")
df_kc1_5c <- calc_interaction_freq(cut_int_5c, c(hepatocytes, "KC1"), "Type_D", "-", "Hep + KC1")
df_kc2_5c <- calc_interaction_freq(cut_int_5c, c(hepatocytes, "KC2"), "Type_D", "-", "Hep + KC2")

all_data_5c <- bind_rows(df_hep_5c, df_kc1_5c, df_kc2_5c)

# Factor Ordering
all_data_5c$Type_D    <- factor(all_data_5c$Type_D, levels = c("-", "+"))
all_data_5c$Condition <- factor(all_data_5c$Condition, levels = c("Hepatocytes All", "Hep + KC1", "Hep + KC2"))

# Plot 5C
p_5c <- ggplot(all_data_5c, aes(x = Condition, y = norm, fill = Type_D)) +
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

png(filename = file.path(OUT_DIR, "figure5_C.png"), width = 10, height = 6, units = "in", res = 400)
print(p_5c)
dev.off()

# ==============================================================================
# Figure 5D
# Note: Can be found in 7.localcomposition.py
# ==============================================================================