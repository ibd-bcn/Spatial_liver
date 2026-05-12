#Read allint
library(readr)
all_int <- read_csv("Maria_v2/Post_analysis/SCOTIA/Results/all_int.csv")
library(Seurat)
library(paletteer)
library(plyr)
library(dplyr)
library(Seurat)
library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
library(ggrepel)
library(tidyr)
library(circlize)
library(devtools)
library(readr)
library(ggplot2)
library(msigdbr)
library(org.Hs.eg.db)
load_all("/home/mmoro/SPATIAL/Spatial_Package/")
options(bitmapType = "cairo")
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)


seu <- readRDS("/home/mmoro/SPATIAL/Maria_CosMx/Maria_v2/Post_analysis/Selecting_regions/seu.RDS")
seu <- JoinLayers(seu)
seu <- NormalizeData(seu)
subset_type <- c(
  "hepatocytes" = "Hepatocytes",
  "myeloids" = "Myeloid_cells",
  "plasmas" = "B_cell_lineage",
  "tcells" = "T_cell_lineage",
  "non-parenquimal" = "Non_parenchymal",
  "IM-hepatocytes" = "Im_hepatocytes"
)


grouped_anot <- c(
  "Non_parenchymal" = "Hepatocytes_Non_Parenchymal",
  "Hepatocytes" = "Hepatocytes_Non_Parenchymal",
  "Myeloid_cells" = "Immune_cells",
  "B_cell_lineage" = "Immune_cells",
  "T_cell_lineage" = "Immune_cells",
  "Im_hepatocytes" = "Hepatocytes_Non_Parenchymal"
)


cell_types <- c(
  "Gamma-delta T cells"       = "Gamma_delta_T_cells",
  "NK cells"                  = "NK_cells",
  "NKT cells"                 = "NKT_cells",
  "Regulatory T cells"        = "Regulatory_T_cells",
  "Tem/Trm cytotoxic T cells" = "Tem_Trm_cytotoxic_T_cells",
  "Effector helper T cells"   = "Effector_helper_T_cells",
  "Naive T cells"             = "Naive_T_cells",
  "Memory B cells"            = "Memory_B_cells",
  "Plasma cells"              = "Plasma_cells",
  "Naive B cells"             = "Naive_B_cells",
  "Monocytes"                 = "Monocytes",
  "KC1"                       = "KC1",
  "M2-LYVE1"                  = "M2_LYVE1",
  "KC2"                       = "KC2",
  "DCs CD1C"                  = "DCs_CD1C",
  "M1"                        = "M1",
  "Endothelial cells 2"       = "Endothelial_cells_2",
  "Fibroblasts"               = "Fibroblasts",
  "Endothelial cells 1"       = "Endothelial_cells_1",
  "Smooth muscle cells"       = "Smooth_muscle_cells",
  "Endothelial cells 4"       = "Endothelial_cells_4",
  "Endothelial cells 3"       = "Endothelial_cells_3",
  "Hepatocyte 2"              = "Hepatocyte_2",
  "Hepatocyte 1"              = "Hepatocyte_1",
  "Hepatocyte 6"              = "Hepatocyte_6",
  "Hepatocyte 3"              = "Hepatocyte_3",
  "Hepatocyte 5"              = "Hepatocyte_5",
  "Hepatocyte 4"              = "Hepatocyte_4",
  "Cholangiocytes"            = "Cholangiocytes",
  "Cycling B lineage cells"   = "Cycling_B_lineage_cells"
)


new_anot <- c(
  "Gamma_delta_T_cells"       = "T_cells",
  "NK_cells"                  = "T_cells",
  "NKT_cells"                 = "T_cells",
  "Regulatory_T_cells"        = "T_cells",
  "Tem_Trm_cytotoxic_T_cells" = "T_cells",
  "Effector_helper_T_cells"   = "T_cells",
  "Naive_T_cells"             = "T_cells",
  "Memory_B_cells"            = "B_cells",
  "Plasma_cells"              = "B_cells",
  "Naive_B_cells"             = "B_cells",
  "Monocytes"                 = "Myeloid_cells",
  "KC1"                       = "Kupffer_cells",
  "M2_LYVE1"                  = "Myeloid_cells",
  "KC2"                       = "Kupffer_cells",
  "DCs_CD1C"                  = "Myeloid_cells",
  "M1"                        = "Myeloid_cells",
  "Endothelial_cells_2"       = "Non_parenchymal",
  "Fibroblasts"               = "Non_parenchymal",
  "Endothelial_cells_1"       = "Non_parenchymal",
  "Smooth_muscle_cells"       = "Non_parenchymal",
  "Endothelial_cells_4"       = "Non_parenchymal",
  "Endothelial_cells_3"       = "Non_parenchymal",
  "Hepatocyte_2"              = "Hepatocytes",
  "Hepatocyte_1"              = "Hepatocytes",
  "Hepatocyte_6"              = "Hepatocytes",
  "Hepatocyte_3"              = "Hepatocytes",
  "Hepatocyte_5"              = "Hepatocytes",
  "Hepatocyte_4"              = "Hepatocytes",
  "Cholangiocytes"            = "Non_parenchymal",
  "Cycling_B_lineage_cells"   = "B_cells"
)


#Red FOVs
seu@meta.data$red <- FALSE
seu@meta.data[seu@meta.data$tissue == "Slide_2" & seu@meta.data$fov %in% c(21,22,23,24),]$red <- TRUE
#Orange FOVs
seu@meta.data$orange <- FALSE
seu@meta.data[seu@meta.data$tissue == "Slide_1" & seu@meta.data$fov %in% c(25),]$orange <- TRUE
seu@meta.data[seu@meta.data$tissue == "Slide_2" & seu@meta.data$fov %in% c(2,3,4,6,7,10),]$orange <- TRUE
seu@meta.data[seu@meta.data$tissue == "Slide_3" & seu@meta.data$fov %in% c(6,12,14,23),]$orange <- TRUE
#Antigen level
seu@meta.data$antigen_level <- NA
#neg
seu@meta.data[seu@meta.data$tissue == "Slide_1" & seu@meta.data$fov %in% c(1,2,3,4,5,6,7,8,25), "antigen_level"] <- "neg"
seu@meta.data[seu@meta.data$tissue == "Slide_3" & seu@meta.data$fov %in% c(1,4,11,16,24,25), "antigen_level"] <- "neg"
#'S_low'
seu@meta.data[seu@meta.data$tissue == "Slide_2" & seu@meta.data$fov %in% c(13,14,15,17,19,20), "antigen_level"] <- "S_low"
seu@meta.data[seu@meta.data$tissue == "Slide_3" & seu@meta.data$fov %in% c(13), "antigen_level"] <- "S_low"
# 'S_high' 
seu@meta.data[seu@meta.data$tissue == "Slide_2" & seu@meta.data$fov %in% c(16,18), "antigen_level"] <- "S_high"
seu@meta.data[seu@meta.data$tissue == "Slide_3" & seu@meta.data$fov %in% c(2,3,5,6,12,14,23), "antigen_level"] <- "S_high"
# 'S+/-'
seu@meta.data[seu@meta.data$tissue == "Slide_3" & seu@meta.data$fov %in% c(7,8,9,10,15,17,18,19,20,21,22), "antigen_level"] <- "S+/-"
#S&D_low' 
seu@meta.data[seu@meta.data$tissue == "Slide_2" & seu@meta.data$fov %in% c(1,5,6,7,9,11,12), "antigen_level"] <- "S&D_low"
#'S&D' 
seu@meta.data[seu@meta.data$tissue == "Slide_2" & seu@meta.data$fov %in% c(2,3,4,8,10), "antigen_level"] <- "S&D"

#etiology
seu@meta.data$etiology <- NA
# HBV
seu@meta.data[seu@meta.data$tissue == "Slide_2" & seu@meta.data$fov %in% c(13,14,15,16,17,18,19,20,21), "etiology"] <- "HBV"
seu@meta.data[seu@meta.data$tissue == "Slide_3" & seu@meta.data$fov %in% c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,23,24,25), "etiology"] <- "HBV"
# HDV RNA+
seu@meta.data[seu@meta.data$tissue == "Slide_2" & seu@meta.data$fov %in% c(1,2,3,4,5,6,7,8,9,10,11,12,22,23,24), "etiology"] <- "HDV RNA+"
# HDV RNA-
seu@meta.data[seu@meta.data$tissue == "Slide_1" & seu@meta.data$fov %in% c(1,2,3,4,5,6,7,8,25), "etiology"] <- "HDV RNA-"
seu@meta.data[seu@meta.data$tissue == "Slide_3" & seu@meta.data$fov %in% c(17,18,19,20,21,22), "etiology"] <- "HDV RNA-"

#patient_type
# Initialize 'patient_type' column
seu@meta.data$patient_type <- NA
# B07
seu@meta.data[seu@meta.data$tissue == "Slide_2" & seu@meta.data$fov %in% c(13,14,15,16,17,18,19,20,21), "patient_type"] <- "B07"
# BH129
seu@meta.data[seu@meta.data$tissue == "Slide_3" & seu@meta.data$fov %in% c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,23,24,25), "patient_type"] <- "BH129"
# D03
seu@meta.data[seu@meta.data$tissue == "Slide_2" & seu@meta.data$fov %in% c(1,2,3,4,5), "patient_type"] <- "D03"
# D13
seu@meta.data[seu@meta.data$tissue == "Slide_2" & seu@meta.data$fov %in% c(6,7,8,9,10,11,12,22,23,24), "patient_type"] <- "D13"
# N10
seu@meta.data[seu@meta.data$tissue == "Slide_1" & seu@meta.data$fov %in% c(1,2,3,4,5,6,7,8,25), "patient_type"] <- "N10"
# N02
seu@meta.data[seu@meta.data$tissue == "Slide_3" & seu@meta.data$fov %in% c(17,18,19,20,21,22), "patient_type"] <- "N02"


# Initialize the new columns
seu@meta.data$Ag_etiology <- NA
seu@meta.data$Ag_general_etiology <- NA
seu@meta.data$Ag_general <- NA

### Ag_etiology
# Neg_HDVn
seu@meta.data[seu@meta.data$tissue == "Slide_1" & seu@meta.data$fov %in% c(1,2,3,4,5,6,7,8,25), "Ag_etiology"] <- "Neg_HDVn"
# Neg_HBV
seu@meta.data[seu@meta.data$tissue == "Slide_3" & seu@meta.data$fov %in% c(1,4,11,16,24,25), "Ag_etiology"] <- "Neg_HBV"
# S_low_HBV
seu@meta.data[seu@meta.data$tissue == "Slide_2" & seu@meta.data$fov %in% c(13,14,15,17,19,20), "Ag_etiology"] <- "S_low_HBV"
seu@meta.data[seu@meta.data$tissue == "Slide_3" & seu@meta.data$fov %in% c(13), "Ag_etiology"] <- "S_low_HBV"
# S_high_HBV
seu@meta.data[seu@meta.data$tissue == "Slide_2" & seu@meta.data$fov %in% c(16,18), "Ag_etiology"] <- "S_high_HBV"
seu@meta.data[seu@meta.data$tissue == "Slide_3" & seu@meta.data$fov %in% c(2,3,5,6,12,14,23), "Ag_etiology"] <- "S_high_HBV"
# S+/-_HBV
seu@meta.data[seu@meta.data$tissue == "Slide_3" & seu@meta.data$fov %in% c(7,8,9,10,15), "Ag_etiology"] <- "S+/-_HBV"
# S+/-_HDVn
seu@meta.data[seu@meta.data$tissue == "Slide_3" & seu@meta.data$fov %in% c(17,18,19,20,21,22), "Ag_etiology"] <- "S+/-_HDVn"
# S&D_low_HDVp
seu@meta.data[seu@meta.data$tissue == "Slide_2" & seu@meta.data$fov %in% c(1,5,6,7,9,11,12), "Ag_etiology"] <- "S&D_low_HDVp"
# S&D_HDVp
seu@meta.data[seu@meta.data$tissue == "Slide_2" & seu@meta.data$fov %in% c(2,3,4,8,10), "Ag_etiology"] <- "S&D_HDVp"
#HC
#seu@meta.data[seu@meta.data$tissue == "Slide_healthy","Ag_etiology" ] <- "HC"
seu@meta.data[seu@meta.data$tissue == "Slide_healthy","etiology" ] <- "HC"

### Ag_general_etiology
# Neg_HDVn
seu@meta.data[seu@meta.data$tissue == "Slide_1" & seu@meta.data$fov %in% c(1,2,3,4,5,6,7,8,25), "Ag_general_etiology"] <- "Neg_HDVn"
# Neg_HBV
seu@meta.data[seu@meta.data$tissue == "Slide_3" & seu@meta.data$fov %in% c(1,4,11,16,24,25), "Ag_general_etiology"] <- "Neg_HBV"
# S+_HBV
seu@meta.data[seu@meta.data$tissue == "Slide_2" & seu@meta.data$fov %in% c(13,14,15,16,17,18,19,20), "Ag_general_etiology"] <- "S+_HBV"
seu@meta.data[seu@meta.data$tissue == "Slide_3" & seu@meta.data$fov %in% c(2,3,5,6,12,13,14,23), "Ag_general_etiology"] <- "S+_HBV"
# S+/-_HBV
seu@meta.data[seu@meta.data$tissue == "Slide_3" & seu@meta.data$fov %in% c(7,8,9,10,15), "Ag_general_etiology"] <- "S+/-_HBV"
# S+/-_HDVn
seu@meta.data[seu@meta.data$tissue == "Slide_3" & seu@meta.data$fov %in% c(17,18,19,20,21,22), "Ag_general_etiology"] <- "S+/-_HDVn"
# S&D_HDVp
seu@meta.data[seu@meta.data$tissue == "Slide_2" & seu@meta.data$fov %in% c(1,2,3,4,5,6,7,8,9,10,11,12), "Ag_general_etiology"] <- "S&D_HDVp"

### Ag_general
# Neg
seu@meta.data[seu@meta.data$tissue == "Slide_1" & seu@meta.data$fov %in% c(1,2,3,4,5,6,7,8,25), "Ag_general"] <- "Neg"
seu@meta.data[seu@meta.data$tissue == "Slide_3" & seu@meta.data$fov %in% c(1,4,11,16,24,25), "Ag_general"] <- "Neg"
# S+
seu@meta.data[seu@meta.data$tissue == "Slide_2" & seu@meta.data$fov %in% c(13,14,15,16,17,18,19,20), "Ag_general"] <- "S+"
seu@meta.data[seu@meta.data$tissue == "Slide_3" & seu@meta.data$fov %in% c(2,3,5,6,12,13,14,23), "Ag_general"] <- "S+"
# S+/-
seu@meta.data[seu@meta.data$tissue == "Slide_3" & seu@meta.data$fov %in% c(7,8,9,10,15,17,18,19,20,21,22), "Ag_general"] <- "S+/-"
# S&D
seu@meta.data[seu@meta.data$tissue == "Slide_2" & seu@meta.data$fov %in% c(1,2,3,4,5,6,7,8,9,10,11,12), "Ag_general"] <- "S&D"


#Cell types
seu@meta.data$subset <- subset_type[seu@meta.data$subset]
seu@meta.data$refined <- cell_types[seu@meta.data$refined]
seu@meta.data$new_anot <- new_anot[seu@meta.data$refined]
seu@meta.data$grouped_anot <- grouped_anot[seu@meta.data$subset]


#Color palettes
color_palette <- c(
  paletteer_d("ggsci::default_igv"),
  paletteer_d("ggsci::category20_d3"),
  paletteer_d("ggsci::default_ucscgb")
)

meta <- seu@meta.data
#Tissue pal
tissue_pal <- color_palette
names(tissue_pal) <- unique(meta$tissue)
#Etiology pal
etiology_pal <- color_palette
names(etiology_pal) <- unique(meta$etiology)
patype_col <- color_palette
names(patype_col) <- unique(meta$patient_type)
sub_col <- c(
  
  Hepatocytes     = "#924822",  # mery_brown
  
  Myeloid_cells   = "#6BD76B",  # mery_green
  
  B_cell_lineage  = "#73c0eD",  # mery_blue
  
  T_cell_lineage  = "#EEDE5B",  # mery_yellow
  
  Non_parenchymal     = "#f5ccda",   # mery_pink
  
  Im_hepatocytes      =  "coral"
  
)

#Wide annotation
wide <- c(
  
  Hepatocytes     = "#ED9824",  # mery_brown
  
  Myeloid_cells   = "#5387C9",  # mery_green
  
  B_cell_lineage  = "#772D8B",  # mery_blue
  
  T_cell_lineage  = "#A8C686",  # mery_yellow
  
  Non_parenchymal     = "#F05365"
  
)

meta$wide <- gsub(pattern = "Im_hepatocytes" ,replacement = "Hepatocytes" ,meta$subset)

refined_col <- c(
  Gamma_delta_T_cells       = "#304796",
  NK_cells                  = "#C2F7F5",
  NKT_cells                 = "#8DD1EA",
  Regulatory_T_cells        = "#3571FA",
  Tem_Trm_cytotoxic_T_cells = "#4ae9ff",
  Effector_helper_T_cells   = "#0091AB",
  Naive_T_cells             = "#BBD6DB",
  Memory_B_cells            = "#E34183",
  Plasma_cells              = "#F1B8EF",
  Naive_B_cells             = "#FE64F9",
  ## MYELOIDS
  Monocytes                 = "#a11191",
  KC1                       = "#ff0077",
  M2_LYVE1                  = "#ffbf00",
  KC2                       = "#56c9f2",
  DCs_CD1C                  = "#02fa8d",
  M1                        = "#9c78fe",
  ###
  Endothelial_cells_2       = "#FFDC5F",
  Fibroblasts               = "#DB9925",
  Endothelial_cells_1       = "#F9F452",
  Smooth_muscle_cells       = "#FF8D08",
  Endothelial_cells_4       = "#CCC618",
  Endothelial_cells_3       = "#EADE8D",
  Hepatocyte_2              = "#c315f9",
  Hepatocyte_1              = "#E23F36",
  Hepatocyte_6              = "#3571FA",
  Hepatocyte_3              = "#DB9925",
  Hepatocyte_5              = "#20aa87",
  Hepatocyte_4              = "#CCC618",
  Cholangiocytes            = "#E23F36",
  Hepatocytes               = "#96307B",
  Cycling_B_lineage_cells   = "coral",
  Hepatocyte                = "#96307B"
)


##Figure 4B
logy <- meta$etiology
names(logy) <- meta$cell_names
all_int$etiology <- logy[all_int$id_source]
unique(all_int$etiology)

hepatocytes <- c("Hepatocyte 1", "Hepatocyte 2", "Hepatocyte 3", "Hepatocyte 5", "Hepatocyte 6")

cut_int <- all_int[all_int$etiology %in% c("HDV RNA+","HBV","HC") & all_int$likelihood > 0.5,]

#n_cells <- table(meta[meta$etiology %in% c("HDV RNA+","HBV","HC") & meta$new_anot %in% c("Kupffer_cells","Hepatocytes"),]$etiology)

# hep_kc1 <- cut_int[cut_int$refined_receptor %in% c("KC2", "KC1", hepatocytes) & cut_int$refined_source %in% c(hepatocytes,"KC2", "KC1") ,]
# table(hep_kc1$etiology)

### Hepatocytes _ all
int_table_all <- as.data.frame(table(cut_int$etiology))
hep_all <- cut_int[cut_int$refined_receptor %in% c( hepatocytes) | cut_int$refined_source %in% c(hepatocytes) ,]
int_hep_all <- as.data.frame(table(hep_all$etiology))

final_table <- int_table_all
names(final_table) <- c("etiology", "all_int")
final_table$hep_all <- int_hep_all$Freq
final_table$perc <- (final_table$hep_all / final_table$all_int )*100
final_table$norm <- final_table$perc / final_table[2,4]
df_hep <- final_table  # Save it
df_hep$Condition <- "Hepatocytes All"

#Hep with KC2
hep_all <- cut_int[cut_int$refined_receptor %in% c( hepatocytes, "KC2") & cut_int$refined_source %in% c(hepatocytes,"KC2") ,]
int_hep_all <- as.data.frame(table(hep_all$etiology))

final_table <- int_table_all
names(final_table) <- c("etiology", "all_int")
final_table$hep_all <- int_hep_all$Freq
final_table$perc <- (final_table$hep_all / final_table$all_int )*100
final_table$norm <- final_table$perc / final_table[2,4]
df_kc2 <- final_table # Save it
df_kc2$Condition <- "Hep + KC2"


# --- Combine them ---
all_data <- rbind(df_hep, df_kc2)

# ---------------------------------------------------------
# 2. FIX ORDERING (Very Important)
# ---------------------------------------------------------

library(ggplot2)

# 1. Ensure Etiology is ordered
all_data$etiology <- factor(all_data$etiology, 
                            levels = c("HC", "HBV", "HDV RNA+"))

# 2. Ensure Tables (Condition) are ordered
all_data$Condition <- factor(all_data$Condition, 
                             levels = c("Hepatocytes All", "Hep + KC1", "Hep + KC2"))

# 3. Create the Plot
p <- ggplot(all_data, aes(x = Condition, y = norm, fill = etiology)) +
  # Create bars
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black", width = 0.7) +
  
  # Custom Colors
  scale_fill_manual(values = c("HC" = "#fe4a49", 
                               "HBV" = "#2ab7ca", 
                               "HDV RNA+" = "#fed766")) + 
  
  # Theme settings for Publication
  theme_minimal() +
  theme(
    # 1. Remove all grid lines
    panel.grid = element_blank(),
    
    # 2. Strong X and Y Axis Lines (Black and Thicker)
    axis.line = element_line(color = "black", linewidth = 1),
    
    # 3. Add Axis Ticks (Essential for "perfect" visibility)
    axis.ticks = element_line(color = "black", linewidth = 1),
    axis.ticks.length = unit(0.2, "cm"),
    
    # 4. Text formatting (Ensure it is solid black, not gray)
    axis.text.x = element_text(size = 12, face = "bold", color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 12, face = "bold", color = "black"),
    
    # Hide "Condition" label
    axis.title.x = element_blank(), 
    
    # Legend formatting
    legend.position = "top"
  ) +
  
  # Labels
  labs(title = NULL,
       y = "Normalized Frequency (Ref: HC)",
       fill = "Etiology")

png(filename = "~/SPATIAL/Maria_CosMx/Spatial_liver/figures/outs/figure4_B.png", width = 8, height = 6, units = "in", res = 400)
print(p)
dev.off

##Cell enrichment
#Figure 4C
library(readr)
all <- read_csv("Maria_v2/Post_analysis/Celltype_enrichment/enrichment_files/all.csv")

all_kc1 <- all[all$from == "Hepatocyte" & all$to %in% c("KC1","KC2") & all$etiology != "HDV RNA-",]
#all_rich$to <- gsub(pattern = "MT-macs",replacement = "SPP1-macs",x = all_rich$to)
column_from = "from"
x = "bin"
y = "enrichment"
group = "to"
color = "to"
palette_color <- refined_col
method = "loess"
linetype = "dashed"
linewidth = 1
all_kc1$interval_numeric <- factor(all_kc1$bin, 
                                    levels = sort(unique(all_kc1$bin)))

all_kc1$etiology <- factor(all_kc1$etiology, 
                           levels = c("HC", "HBV", "HDV RNA+"))
png(filename = "~/SPATIAL/Maria_CosMx/Spatial_liver/figures/outs/coocurrence_hep_kc1.png",
    width = 14, height = 6, units = "in", res = 1200)
# 1. Ensure Etiology is ordered
library(ggplot2)

# 1. Ensure Etiology is ordered
all_kc1$etiology <- factor(all_kc1$etiology, 
                           levels = c("HC", "HBV", "HDV RNA+"))

# 2. Plot
p <- ggplot(all_kc1, aes(
  x = .data[[x]],
  y = .data[[y]],
  group = .data[[group]],
  color = .data[[color]]
)) +
  # 1. Data Lines: Reduced from 3 -> 1.5
  geom_smooth(alpha = 0.05, linewidth = 1.5) + 
  
  # 2. Horizontal Line: Reduced from 2 -> 1
  geom_hline(
    yintercept = 0,
    color = "black",
    linewidth = 1, 
    linetype = "dashed"
  ) +
  labs(x = NULL, y = NULL) + 
  facet_wrap(~ etiology) + 
  theme_linedraw() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    
    # 3. Square Border: Reduced from 3 -> 1.5
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
    
    # Hide axis lines (border handles it)
    axis.line = element_blank(),
    
    legend.position = "none",
    
    # No text labels
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    
    # 4. Axis Ticks: Reduced from 3 -> 1.5
    axis.ticks.x = element_line(linewidth = 1.5, color = "black"),
    axis.ticks.y = element_line(linewidth = 1.5, color = "black"),
    axis.ticks.length = unit(6, "pt"), # Reduced length slightly to match
    
    # Text size
    text = element_text(size = 20, face = "bold")
  ) +
  scale_color_manual(values = palette_color) +
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 10),   
    labels = NULL 
  )

p

dev.off()




## PErform same aalysis but only EXPLAN s+ and s-
cut_int <- all_int[all_int$fov %in% c(1,4,11,16,24,25,13,2,3,5,6,12,14,23) & all_int$tissue == "Slide_3"  & all_int$likelihood > 0.5,]
explant_anot <- seu$Ag_general
names(explant_anot) <- seu$cell_names
cut_int$ag_general <- explant_anot[cut_int$id_source]


### Hepatocytes _ all
int_table_all <- as.data.frame(table(cut_int$ag_general))
hep_all <- cut_int[cut_int$refined_receptor %in% c( hepatocytes) | cut_int$refined_source %in% c(hepatocytes) ,]
int_hep_all <- as.data.frame(table(hep_all$ag_general))

final_table <- int_table_all
names(final_table) <- c("ag_general", "all_int")
final_table$hep_all <- int_hep_all$Freq
final_table$perc <- (final_table$hep_all / final_table$all_int )*100
final_table$norm <- final_table$perc / final_table[1,4]
df_hep <- final_table  # Save it
df_hep$Condition <- "Hepatocytes All"

#Hep with KC1
hep_all <- cut_int[cut_int$refined_receptor %in% c( hepatocytes, "KC1") & cut_int$refined_source %in% c(hepatocytes,"KC1") ,]
int_hep_all <- as.data.frame(table(hep_all$ag_general))

final_table <- int_table_all
names(final_table) <- c("ag_general", "all_int")
final_table$hep_all <- int_hep_all$Freq
final_table$perc <- (final_table$hep_all / final_table$all_int )*100
final_table$norm <- final_table$perc / final_table[1,4]
df_kc1 <- final_table # Save it
df_kc1$Condition <- "Hep + KC1"


#Hep with KC2
hep_all <- cut_int[cut_int$refined_receptor %in% c( hepatocytes, "KC2") & cut_int$refined_source %in% c(hepatocytes,"KC2") ,]
int_hep_all <- as.data.frame(table(hep_all$ag_general))

final_table <- int_table_all
names(final_table) <- c("ag_general", "all_int")
final_table$hep_all <- int_hep_all$Freq
final_table$perc <- (final_table$hep_all / final_table$all_int )*100
final_table$norm <- final_table$perc / final_table[1,4]
df_kc2 <- final_table # Save it
df_kc2$Condition <- "Hep + KC2"

# --- Combine them ---
all_data <- rbind(df_hep, df_kc1, df_kc2)

# ---------------------------------------------------------
# 2. FIX ORDERING (Very Important)
# ---------------------------------------------------------

library(ggplot2)

# 1. Ensure Etiology is ordered
all_data$ag_general  <- factor(all_data$ag_general, 
                            levels = c("Neg", "S+"))

# 2. Ensure Tables (Condition) are ordered
all_data$Condition <- factor(all_data$Condition, 
                             levels = c("Hepatocytes All", "Hep + KC1", "Hep + KC2"))

# 3. Create the Plot
p <- ggplot(all_data, aes(x = Condition, y = norm, fill = ag_general)) +
  # Create bars
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black", width = 0.7) +
  
  # Custom Colors
  scale_fill_manual(values = c("Neg" = "#02876f", 
                               "S+" = "#f07801")) + 
  
  # Theme settings for Publication
  theme_minimal() +
  theme(
    # 1. Remove all grid lines
    panel.grid = element_blank(),
    
    # 2. Strong X and Y Axis Lines (Black and Thicker)
    axis.line = element_line(color = "black", linewidth = 1),
    
    # 3. Add Axis Ticks (Essential for "perfect" visibility)
    axis.ticks = element_line(color = "black", linewidth = 1),
    axis.ticks.length = unit(0.2, "cm"),
    
    # 4. Text formatting (Ensure it is solid black, not gray)
    axis.text.x = element_text(size = 12, face = "bold", color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 12, face = "bold", color = "black"),
    
    # Hide "Condition" label
    axis.title.x = element_blank(), 
    
    # Legend formatting
    legend.position = "top"
  ) +
  
  # Labels
  labs(title = NULL,
       y = "Normalized Frequency (Ref: HC)",
       fill = "Etiology")

png(filename = "~/SPATIAL/Maria_CosMx/Spatial_liver/figures/outs/figure4_E.png", width = 10, height = 6, units = "in", res = 400)
print(p)
dev.off()




#### Figure 4 E
all_kc1 <- all[all$from == "Hepatocyte" & all$to %in% c("KC1","KC2") & all$fov %in% c(1,4,11,16,24,25,13,2,3,5,6,12,14,23) & all$tissue == "Slide_3" ,]
#all_rich$to <- gsub(pattern = "MT-macs",replacement = "SPP1-macs",x = all_rich$to)
column_from = "from"
x = "bin"
y = "enrichment"
group = "to"
color = "to"
palette_color <- refined_col
method = "loess"
linetype = "dashed"
linewidth = 1
all_kc1$interval_numeric <- factor(all_kc1$bin, 
                                   levels = sort(unique(all_kc1$bin)))

all_kc1$Ag_general <- factor(all_kc1$Ag_general, 
                           levels = c("Neg", "S+"))
png(filename = "~/SPATIAL/Maria_CosMx/Spatial_liver/figures/outs/coocurrence_hep_kc_explant.png",
    width = 9.2, height = 6, units = "in", res = 1200)
# 1. Ensure Etiology is ordered
library(ggplot2)

# 1. Ensure Etiology is ordered
all_kc1$Ag_general <- factor(all_kc1$Ag_general, 
                           levels = c("Neg", "S+"))

# 2. Plot
p <- ggplot(all_kc1, aes(
  x = .data[[x]],
  y = .data[[y]],
  group = .data[[group]],
  color = .data[[color]]
)) +
  # 1. Data Lines: Reduced from 3 -> 1.5
  geom_smooth(alpha = 0.05, linewidth = 1.5) + 
  
  # 2. Horizontal Line: Reduced from 2 -> 1
  geom_hline(
    yintercept = 0,
    color = "black",
    linewidth = 1, 
    linetype = "dashed"
  ) +
  labs(x = NULL, y = NULL) + 
  facet_wrap(~ Ag_general) + 
  theme_linedraw() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    
    # 3. Square Border: Reduced from 3 -> 1.5
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
    
    # Hide axis lines (border handles it)
    axis.line = element_blank(),
    
    legend.position = "none",
    
    # No text labels
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    
    # 4. Axis Ticks: Reduced from 3 -> 1.5
    axis.ticks.x = element_line(linewidth = 1.5, color = "black"),
    axis.ticks.y = element_line(linewidth = 1.5, color = "black"),
    axis.ticks.length = unit(6, "pt"), # Reduced length slightly to match
    
    # Text size
    text = element_text(size = 20, face = "bold")
  ) +
  scale_color_manual(values = palette_color) +
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 10),   
    labels = NULL 
  )

p

dev.off()

