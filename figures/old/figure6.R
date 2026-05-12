library(Seurat)
library(paletteer)
library(plyr)
library(dplyr)
library(Seurat)
#library(clusterProfiler)
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


seu <-  readRDS("/home/mmoro/SPATIAL/Maria_CosMx/Maria_v2/Objects/seurats_all.RDS")
seu <- JoinLayers(seu)

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
volcano <- function(anot = "subset", ct, dif_col = "tissue", seu,id1,id2){
  metadata <- seu@meta.data
  # Filter for myeloid cells
  if(ct != "all"){
    myeloid_metadata <- metadata[metadata[[anot]] == ct, ]
  }else{
    myeloid_metadata <- metadata
  }
  # Use table to count cells in each tissue
  myeloid_counts_table <- table(myeloid_metadata[[dif_col]])
  #seu
  if(ct != "all"){
    object <- seu[,rownames(seu@meta.data[seu@meta.data[[anot]]== ct,])]
  }else{
    object <- seu
  }
  #object <- object[list_of_genes,]
  object <- NormalizeData(object)
  object <- ScaleData(object)
  #object@active.ident <- as.factor(object@meta.data[[dif_col]])
  
  # Convert the specified column to a factor
  object@meta.data[[dif_col]] <- as.factor(object@meta.data[[dif_col]])
  # Update the active identity with the factor column
  object <- SetIdent(object, value = object@meta.data[[dif_col]])
  
  #IBD vs NHC
  deg_results <- FindMarkers(object, ident.1 = id1, ident.2 = id2)
  deg_results <- na.omit(deg_results)
  deg_results$genes <- rownames(deg_results)
  #Diff expressed
  deg_results$diffexpressed <- "NO"
  deg_results$diffexpressed[deg_results$avg_log2FC > log2(1.2) & deg_results$p_val < 0.05] <- "UP"
  deg_results$diffexpressed[deg_results$avg_log2FC < -log2(1.2) & deg_results$p_val < 0.05] <- "DOWN"
  deg_results$diffexpressed[deg_results$avg_log2FC > log2(1.2) & deg_results$p_val_adj < 0.05] <- "UPP"
  deg_results$diffexpressed[deg_results$avg_log2FC < -log2(1.2) & deg_results$p_val_adj < 0.05] <- "DOWNN"
  #Label of most down-up regulated genes
  deg_results$delabel <- NA
  deg_results$delabel[deg_results$diffexpressed != "NO"] <- deg_results$genes[deg_results$diffexpressed != "NO"]
  deg_results$p_val <- ifelse(deg_results$p_val < 1e-300, 1e-300, deg_results$p_val)
  diff_colors <- c("UPP" = "#803800", "DOWNN" = "#003F54", "UP" = "#B47846", "DOWN" = "steelblue")
  
  genes_to_label <- deg_results[deg_results$diffexpressed %in% c("UPP","DOWNN"),]$delabel
  
  deg_results$delabel <- ifelse(deg_results$gene %in% genes_to_label,
                                deg_results$gene, NA)
  
  #Plot
  p <- ggplot(data = deg_results, aes(x = avg_log2FC, y = -log10(p_val), col = diffexpressed)) +
    geom_point(size = 2) +
    scale_color_manual(values = diff_colors) +
    theme_classic() +
    theme(text = element_text(family = "Helvetica")) +
    guides(color = guide_legend(override.aes = list(shape = 1))) +
    theme(legend.position = "none") +
    geom_label_repel(
      aes(label = delabel),      # ONLY shows labels for the selected genes
      size = 18 / .pt,
      segment.color = "black",
      segment.size = 2,
      label.padding = unit(0.4, "lines"),
      label.size = 1,
      fontface = "bold"
    ) +
    theme(
      plot.title = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.line = element_line(linewidth = 2),
      axis.ticks.length = unit(0.5, "cm"),
      axis.ticks = element_line(linewidth = 2),  # Increased the width of the axis ticks
      axis.text = element_text(size = 22)
    )
  
  # Print the plot
  return(list(p,deg_results))
  
}

#Polygons Slide 3 FOV 9
Slide2_pols <- read_csv("/home/mmoro/SPATIAL/Maria_CosMx/Maria_v2/Polygons/Slide_3.csv")
Slide2_pols <- Slide2_pols[Slide2_pols$fov == "9",]
meta_cut <- meta[meta$tissue == "Slide_3" & meta$fov == 9,]
Slide2_pols <- Slide2_pols[Slide2_pols$cell_names %in% meta_cut$cell_names,]
Slide2_pols[["wide"]] <- plyr::mapvalues(x = Slide2_pols[["cell_names"]],
                                         from = meta_cut[["cell_names"]],
                                         to = meta_cut[["wide"]])
Slide2_pols[["wide"]] <- as.factor(Slide2_pols[["wide"]])


png(filename = "~/SPATIAL/Maria_CosMx/Spatial_liver/figures/outs/fig6_wide_S3F9.png",width = 8,height = 8,units = "in",res = 1200)
p <-
  ggplot2::ggplot(Slide2_pols, ggplot2::aes(x = .data[["x_global_px"]], y = .data[["y_global_px"]])) +
  ggplot2::geom_polygon(ggplot2::aes(group = .data[["cell_names"]], fill = .data[["wide"]]),
                        color = "#000000") +
  ggdark::dark_theme_gray(base_family = "Fira Sans Condensed Light", base_size = 20) +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    axis.title = ggplot2::element_blank(),
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank()
  )  + ggplot2::scale_fill_manual(values = wide) +guides(fill="none")
p
dev.off()

#Polygons Slide 2 FOV 7 
Slide2_pols <- read_csv("/home/mmoro/SPATIAL/Maria_CosMx/Maria_v2/Polygons/Slide_2.csv")
Slide2_pols <- Slide2_pols[Slide2_pols$fov == "7",]
meta_cut <- meta[meta$tissue == "Slide_2" & meta$fov == 7,]
Slide2_pols <- Slide2_pols[Slide2_pols$cell_names %in% meta_cut$cell_names,]
Slide2_pols[["wide"]] <- plyr::mapvalues(x = Slide2_pols[["cell_names"]],
                                         from = meta_cut[["cell_names"]],
                                         to = meta_cut[["wide"]])
Slide2_pols[["wide"]] <- as.factor(Slide2_pols[["wide"]])


png(filename = "~/SPATIAL/Maria_CosMx/Spatial_liver/figures/outs/fig6_wide_S2F7.png",width = 8,height = 8,units = "in",res = 1200)
p <-
  ggplot2::ggplot(Slide2_pols, ggplot2::aes(x = .data[["x_global_px"]], y = .data[["y_global_px"]])) +
  ggplot2::geom_polygon(ggplot2::aes(group = .data[["cell_names"]], fill = .data[["wide"]]),
                        color = "#000000") +
  ggdark::dark_theme_gray(base_family = "Fira Sans Condensed Light", base_size = 20) +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    axis.title = ggplot2::element_blank(),
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank()
  )  + ggplot2::scale_fill_manual(values = wide) +guides(fill="none")
p
dev.off()


#Polygons Slide 1 FOV 25 
Slide2_pols <- read_csv("/home/mmoro/SPATIAL/Maria_CosMx/Maria_v2/Polygons/Slide_1.csv")
Slide2_pols <- Slide2_pols[Slide2_pols$fov == "25",]
meta_cut <- meta[meta$tissue == "Slide_1" & meta$fov == 25,]
Slide2_pols <- Slide2_pols[Slide2_pols$cell_names %in% meta_cut$cell_names,]
Slide2_pols[["wide"]] <- plyr::mapvalues(x = Slide2_pols[["cell_names"]],
                                         from = meta_cut[["cell_names"]],
                                         to = meta_cut[["wide"]])
Slide2_pols[["wide"]] <- as.factor(Slide2_pols[["wide"]])


png(filename = "~/SPATIAL/Maria_CosMx/Spatial_liver/figures/outs/fig6_wide_S1F25.png",width = 8,height = 8,units = "in",res = 1200)
p <-
  ggplot2::ggplot(Slide2_pols, ggplot2::aes(x = .data[["x_global_px"]], y = .data[["y_global_px"]])) +
  ggplot2::geom_polygon(ggplot2::aes(group = .data[["cell_names"]], fill = .data[["wide"]]),
                        color = "#000000") +
  ggdark::dark_theme_gray(base_family = "Fira Sans Condensed Light", base_size = 20) +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    axis.title = ggplot2::element_blank(),
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank()
  )  + ggplot2::scale_fill_manual(values = wide) +guides(fill="none")
p
dev.off()

#Polygons Slide 4 FOV 4 
#Slide2_pols <- read_csv("/home/mmoro/SPATIAL/Maria_CosMx/Maria_v2/Polygons/Slide_4.csv")
#Slide2_pols <- Slide2_pols[Slide2_pols$fov == "25",]
meta_cut <- meta[meta$tissue == "Slide_healthy" & meta$fov ==4,]

png(filename = "~/SPATIAL/Maria_CosMx/Spatial_liver/figures/outs/fig6_wide_SHCF4.png", 
    width = 8, height = 8, units = "in", res = 1200)

p <- ggplot2::ggplot(meta_cut, ggplot2::aes(x = .data[["CenterX_global_px"]], y = .data[["CenterY_global_px"]], color = .data[["wide"]])) +
  ggplot2::geom_point(size = 3) +
  ggdark::dark_theme_gray(base_family = "Fira Sans Condensed Light", base_size = 20) +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    axis.title = ggplot2::element_blank(),
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank()
  ) + 
  ggplot2::scale_color_manual(values = wide) +
  ggplot2::guides(color = "none")

print(p)
dev.off()


## Only showing hepatocytes
## refineD hepatos
#Polygons Slide 3 FOV 9, Showing refined
Slide2_pols <- read_csv("/home/mmoro/SPATIAL/Maria_CosMx/Maria_v2/Polygons/Slide_3.csv")
Slide2_pols <- Slide2_pols[Slide2_pols$fov == "9",]
meta_cut <- meta[meta$tissue == "Slide_3" & meta$fov == 9,]
Slide2_pols <- Slide2_pols[Slide2_pols$cell_names %in% meta_cut$cell_names,]
Slide2_pols[["refined"]] <- plyr::mapvalues(x = Slide2_pols[["cell_names"]],
                                            from = meta_cut[["cell_names"]],
                                            to = meta_cut[["refined"]])
Slide2_pols[["refined"]] <- as.factor(Slide2_pols[["refined"]])

Slide2_pols <- Slide2_pols %>%
  mutate(refined = if_else(refined %in% c("Hepatocyte_1" , "Hepatocyte_2" ,"Hepatocyte_3" ,"Hepatocyte_5" ,"Hepatocyte_6"), refined, "Other"))
refined_col["Other"] <- "#393939"

png(filename = "~/SPATIAL/Maria_CosMx/Spatial_liver/figures/outs/fig6S3F9_refined.png",width = 8,height = 8,units = "in",res = 1200)
p <-
  ggplot2::ggplot(Slide2_pols, ggplot2::aes(x = .data[["x_global_px"]], y = .data[["y_global_px"]])) +
  ggplot2::geom_polygon(ggplot2::aes(group = .data[["cell_names"]], fill = .data[["refined"]]),
                        color = "#000000") +
  ggdark::dark_theme_gray(base_family = "Fira Sans Condensed Light", base_size = 20) +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    axis.title = ggplot2::element_blank(),
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank()
  )  + ggplot2::scale_fill_manual(values = refined_col) +guides(fill="none")
p
dev.off()


#Polygons Slide 2 FOV 7, Showing refined
Slide2_pols <- read_csv("/home/mmoro/SPATIAL/Maria_CosMx/Maria_v2/Polygons/Slide_2.csv")
Slide2_pols <- Slide2_pols[Slide2_pols$fov == "7",]
meta_cut <- meta[meta$tissue == "Slide_2" & meta$fov == 7,]
Slide2_pols <- Slide2_pols[Slide2_pols$cell_names %in% meta_cut$cell_names,]
Slide2_pols[["refined"]] <- plyr::mapvalues(x = Slide2_pols[["cell_names"]],
                                            from = meta_cut[["cell_names"]],
                                            to = meta_cut[["refined"]])
Slide2_pols[["refined"]] <- as.factor(Slide2_pols[["refined"]])

Slide2_pols <- Slide2_pols %>%
  mutate(refined = if_else(refined %in% c("Hepatocyte_1" , "Hepatocyte_2" ,"Hepatocyte_3" ,"Hepatocyte_5" ,"Hepatocyte_6"), refined, "Other"))
refined_col["Other"] <- "#393939"

png(filename = "~/SPATIAL/Maria_CosMx/Spatial_liver/figures/outs/fig6S2F7_refined.png",width = 8,height = 8,units = "in",res = 1200)
p <-
  ggplot2::ggplot(Slide2_pols, ggplot2::aes(x = .data[["x_global_px"]], y = .data[["y_global_px"]])) +
  ggplot2::geom_polygon(ggplot2::aes(group = .data[["cell_names"]], fill = .data[["refined"]]),
                        color = "#000000") +
  ggdark::dark_theme_gray(base_family = "Fira Sans Condensed Light", base_size = 20) +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    axis.title = ggplot2::element_blank(),
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank()
  )  + ggplot2::scale_fill_manual(values = refined_col) +guides(fill="none")
p
dev.off()


#Polygons Slide 1 FOV 25, Showing refined
Slide2_pols <- read_csv("/home/mmoro/SPATIAL/Maria_CosMx/Maria_v2/Polygons/Slide_1.csv")
Slide2_pols <- Slide2_pols[Slide2_pols$fov == "25",]
meta_cut <- meta[meta$tissue == "Slide_1" & meta$fov == 25,]
Slide2_pols <- Slide2_pols[Slide2_pols$cell_names %in% meta_cut$cell_names,]
Slide2_pols[["refined"]] <- plyr::mapvalues(x = Slide2_pols[["cell_names"]],
                                            from = meta_cut[["cell_names"]],
                                            to = meta_cut[["refined"]])
Slide2_pols[["refined"]] <- as.factor(Slide2_pols[["refined"]])

Slide2_pols <- Slide2_pols %>%
  mutate(refined = if_else(refined %in% c("Hepatocyte_1" , "Hepatocyte_2" ,"Hepatocyte_3" ,"Hepatocyte_5" ,"Hepatocyte_6"), refined, "Other"))
refined_col["Other"] <- "#393939"

png(filename = "~/SPATIAL/Maria_CosMx/Spatial_liver/figures/outs/fig6S1F25_refined.png",width = 8,height = 8,units = "in",res = 1200)
p <-
  ggplot2::ggplot(Slide2_pols, ggplot2::aes(x = .data[["x_global_px"]], y = .data[["y_global_px"]])) +
  ggplot2::geom_polygon(ggplot2::aes(group = .data[["cell_names"]], fill = .data[["refined"]]),
                        color = "#000000") +
  ggdark::dark_theme_gray(base_family = "Fira Sans Condensed Light", base_size = 20) +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    axis.title = ggplot2::element_blank(),
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank()
  )  + ggplot2::scale_fill_manual(values = refined_col) +guides(fill="none")
p
dev.off()



##
meta_cut <- meta[meta$tissue == "Slide_healthy" & meta$fov ==4,]
meta_cut <- meta_cut %>%
  mutate(refined = if_else(refined %in% c("Hepatocyte_1" , "Hepatocyte_2" ,"Hepatocyte_3" ,"Hepatocyte_5" ,"Hepatocyte_6"), refined, "Other"))
refined_col["Other"] <- "#393939"


png(filename = "~/SPATIAL/Maria_CosMx/Spatial_liver/figures/outs/fig6_refined_SHCF4.png", 
    width = 8, height = 8, units = "in", res = 1200)

p <- ggplot2::ggplot(meta_cut, ggplot2::aes(x = .data[["CenterX_global_px"]], y = .data[["CenterY_global_px"]], color = .data[["refined"]])) +
  ggplot2::geom_point(size = 3) +
  ggdark::dark_theme_gray(base_family = "Fira Sans Condensed Light", base_size = 20) +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    axis.title = ggplot2::element_blank(),
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank()
  ) + 
  ggplot2::scale_color_manual(values = refined_col) +
  ggplot2::guides(color = "none")

print(p)
dev.off()

## SEU NEIGH
umap_csv    <- "/home/mmoro/SPATIAL/Maria_CosMx/Maria_v2/Post_analysis/Nieghborhood/Run/30/umap.csv"
seu_neigh   <- "/home/mmoro/SPATIAL/Maria_CosMx/Maria_v2/Post_analysis/Nieghborhood/Run/30/seu_neigh.csv"

umap        <- read_csv(umap_csv)
seu_neigh   <- read_csv(seu_neigh)

# Organized 12-color palette (sorted from warm to cool hues)
neigh_cols <- c(
  "#CE3D32FF", # 1. Red
  "#F24B91", # 2. Dusty Pink
  "#BA6338FF", # 3. Burnt Orange
  "#F0E685FF", # 4. Pale Yellow
  "#6BD76BFF", # 5. Bright Green
  "#749B58FF", # 6. Olive Green
  "#5DB1DDFF", # 7. Sky Blue
  "#5050FFFF", # 8. Royal Blue
  "#466983FF", # 9. Steel Blue
  "#802268FF", # 10. Plum Purple
  "#F88A4C", # 11. Muted Purple/Grey
  "#D5A5DA"  # 12. Dark Brown
)

names(neigh_cols) <- unique(seu_neigh$leiden_res0.2)



### Plots neigh
Slide2_pols <- read_csv("/home/mmoro/SPATIAL/Maria_CosMx/Maria_v2/Polygons/Slide_1.csv")
Slide2_pols <- Slide2_pols[Slide2_pols$fov == "25",]
meta_cut <- seu_neigh[seu_neigh$tissue == "Slide_1" & seu_neigh$fov == 25,]
Slide2_pols <- Slide2_pols[Slide2_pols$cell_names %in% meta_cut$cell_names,]
Slide2_pols[["neigh"]] <- plyr::mapvalues(x = Slide2_pols[["cell_names"]],
                                            from = meta_cut[["cell_names"]],
                                            to = meta_cut[["leiden_res0.2"]])
Slide2_pols[["neigh"]] <- as.factor(Slide2_pols[["neigh"]])

png(filename = "~/SPATIAL/Maria_CosMx/Spatial_liver/figures/outs/fig6S1F25_neigh.png",width = 8,height = 8,units = "in",res = 1200)
p <-
  ggplot2::ggplot(Slide2_pols, ggplot2::aes(x = .data[["x_global_px"]], y = .data[["y_global_px"]])) +
  ggplot2::geom_polygon(ggplot2::aes(group = .data[["cell_names"]], fill = .data[["neigh"]]),
                        color = "#000000") +
  ggdark::dark_theme_gray(base_family = "Fira Sans Condensed Light", base_size = 20) +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    axis.title = ggplot2::element_blank(),
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank()
  )  + ggplot2::scale_fill_manual(values = neigh_cols) +guides(fill="none")
p
dev.off()

##
Slide2_pols <- read_csv("/home/mmoro/SPATIAL/Maria_CosMx/Maria_v2/Polygons/Slide_3.csv")
Slide2_pols <- Slide2_pols[Slide2_pols$fov == "9",]
meta_cut <- seu_neigh[seu_neigh$tissue == "Slide_3" & seu_neigh$fov == 9,]
Slide2_pols <- Slide2_pols[Slide2_pols$cell_names %in% meta_cut$cell_names,]
Slide2_pols[["neigh"]] <- plyr::mapvalues(x = Slide2_pols[["cell_names"]],
                                          from = meta_cut[["cell_names"]],
                                          to = meta_cut[["leiden_res0.2"]])
Slide2_pols[["neigh"]] <- as.factor(Slide2_pols[["neigh"]])

png(filename = "~/SPATIAL/Maria_CosMx/Spatial_liver/figures/outs/fig6S3F9_neigh.png",width = 8,height = 8,units = "in",res = 1200)
p <-
  ggplot2::ggplot(Slide2_pols, ggplot2::aes(x = .data[["x_global_px"]], y = .data[["y_global_px"]])) +
  ggplot2::geom_polygon(ggplot2::aes(group = .data[["cell_names"]], fill = .data[["neigh"]]),
                        color = "#000000") +
  ggdark::dark_theme_gray(base_family = "Fira Sans Condensed Light", base_size = 20) +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    axis.title = ggplot2::element_blank(),
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank()
  )  + ggplot2::scale_fill_manual(values = neigh_cols) +guides(fill="none")
p
dev.off()

##
Slide2_pols <- read_csv("/home/mmoro/SPATIAL/Maria_CosMx/Maria_v2/Polygons/Slide_2.csv")
Slide2_pols <- Slide2_pols[Slide2_pols$fov == "7",]
meta_cut <- seu_neigh[seu_neigh$tissue == "Slide_2" & seu_neigh$fov == 7,]
Slide2_pols <- Slide2_pols[Slide2_pols$cell_names %in% meta_cut$cell_names,]
Slide2_pols[["neigh"]] <- plyr::mapvalues(x = Slide2_pols[["cell_names"]],
                                          from = meta_cut[["cell_names"]],
                                          to = meta_cut[["leiden_res0.2"]])
Slide2_pols[["neigh"]] <- as.factor(Slide2_pols[["neigh"]])

png(filename = "~/SPATIAL/Maria_CosMx/Spatial_liver/figures/outs/fig6S2F7_neigh.png",width = 8,height = 8,units = "in",res = 1200)
p <-
  ggplot2::ggplot(Slide2_pols, ggplot2::aes(x = .data[["x_global_px"]], y = .data[["y_global_px"]])) +
  ggplot2::geom_polygon(ggplot2::aes(group = .data[["cell_names"]], fill = .data[["neigh"]]),
                        color = "#000000") +
  ggdark::dark_theme_gray(base_family = "Fira Sans Condensed Light", base_size = 20) +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    axis.title = ggplot2::element_blank(),
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank()
  )  + ggplot2::scale_fill_manual(values = neigh_cols) +guides(fill="none")
p
dev.off()

## Niegh health
meta_cut <- seu_neigh[seu_neigh$tissue == "Slide_healthy" & seu_neigh$fov ==4,]
meta_cut$leiden_res0.2 <- as.factor(meta_cut$leiden_res0.2)


png(filename = "~/SPATIAL/Maria_CosMx/Spatial_liver/figures/outs/fig6_neigh_SHCF4.png", 
    width = 8, height = 8, units = "in", res = 1200)

p <- ggplot2::ggplot(meta_cut, ggplot2::aes(x = .data[["CenterX_global_px"]], y = .data[["CenterY_global_px"]], color = .data[["leiden_res0.2"]])) +
  ggplot2::geom_point(size = 3) +
  ggdark::dark_theme_gray(base_family = "Fira Sans Condensed Light", base_size = 20) +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    axis.title = ggplot2::element_blank(),
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank()
  ) + 
  ggplot2::scale_color_manual(values = neigh_cols) +
  ggplot2::guides(color = "none")

print(p)
dev.off()


## Fig6B. Split umap by condition
seu_neigh_cut <- as.data.frame(seu_neigh)[,c("cell_names","leiden_res0.2")]
colnames(umap)[1] <- "cell_names"
umap <- umap %>% left_join(seu_neigh_cut, by = "cell_names" )
meta_cuti <- meta[,c("etiology","cell_names")]
umap <- umap %>% left_join(meta_cuti, by = "cell_names")



umap <- umap %>%
  mutate(etiology = factor(etiology, levels = c("HC", "HBV", "HDV RNA+", "HDV RNA-")))

# 1. Get a list of all unique etiologies in your dataset
etiologies <- unique(umap$etiology)

# 2. Get the global UMAP limits so all plots maintain the same scale/shape
x_limits <- range(umap$UMAP1, na.rm = TRUE)
y_limits <- range(umap$UMAP2, na.rm = TRUE)

# 3. Loop through each etiology, create its plot, and save it
for (eti in etiologies) {
  
  # Filter the data for just this specific etiology
  subset_umap <- subset(umap, etiology == eti)
  
  # Create a clean filename (replaces spaces or slashes with underscores just in case)
  #clean_eti_name <- gsub("[^A-Za-z0-9]", "_", eti)
  current_filename <- paste0("~/SPATIAL/Maria_CosMx/Spatial_liver/figures/outs/UMAP_etiology_", eti, "_minimal.png")
  
  # Open the PNG device
  png(filename = current_filename, width = 8, height = 8, units = "in", res = 1200)
  
  # Create the plot (removed facet_wrap and strip.text/panel.spacing since they aren't needed now)
  p <- ggplot(subset_umap, aes(x = UMAP1, y = UMAP2, color = as.factor(leiden_res0.2))) +
    geom_point(size = 0.5) + 
    scale_color_manual(values = neigh_cols) + 
    coord_cartesian(xlim = x_limits, ylim = y_limits) + # Keeps the UMAP shape consistent
    theme_void() + 
    theme(
      legend.position = "none"
    )
  
  # Print and close
  print(p)
  dev.off()
  
  # Optional: Print a progress message to the console
  cat("Saved:", current_filename, "\n")
}

#Stackplot cute pie
# healthy: 11-6
# active parenchyma: 2-0
# immune active: 3-4-7-8-10-5
# homeostatic/anti-inflamatory: 1-9
meta_cuti <- meta[,c("etiology","cell_names")]
seu_neigh <- seu_neigh %>% left_join(meta_cuti, by = "cell_names")

library(ggplot2)

# 1. Calculate cell counts instantly 
df <- as.data.frame(table(seu_neigh$leiden_res0.2, seu_neigh$etiology))
colnames(df) <- c("leiden_res0.2", "etiology", "value")

# 2. Set factor levels 
df$leiden_res0.2 <- as.factor(df$leiden_res0.2)
df$etiology <- factor(df$etiology, levels = c("HC", "HBV", "HDV RNA+", "HDV RNA-"))

# 3. Define the groups as a named list (this makes looping easy!)
cluster_groups <- list(
  "healthy"           = c(11, 6),
  "active_parenchyma" = c(2, 0),
  "immune_active"     = c(3, 4, 7, 8, 10, 5),
  "homeostatic"       = c(1, 9)
)

# 4. Define colors
etiology_colors <- c(
  "HC" = "#fe4a49", 
  "HBV" = "#2ab7ca", 
  "HDV RNA+" = "#fed766", 
  "HDV RNA-" = "#D866FE"
)

# 5. Loop through each group, filter the data, plot, and save!
for (group_name in names(cluster_groups)) {
  
  # Get the clusters for the current group
  clusters <- cluster_groups[[group_name]]
  
  # Filter and re-level the dataframe
  df_sub <- df[df$leiden_res0.2 %in% clusters, ]
  df_sub$leiden_res0.2 <- factor(df_sub$leiden_res0.2, levels = as.character(clusters))
  
  # Build the plot
  p <- ggplot(df_sub, aes(fill = etiology, y = leiden_res0.2, x = value)) + 
    geom_bar(position = "fill", stat = "identity", width = 1) +
    scale_fill_manual(values = etiology_colors) +
    theme_classic() + 
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      text = element_text(size = 25),
      axis.line = element_line(size = 1.5),
      axis.ticks = element_line(size = 1.2),
      axis.ticks.length = unit(0.3, "cm")
    )
  
  # Define the dynamic file path
  file_out <- paste0("~/SPATIAL/Maria_CosMx/Spatial_liver/figures/outs/stackplot_", group_name, "_fig6.png")
  
  # Save using your exact specifications
  png(filename = file_out, width = 9.2, height = 6, units = "in", res = 1200)
  print(p) # You must explicitly call print() when plotting inside a loop!
  dev.off()
}


#### Panel D
library(slingshot)
library(tradeSeq)
library(viridis)
sce_slingshot <- readRDS("~/SPATIAL/Maria_CosMx/Maria_v2/Post_analysis/Post_neigh/sce_slingshot.rds")
# 2. Create a data frame for plotting cells
traj <- 5
plot_df <- data.frame(
  UMAP_1 = reducedDim(sce_slingshot, 'UMAP')[, 1],
  UMAP_2 = reducedDim(sce_slingshot, 'UMAP')[, 2],
  pseudotime = slingPseudotime(sce_slingshot)[, traj]
)
curves <- slingCurves(sce_slingshot)
curve1_df <- as.data.frame(curves[[traj]]$s) # Get coordinates for the first curve
colnames(curve1_df) <- c("UMAP_1", "UMAP_2")

# 4. Generate the plot
umap_plot_with_line <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2)) +
  # Plot cells, colored by pseudotime
  geom_point(aes(color = pseudotime), size = 0.01, alpha = 0.8) +
  # Add the trajectory line ----------------------------------- NEW
  geom_path(data = curve1_df, aes(x = UMAP_1, y = UMAP_2), size = 1.2, color = "black") +
  # Use a nice color scale for continuous data
  scale_color_viridis(option = "magma", na.value = "lightgrey") +
  labs(
    title = "Slingshot Pseudotime on UMAP",
    subtitle = paste0("Trajectory ",traj," with Principal Curve", sep = ""),
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

# 5. Display and save the plot
png(filename = "~/SPATIAL/Maria_CosMx/Spatial_liver/figures/outs/UMAP1.png", width = 9.2, height = 6, units = "in", res = 800)
print(umap_plot_with_line)
dev.off()


traj <- 6
plot_df <- data.frame(
  UMAP_1 = reducedDim(sce_slingshot, 'UMAP')[, 1],
  UMAP_2 = reducedDim(sce_slingshot, 'UMAP')[, 2],
  pseudotime = slingPseudotime(sce_slingshot)[, traj]
)
curves <- slingCurves(sce_slingshot)
curve1_df <- as.data.frame(curves[[traj]]$s) # Get coordinates for the first curve
colnames(curve1_df) <- c("UMAP_1", "UMAP_2")

# 4. Generate the plot
umap_plot_with_line <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2)) +
  # Plot cells, colored by pseudotime
  geom_point(aes(color = pseudotime), size = 0.01, alpha = 0.8) +
  # Add the trajectory line ----------------------------------- NEW
  geom_path(data = curve1_df, aes(x = UMAP_1, y = UMAP_2), size = 1.2, color = "black") +
  # Use a nice color scale for continuous data
  scale_color_viridis(option = "magma", na.value = "lightgrey") +
  labs(
    title = "Slingshot Pseudotime on UMAP",
    subtitle = paste0("Trajectory ",traj," with Principal Curve", sep = ""),
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

# 5. Display and save the plot
png(filename = "~/SPATIAL/Maria_CosMx/Spatial_liver/figures/outs/UMAP2.png", width = 9.2, height = 6, units = "in", res = 800)
print(umap_plot_with_line)
dev.off()

#TRAJ 1
traj <- 1
plot_df <- data.frame(
  UMAP_1 = reducedDim(sce_slingshot, 'UMAP')[, 1],
  UMAP_2 = reducedDim(sce_slingshot, 'UMAP')[, 2],
  pseudotime = slingPseudotime(sce_slingshot)[, traj]
)
curves <- slingCurves(sce_slingshot)
curve1_df <- as.data.frame(curves[[traj]]$s) # Get coordinates for the first curve
colnames(curve1_df) <- c("UMAP_1", "UMAP_2")

# 4. Generate the plot
umap_plot_with_line <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2)) +
  # Plot cells, colored by pseudotime
  geom_point(aes(color = pseudotime), size = 0.01, alpha = 0.8) +
  # Add the trajectory line ----------------------------------- NEW
  geom_path(data = curve1_df, aes(x = UMAP_1, y = UMAP_2), size = 1.2, color = "black") +
  # Use a nice color scale for continuous data
  scale_color_viridis(option = "magma", na.value = "lightgrey") +
  labs(
    title = "Slingshot Pseudotime on UMAP",
    subtitle = paste0("Trajectory ",traj," with Principal Curve", sep = ""),
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

# 5. Display and save the plot
png(filename = "~/SPATIAL/Maria_CosMx/Spatial_liver/figures/outs/traj1.png", width = 9.2, height = 6, units = "in", res = 800)
print(umap_plot_with_line)
dev.off()

#TRAJ 2
traj <- 2
plot_df <- data.frame(
  UMAP_1 = reducedDim(sce_slingshot, 'UMAP')[, 1],
  UMAP_2 = reducedDim(sce_slingshot, 'UMAP')[, 2],
  pseudotime = slingPseudotime(sce_slingshot)[, traj]
)
curves <- slingCurves(sce_slingshot)
curve1_df <- as.data.frame(curves[[traj]]$s) # Get coordinates for the first curve
colnames(curve1_df) <- c("UMAP_1", "UMAP_2")

# 4. Generate the plot
umap_plot_with_line <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2)) +
  # Plot cells, colored by pseudotime
  geom_point(aes(color = pseudotime), size = 0.01, alpha = 0.8) +
  # Add the trajectory line ----------------------------------- NEW
  geom_path(data = curve1_df, aes(x = UMAP_1, y = UMAP_2), size = 1.2, color = "black") +
  # Use a nice color scale for continuous data
  scale_color_viridis(option = "magma", na.value = "lightgrey") +
  labs(
    title = "Slingshot Pseudotime on UMAP",
    subtitle = paste0("Trajectory ",traj," with Principal Curve", sep = ""),
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

# 5. Display and save the plot
png(filename = "~/SPATIAL/Maria_CosMx/Spatial_liver/figures/outs/traj2.png", width = 9.2, height = 6, units = "in", res = 800)
print(umap_plot_with_line)
dev.off()

#TRAJ 3
traj <- 3
plot_df <- data.frame(
  UMAP_1 = reducedDim(sce_slingshot, 'UMAP')[, 1],
  UMAP_2 = reducedDim(sce_slingshot, 'UMAP')[, 2],
  pseudotime = slingPseudotime(sce_slingshot)[, traj]
)
curves <- slingCurves(sce_slingshot)
curve1_df <- as.data.frame(curves[[traj]]$s) # Get coordinates for the first curve
colnames(curve1_df) <- c("UMAP_1", "UMAP_2")

# 4. Generate the plot
umap_plot_with_line <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2)) +
  # Plot cells, colored by pseudotime
  geom_point(aes(color = pseudotime), size = 0.01, alpha = 0.8) +
  # Add the trajectory line ----------------------------------- NEW
  geom_path(data = curve1_df, aes(x = UMAP_1, y = UMAP_2), size = 1.2, color = "black") +
  # Use a nice color scale for continuous data
  scale_color_viridis(option = "magma", na.value = "lightgrey") +
  labs(
    title = "Slingshot Pseudotime on UMAP",
    subtitle = paste0("Trajectory ",traj," with Principal Curve", sep = ""),
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

# 5. Display and save the plot
png(filename = "~/SPATIAL/Maria_CosMx/Spatial_liver/figures/outs/traj3.png", width = 9.2, height = 6, units = "in", res = 800)
print(umap_plot_with_line)
dev.off()

#TRAJ 4
traj <- 4
plot_df <- data.frame(
  UMAP_1 = reducedDim(sce_slingshot, 'UMAP')[, 1],
  UMAP_2 = reducedDim(sce_slingshot, 'UMAP')[, 2],
  pseudotime = slingPseudotime(sce_slingshot)[, traj]
)
curves <- slingCurves(sce_slingshot)
curve1_df <- as.data.frame(curves[[traj]]$s) # Get coordinates for the first curve
colnames(curve1_df) <- c("UMAP_1", "UMAP_2")

# 4. Generate the plot
umap_plot_with_line <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2)) +
  # Plot cells, colored by pseudotime
  geom_point(aes(color = pseudotime), size = 0.01, alpha = 0.8) +
  # Add the trajectory line ----------------------------------- NEW
  geom_path(data = curve1_df, aes(x = UMAP_1, y = UMAP_2), size = 1.2, color = "black") +
  # Use a nice color scale for continuous data
  scale_color_viridis(option = "magma", na.value = "lightgrey") +
  labs(
    title = "Slingshot Pseudotime on UMAP",
    subtitle = paste0("Trajectory ",traj," with Principal Curve", sep = ""),
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

# 5. Display and save the plot
png(filename = "~/SPATIAL/Maria_CosMx/Spatial_liver/figures/outs/traj4.png", width = 9.2, height = 6, units = "in", res = 800)
print(umap_plot_with_line)
dev.off()

## sEU EXPRESSION
umap
seu$cell_names


library(Seurat)


umap_ordered <- umap[match(colnames(seu), umap$cell_names), ]
umap_matrix <- as.matrix(umap_ordered[, c("UMAP1", "UMAP2")])
rownames(umap_matrix) <- colnames(seu)
colnames(umap_matrix) <- c("UMAP_1", "UMAP_2") 
seu[["umap"]] <- CreateDimReducObject(
  embeddings = umap_matrix,
  key = "UMAP_",
  assay = DefaultAssay(seu)
)
seu$leiden_res0.2 <- umap_ordered$leiden_res0.2
seu$etiology <- umap_ordered$etiology
seu <- NormalizeData(seu)
library(Seurat)
library(ggplot2)

# 1. Define your signatures
signatures <- list(
  Healthy_parenchyma = c("APOA1", "APOC1", "APOE", "TTR", "HPGDS", "FABP4", "CYP1B1", "INSR", "IGF1", "IGF2", "GLUL", "SLC2A1", "HSD17B2"),
  Immune_activation = c("B2M", "TAP1", "TAP2", "HLA-DRA", "HLA-DPA1", "CIITA", "STAT1", "IFIT3", "IFIT1", "ISG15", "MX1", "OAS1"),
  Antigen_presenting = c("CSF1R", "CLEC4A", "CD14", "CD1C", "TLR2", "IL15RA", "NOD2", "IDO1", "HLA-DPA1", "CXCL17", "CCL20"),
  Tissue_remodelling = c("ANGPT2", "COL16A1", "COL15A1", "TWIST1", "MMP7", "TIMP1", "TAGLN", "NR1H2", "NR2F2", "CXCL17"),
  Immune_tolerance = c("FOXP3", "IL10", "TGFB1", "CTLA4", "IL2RA", "ICOS", "CD274", "HAVCR2", "GATA3", "LGALS9", "STAT5B", "AHR", "TGFBR2", "IL1RN", "C1QA"),
  Immune_exhaustion = c("PDCD1", "LAG3", "TIGIT", "TOX", "EOMES", "CXCL13", "CD38", "HIF1A", "BCL2", "GZMK", "IL7R", "STAT3", "NR3C1", "IRF4"),
  Immune_suppression = c("ARG1", "IDO1", "IL10RA", "IL10RB", "TGFB2", "TGFBR1", "MERTK", "CSF1R", "CD163", "MRC1", "VCAN", "CSTB", "S100A9", "S100A8", "LGALS3")
)

out_dir <- "~/SPATIAL/Maria_CosMx/Spatial_liver/figures/outs/"

# 3. Loop through, calculate scores, plot, and save
for (sig_name in names(signatures)) {
  
  genes_to_plot <- signatures[[sig_name]]
  
  seu <- AddModuleScore(
    object = seu,
    features = list(genes_to_plot),
    name = sig_name,
    ctrl = 25 
  )
  
  target_feature <- paste0(sig_name, "1")
  
  # --- THE FIXED PLOT ---
  p <- FeaturePlot(
    object = seu, 
    features = target_feature, 
    pt.size = 0.05,        
    order = TRUE,          
    raster = FALSE,        
    min.cutoff = "q80",    # THE FIX: Forces the bottom 80% of cells to stay pure grey!
    max.cutoff = "q99"     # THE FIX: Only the top 1% get the absolute maximum red!
  ) + 
    scale_colour_gradient(
      low = "grey90",      
      high = "#FF0033"     
    ) + 
    theme(
      aspect.ratio = 1,
      axis.line = element_blank(),   
      axis.text = element_blank(),   
      axis.ticks = element_blank(),  
      axis.title = element_blank()   
    ) +
    ggtitle(gsub("_", " ", sig_name)) 
  
  file_out <- paste0(out_dir, "UMAP_", sig_name, ".png")
  png(filename = file_out, width = 9.2, height = 6, units = "in", res = 1200)
  print(p)
  dev.off()
}








## Stackplots PART B -------------------------------------------------------
library(dplyr)
library(ggplot2)

# 1. Map the clusters to your new broad groups
seu_neigh <- seu_neigh %>%
  mutate(
    broad_group = case_when(
      leiden_res0.2 %in% c(11, 6, 2) ~ "healthy/recovery",
      leiden_res0.2 %in% c(0, 3, 4, 7, 8, 10, 5) ~ "immune_active",
      leiden_res0.2 %in% c(1, 9) ~ "anti-inflamatory",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(broad_group)) # Drops any cells not in these groups

# 2. Calculate OVERALL proportions PER ETIOLOGY (No sample_id needed!)
df_bar <- seu_neigh %>%
  group_by(etiology, broad_group) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  group_by(etiology) %>%
  mutate(
    total_cells = sum(cell_count),
    proportion = cell_count / total_cells
  ) %>%
  ungroup()

# 3. Set your factor levels 
df_bar$etiology <- factor(df_bar$etiology, levels = c("HC", "HBV", "HDV RNA+", "HDV RNA-"))
df_bar$broad_group <- factor(df_bar$broad_group, levels = c("healthy/recovery", "immune_active", "anti-inflamatory"))

# 4. Define colors
group_colors <- c(
  "healthy/recovery" = "#4daf4a",   
  "immune_active"    = "#e41a1c",   
  "anti-inflamatory" = "#377eb8"    
)

# 5. Build the Grouped Bar Plot!
p <- ggplot(df_bar, aes(x = etiology, y = proportion, fill = broad_group)) +
  # Use position_dodge to put the bars side-by-side instead of stacked
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, color = "black") +
  scale_fill_manual(values = group_colors) +
  labs(x = "Etiology", y = "Overall Proportion of Cells", fill = "Cluster Group") +
  theme_classic() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
    text = element_text(size = 25),
    axis.line = element_line(linewidth = 1.5),
    axis.ticks = element_line(linewidth = 1.2),
    axis.ticks.length = unit(0.3, "cm")
  )

# 6. Save using your high-res specifications
file_out <- "~/SPATIAL/Maria_CosMx/Spatial_liver/figures/outs/barplot_overall_proportions_fig6.png"
png(filename = file_out, width = 15, height = 9, units = "in", res = 1200)
print(p)
dev.off()


#sTACKPLOT ---------------------------------------------------------------------
library(dplyr)
library(ggplot2)

# 1. Map the clusters to your new broad groups
seu_neigh <- seu_neigh %>%
  mutate(
    broad_group = case_when(
      leiden_res0.2 %in% c(11, 6, 2) ~ "healthy/recovery",
      leiden_res0.2 %in% c(0, 3, 4, 7, 8, 10, 5) ~ "immune_active",
      leiden_res0.2 %in% c(1, 9) ~ "anti-inflamatory",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(broad_group)) # Drops any cells not in these groups

# 2. Calculate OVERALL proportions PER ETIOLOGY (No sample_id needed!)
df_bar <- seu_neigh %>%
  group_by(etiology, broad_group) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  group_by(etiology) %>%
  mutate(
    total_cells = sum(cell_count),
    proportion = cell_count / total_cells
  ) %>%
  ungroup()

# 3. Set your factor levels 
df_bar$etiology <- factor(df_bar$etiology, levels = c("HC", "HBV", "HDV RNA+", "HDV RNA-"))
df_bar$broad_group <- factor(df_bar$broad_group, levels = c("healthy/recovery", "immune_active", "anti-inflamatory"))

# 4. Define colors
group_colors <- c(
  "healthy/recovery" = "#4daf4a",   
  "immune_active"    = "#e41a1c",   
  "anti-inflamatory" = "#377eb8"    
)

# 5. Build the Grouped Bar Plot!
p <- ggplot(df_bar, aes(x = etiology, y = proportion, fill = broad_group)) +
  # Use position_dodge to put the bars side-by-side instead of stacked
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, color = "black") +
  scale_fill_manual(values = group_colors) +
  labs(x = "Etiology", y = "Overall Proportion of Cells", fill = "Cluster Group") +
  theme_classic() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
    text = element_text(size = 25),
    axis.line = element_line(linewidth = 1.5),
    axis.ticks = element_line(linewidth = 1.2),
    axis.ticks.length = unit(0.3, "cm")
  )

# 6. Save using your high-res specifications
file_out <- "~/SPATIAL/Maria_CosMx/Spatial_liver/figures/outs/barplot_overall_proportions_fig6.png"
png(filename = file_out, width = 9.2, height = 6, units = "in", res = 1200)
print(p)
dev.off()

## Volcano Hepatos HDV RNA - vs HDV RNA +
# Volcano
seu@meta.data <- meta
png(filename = "~/SPATIAL/Maria_CosMx/Spatial_liver/figures/outs/volcano_HDVRNAneg_HDVRNApos.png",width = 11,height = 8,units = "in",res = 1200)
volcano(anot = "wide",ct = "Hepatocytes",dif_col = "etiology",seu = seu,id1 = "HDV RNA-",id2 = "HDV RNA+")[[1]]
dev.off()
