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


#Polygons Slide 2 FOV 2, Showing subset  HDV RNA+ 
Slide2_pols <- read_csv("/home/mmoro/SPATIAL/Maria_CosMx/Maria_v2/Polygons/Slide_2.csv")
Slide2_pols <- Slide2_pols[Slide2_pols$fov == "2",]
meta_cut <- meta[meta$tissue == "Slide_2" & meta$fov == 2,]
Slide2_pols <- Slide2_pols[Slide2_pols$cell_names %in% meta_cut$cell_names,]
Slide2_pols[["wide"]] <- plyr::mapvalues(x = Slide2_pols[["cell_names"]],
                                                       from = meta_cut[["cell_names"]],
                                                       to = meta_cut[["wide"]])
Slide2_pols[["wide"]] <- as.factor(Slide2_pols[["wide"]])

# Control_530_pols <- Control_530_pols %>%
#   mutate(reduced_clust = if_else(reduced_clust %in% c("12", "20"), reduced_clust, "Other"))
# resols_neigh["Other"] <- "#393939"

png(filename = "~/SPATIAL/Maria_CosMx/Spatial_liver/figures/outs/cslide2_fov2_wide.png",width = 8,height = 8,units = "in",res = 1200)
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


#Polygons Slide 4 FOV 4, Showing subset
Slide1_pols <- read_csv("/home/mmoro/SPATIAL/Maria_CosMx/Maria_v2/Polygons/Slide_1.csv")
Slide1_pols <- Slide1_pols[Slide1_pols$fov == "1",]
meta_cut <- meta[meta$tissue == "Slide_1" & meta$fov == 1,]
Slide1_pols <- Slide1_pols[Slide1_pols$cell_names %in% meta_cut$cell_names,]
Slide1_pols[["wide"]] <- plyr::mapvalues(x = Slide1_pols[["cell_names"]],
                                           from = meta_cut[["cell_names"]],
                                           to = meta_cut[["wide"]])
Slide1_pols[["wide"]] <- as.factor(Slide1_pols[["wide"]])

# Control_530_pols <- Control_530_pols %>%
#   mutate(reduced_clust = if_else(reduced_clust %in% c("12", "20"), reduced_clust, "Other"))
# resols_neigh["Other"] <- "#393939"

png(filename = "~/SPATIAL/Maria_CosMx/Spatial_liver/figures/outs/cslide1_fov1_wide.png",width = 8,height = 8,units = "in",res = 1200)
p <-
  ggplot2::ggplot(Slide1_pols, ggplot2::aes(x = .data[["x_global_px"]], y = .data[["y_global_px"]])) +
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


#Polygons Slide 2 FOV 16, Showing subset  HBV
Slide2_pols <- read_csv("/home/mmoro/SPATIAL/Maria_CosMx/Maria_v2/Polygons/Slide_2.csv")
Slide2_pols <- Slide2_pols[Slide2_pols$fov == "16",]
meta_cut <- meta[meta$tissue == "Slide_2" & meta$fov == 16,]
Slide2_pols <- Slide2_pols[Slide2_pols$cell_names %in% meta_cut$cell_names,]
Slide2_pols[["wide"]] <- plyr::mapvalues(x = Slide2_pols[["cell_names"]],
                                         from = meta_cut[["cell_names"]],
                                         to = meta_cut[["wide"]])
Slide2_pols[["wide"]] <- as.factor(Slide2_pols[["wide"]])

# Control_530_pols <- Control_530_pols %>%
#   mutate(reduced_clust = if_else(reduced_clust %in% c("12", "20"), reduced_clust, "Other"))
# resols_neigh["Other"] <- "#393939"

png(filename = "~/SPATIAL/Maria_CosMx/Spatial_liver/figures/outs/cslide2_fov16_wide.png",width = 8,height = 8,units = "in",res = 1200)
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

# 
# 
# ## refineD
# #Polygons Slide 2 FOV 2, Showing refined
Slide2_pols <- read_csv("/home/mmoro/SPATIAL/Maria_CosMx/Maria_v2/Polygons/Slide_2.csv")
Slide2_pols <- Slide2_pols[Slide2_pols$fov == "2",]
meta_cut <- meta[meta$tissue == "Slide_2" & meta$fov == 2,]
Slide2_pols <- Slide2_pols[Slide2_pols$cell_names %in% meta_cut$cell_names,]
Slide2_pols[["refined"]] <- plyr::mapvalues(x = Slide2_pols[["cell_names"]],
                                           from = meta_cut[["cell_names"]],
                                           to = meta_cut[["refined"]])
Slide2_pols[["refined"]] <- as.factor(Slide2_pols[["refined"]])

Slide2_pols <- Slide2_pols %>%
  mutate(refined = if_else(refined %in% c("Hepatocyte_1" , "Hepatocyte_2" ,"Hepatocyte_3" ,"Hepatocyte_5" ,"Hepatocyte_5"), refined, "Other"))
refined_col["Other"] <- "#393939"

png(filename = "~/SPATIAL/Maria_CosMx/Spatial_liver/figures/outs/cslide2_fov2_refined.png",width = 8,height = 8,units = "in",res = 1200)
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


## refineD hepatos
#Polygons Slide 2 FOV 16, Showing refined
Slide2_pols <- read_csv("/home/mmoro/SPATIAL/Maria_CosMx/Maria_v2/Polygons/Slide_2.csv")
Slide2_pols <- Slide2_pols[Slide2_pols$fov == "16",]
meta_cut <- meta[meta$tissue == "Slide_2" & meta$fov == 16,]
Slide2_pols <- Slide2_pols[Slide2_pols$cell_names %in% meta_cut$cell_names,]
Slide2_pols[["refined"]] <- plyr::mapvalues(x = Slide2_pols[["cell_names"]],
                                            from = meta_cut[["cell_names"]],
                                            to = meta_cut[["refined"]])
Slide2_pols[["refined"]] <- as.factor(Slide2_pols[["refined"]])

Slide2_pols <- Slide2_pols %>%
  mutate(refined = if_else(refined %in% c("Hepatocyte_1" , "Hepatocyte_2" ,"Hepatocyte_3" ,"Hepatocyte_5" ,"Hepatocyte_5"), refined, "Other"))
refined_col["Other"] <- "#393939"

png(filename = "~/SPATIAL/Maria_CosMx/Spatial_liver/figures/outs/cslide2_fov16_refined.png",width = 8,height = 8,units = "in",res = 1200)
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


## refineD but myeloids
#Polygons Slide 1 FOV 1, Showing refined
Slide2_pols <- read_csv("/home/mmoro/SPATIAL/Maria_CosMx/Maria_v2/Polygons/Slide_1.csv")
Slide2_pols <- Slide2_pols[Slide2_pols$fov == "1",]
meta_cut <- meta[meta$tissue == "Slide_1" & meta$fov == 1,]
Slide2_pols <- Slide2_pols[Slide2_pols$cell_names %in% meta_cut$cell_names,]
Slide2_pols[["refined"]] <- plyr::mapvalues(x = Slide2_pols[["cell_names"]],
                                            from = meta_cut[["cell_names"]],
                                            to = meta_cut[["refined"]])
Slide2_pols[["refined"]] <- as.factor(Slide2_pols[["refined"]])

Slide2_pols <- Slide2_pols %>%
  mutate(refined = if_else(refined %in% c("M2_LYVE1",  "KC1"    ,   "KC2"    ,   "Monocytes" ,"DCs_CD1C" , "M1"  ), refined, "Other"))
refined_col["Other"] <- "#393939"

png(filename = "~/SPATIAL/Maria_CosMx/Spatial_liver/figures/outs/cslide1_fov1_refined_myel.png",width = 8,height = 8,units = "in",res = 1200)
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


## REFINED BUT MYELOIDSS
## refineD
#Polygons Slide 2 FOV 16, Showing refined
Slide2_pols <- read_csv("/home/mmoro/SPATIAL/Maria_CosMx/Maria_v2/Polygons/Slide_2.csv")
Slide2_pols <- Slide2_pols[Slide2_pols$fov == "2",]
meta_cut <- meta[meta$tissue == "Slide_2" & meta$fov == 2,]
Slide2_pols <- Slide2_pols[Slide2_pols$cell_names %in% meta_cut$cell_names,]
Slide2_pols[["refined"]] <- plyr::mapvalues(x = Slide2_pols[["cell_names"]],
                                            from = meta_cut[["cell_names"]],
                                            to = meta_cut[["refined"]])
Slide2_pols[["refined"]] <- as.factor(Slide2_pols[["refined"]])

Slide2_pols <- Slide2_pols %>%
  mutate(refined = if_else(refined %in%c("M2_LYVE1",  "KC1"    ,   "KC2"    ,   "Monocytes" ,"DCs_CD1C" , "M1"  ), refined, "Other"))
refined_col["Other"] <- "#393939"

png(filename = "~/SPATIAL/Maria_CosMx/Spatial_liver/figures/outs/cslide2_fov2_refined_myel.png",width = 8,height = 8,units = "in",res = 1200)
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


## refineD
#Polygons Slide 1 FOV 1, Showing refined
Slide2_pols <- read_csv("/home/mmoro/SPATIAL/Maria_CosMx/Maria_v2/Polygons/Slide_1.csv")
Slide2_pols <- Slide2_pols[Slide2_pols$fov == "1",]
meta_cut <- meta[meta$tissue == "Slide_1" & meta$fov == 1,]
Slide2_pols <- Slide2_pols[Slide2_pols$cell_names %in% meta_cut$cell_names,]
Slide2_pols[["refined"]] <- plyr::mapvalues(x = Slide2_pols[["cell_names"]],
                                            from = meta_cut[["cell_names"]],
                                            to = meta_cut[["refined"]])
Slide2_pols[["refined"]] <- as.factor(Slide2_pols[["refined"]])

Slide2_pols <- Slide2_pols %>%
  mutate(refined = if_else(refined %in% c("M2_LYVE1",  "KC1"    ,   "KC2"    ,   "Monocytes" ,"DCs_CD1C" , "M1"  ), refined, "Other"))
refined_col["Other"] <- "#393939"

png(filename = "~/SPATIAL/Maria_CosMx/Spatial_liver/figures/outs/cslide1_fov1_refined.png",width = 8,height = 8,units = "in",res = 1200)
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



## Perform pathway analysis
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
  deg_results$diffexpressed[deg_results$avg_log2FC > log2(1.2) & deg_results$p_val < 0.05] <- "p.val<0.05 & FC>1.2"
  deg_results$diffexpressed[deg_results$avg_log2FC < -log2(1.2) & deg_results$p_val < 0.05] <- "p.val<0.05 & FC<0.83"
  deg_results$diffexpressed[deg_results$avg_log2FC > log2(1.2) & deg_results$p_val_adj < 0.05] <- "p.adj<0.05 & FC>1.2"
  deg_results$diffexpressed[deg_results$avg_log2FC < -log2(1.2) & deg_results$p_val_adj < 0.05] <- "p.adj<0.05 & FC<0.83"
  #Label of most down-up regulated genes
  deg_results$delabel <- NA
  deg_results$delabel[deg_results$diffexpressed != "NO"] <- deg_results$genes[deg_results$diffexpressed != "NO"]
  deg_results$p_val <- ifelse(deg_results$p_val < 1e-300, 1e-300, deg_results$p_val)
  
  #Plot
  p <- ggplot(data = deg_results, aes(x = avg_log2FC, y = -log10(p_val), col = diffexpressed, label = delabel)) + 
    geom_point() + 
    theme_bw() +
    geom_text_repel() +
    scale_color_manual(values = c("p.val<0.05 & FC<0.83" = "green", "p.adj<0.05 & FC<0.83" = "darkgreen", "p.val<0.05 & FC>1.2" = "red", "p.adj<0.05 & FC>1.2" = "darkred")) +
    geom_vline(xintercept=c(-log2(1.2),log2(1.2) ), col="black",linetype = "dashed") +
    geom_hline(yintercept=-log10(0.05), col="black",linetype = "dashed") +
    theme(text = element_text(size = 18)) +
    ggtitle(paste0(ct, " : ",id1, " vs ", id2, sep = ""))
  
  # Print the plot
  return(list(p,deg_results))
  
}

msigdb_hallmark <- msigdbr(species = "Homo sapiens", category = "C5" , subcategory = "BP")

#HC vs HBV
hc_vs_hbv <- volcano(anot = "new_anot", ct = "Hepatocytes", seu = seu, dif_col = "etiology",id1 = "HDV RNA+", id2 = "HBV")[[2]]


up <- hc_vs_hbv[hc_vs_hbv$diffexpressed %in%  c("p.adj<0.05 & FC>1.2","p.val<0.05 & FC>1.2"),]$genes
gene_list_up <- bitr(up, fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)

down <- hc_vs_hbv[hc_vs_hbv$diffexpressed %in% c("p.adj<0.05 & FC<0.83","p.val<0.05 & FC<0.83"),]$genes
gene_list_down <- bitr(down, fromType = "SYMBOL",
                       toType = "ENTREZID",
                       OrgDb = org.Hs.eg.db)

gene_sets <- msigdb_hallmark %>%
  dplyr::select(gs_name, entrez_gene)

# Perform enrichment analysis
up_enrichment_results <- as.data.frame(enricher(gene = gene_list_up$ENTREZID,
                                                TERM2GENE = gene_sets))
up_enrichment_results <- up_enrichment_results %>%
  mutate(GeneRatio = as.numeric(sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))))


down_enrichment_results <- as.data.frame(enricher(gene = gene_list_down$ENTREZID,
                                                  TERM2GENE = gene_sets))

down_enrichment_results <- down_enrichment_results %>%
  mutate(GeneRatio = as.numeric(sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))))


up_enrichment_results <- up_enrichment_results[order(up_enrichment_results$qvalue), ]

# # Order the dataframe by the q_value column in ascending order
down_enrichment_results <- down_enrichment_results[order(down_enrichment_results$qvalue), ]

up_enrichment_results$s1 <- "Upregulated"
down_enrichment_results$s1 <- "Downregulated"

final_enrichment_result <- rbind(up_enrichment_results,down_enrichment_results)

final_enrichment_result$log10pval <- -log10(final_enrichment_result$pvalue)

final_enrichment_result <- final_enrichment_result %>%
  arrange(GeneRatio)
final_enrichment_result1 <- final_enrichment_result[final_enrichment_result$ID %in% pathways,]


  
library(ggplot2)

p <- ggplot(final_enrichment_result1, aes(x = s1 , y = Description, size = GeneRatio, color = s1)) +
  geom_point() +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +
  scale_size_continuous(name = "Gene ratio", range = c(3, 10),
                        guide = guide_legend(override.aes = list(color = "black", fill = "black"))) +
  theme_bw() +
  labs(title = "", x = "", y = "Pathway description", size = "Gene ratio", color = "Pathway")+
  theme(
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 15),
    legend.text = element_text(size = 15),
    legend.title = element_text(face = "bold", size = 15),
    
    # use size for older ggplot2, linewidth also works on new versions
    axis.line  = element_line(size = 1.2, colour = "black"),
    axis.ticks = element_line(size = 1.2, colour = "black"),
    axis.ticks.length = unit(6, "pt"),
    
    # Make sure the panel background is transparent so the border is visible
    panel.background = element_rect(fill = NA, colour = NA),
    panel.border     = element_rect(colour = "black", fill = NA, size = 1.2)
  )
png(filename = "~/SPATIAL/Maria_CosMx/Spatial_liver/figures/outs/cslide2_fov16_wide.png",width = 11,height = 8,units = "in",res = 1200)
p
dev.off()

write.csv(final_enrichment_result1,"~/SPATIAL/Maria_CosMx/Spatial_liver/figures/outs/pathway_cosmx.csv")