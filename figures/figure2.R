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

## Functions ----------------------------------
align_image <- function(patient, fov) {
  # Mapping
  dig_alignment <- c("Slide_1" = "S1", 
                     "Slide_2" = "S2", 
                     "Slide_3" = "S3")
  
  # Load image
  img_path <- paste0(
    "/home/mmoro/SPATIAL/Maria_CosMx/Maria_v1/Post_analysis/Aligment/aligment/",
    dig_alignment[patient], "_F", fov, ".png"
  )
  
  img <- png::readPNG(img_path, native = FALSE)
  img_dims <- dim(img)[1:2]   # height, width
  
  # Scaling parameters depending on slide
  if (patient == "Slide_3") {
    x_scale <- 0.149; y_scale <- 0.13
    x_shift <- 2;     y_shift <- 43
    
  } else if (patient == "Slide_2") {
    x_scale <- 0.102; y_scale <- 0.146
    x_shift <- 98;    y_shift <- 0
    
  } else if (patient == "Slide_1") {
    x_scale <- 0.15;  y_scale <- 0.13
    x_shift <- 8;     y_shift <- 50
  }
  
  # Adjust coordinates
  xmin_adj <- (-x_shift) / x_scale
  xmax_adj <- (img_dims[2] - x_shift) / x_scale
  ymin_adj <- (-y_shift) / y_scale
  ymax_adj <- (img_dims[1] - y_shift) / y_scale
  
  # Return plotted image
  ggplot2::ggplot() + 
    ggplot2::annotation_custom(
      grid::rasterGrob(img),
      xmin = xmin_adj,
      xmax = xmax_adj,
      ymin = ymin_adj,
      ymax = ymax_adj
    )
}

 
#Volcano function
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



## Figure 3A ----------------------------------

#HDAg+/-	
#S2F8 (for D+)
final_pal <- c(
  "none" = "#858383",
  "+" = "#B47846",
  "-" = "steelblue"
)

meta$Type_D <- ifelse(is.na(meta$Type_D), yes = "none",no = meta$Type_D )
alig <- align_image(patient = "Slide_2", fov = 8)
slide = "Slide_2"
fov = "8"
pols <- read_csv(paste0("/home/mmoro/SPATIAL/Maria_CosMx/Maria_v1/Post_analysis/Selection_regions/excels_polygons_HDAg/polygons_","Slide_2","_fov_","8",".csv",sep = ""))
library(ggforce)
pols$group <- ifelse(test = pols$polygon_id == 4,yes =  "-",no = "+" )

p <- alig +
  ggplot2::coord_fixed() +
  geom_bspline_closed(
    data = pols,
    aes(x = x, y = y, group = polygon_id, fill = group),
    color = "black",
    size = 1,
    alpha = 0.15
  ) +
  scale_fill_manual(values = final_pal) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    axis.line  = element_blank(), legend.position = "none"
  )

png(filename = "~/SPATIAL/Maria_CosMx/Spatial_liver/figures/outs/fig3B_S2F8_alig.png",width = 8,height = 8,units = "in",res = 800)
p
dev.off()


#Plot FOV!
meta_cut <- meta[meta$tissue == "Slide_2" & meta$fov == 8,]
cell_pols <- read.csv("/home/mmoro/SPATIAL/Maria_CosMx/Maria_v2/Polygons/Slide_2.csv")
cells_pols <- cell_pols[cell_pols$fov == 8,]
type_d <- meta_cut$Type_D
names(type_d) <- meta_cut$cell_names
cells_pols$type_d <- type_d[cells_pols$cell_names]

p <- plot.spatial(
  meta = meta,
  polygon_file = read.csv(
    paste0(
      "/home/mmoro/SPATIAL/Maria_CosMx/Maria_v1/Polygons/",
      "Slide_2",
      ".csv",
      sep = ""
    )
  ),
  x = "x_local_px",
  y = "y_local_px",
  col_fov = "fov",
  fov =  8,
  palette = final_pal,
  color = "Type_D",
  per = "black",
  alpha = 1,
  dark = F,
  aligment = alig,
  ptsize = 0.1
) 

png(filename = "~/SPATIAL/Maria_CosMx/Spatial_liver/figures/outs/fig3B_S2F8.png",width = 10,height = 8,units = "in",res = 800)
p
dev.off()


#Volcanos

#Volcano function
seu@meta.data$wide <-
  ifelse(
    test = seu@meta.data$subset == "Im_hepatocytes",
    yes = "Hepatocytes",
    no = seu@meta.data$subset
  )

png(filename = "~/SPATIAL/Maria_CosMx/Spatial_liver/figures/outs/volc_typeD.png",width = 10,height = 7,units = "in",res = 800)
volcano(anot = "new_anot", ct = "Hepatocytes", dif_col = "Type_D", seu = seu,id1 = "+",id2 = "-")
dev.off()


##Figure 2C
#S3F3 (for S+)
#Plot FOV!

#META
meta <- seu@meta.data
alig <- align_image(patient = "Slide_3", fov = 3)
meta_cut <- meta[meta$tissue == "Slide_3" & meta$fov == 3,]

p <- plot.spatial(
  meta = meta,
  polygon_file = read.csv(
    paste0(
      "/home/mmoro/SPATIAL/Maria_CosMx/Maria_v1/Polygons/",
      "Slide_3",
      ".csv",
      sep = ""
    )
  ),
  x = "x_local_px",
  y = "y_local_px",
  col_fov = "fov",
  fov =  3,
  palette = wide,
  color = "wide",
  per = "black",
  alpha = 0,
  dark = F,
  aligment = alig,
  ptsize = 0.1
) 

png(filename = "~/SPATIAL/Maria_CosMx/Spatial_liver/figures/outs/fig3B_S3F3.png",width = 10,height = 8,units = "in",res = 600)
p
dev.off()


#S3F1 (for S-)
alig <- align_image(patient = "Slide_3", fov = 1)
meta_cut <- meta[meta$tissue == "Slide_3" & meta$fov == 1,]

p <- plot.spatial(
  meta = meta,
  polygon_file = read.csv(
    paste0(
      "/home/mmoro/SPATIAL/Maria_CosMx/Maria_v1/Polygons/",
      "Slide_3",
      ".csv",
      sep = ""
    )
  ),
  x = "x_local_px",
  y = "y_local_px",
  col_fov = "fov",
  fov =  1,
  palette = wide,
  color = "wide",
  per = "black",
  alpha = 0,
  dark = F,
  aligment = alig,
  ptsize = 0.1
) 
png(filename = "~/SPATIAL/Maria_CosMx/Spatial_liver/figures/outs/fig3C_S3F1.png",width = 10,height = 8,units = "in",res = 600)
p
dev.off()


#Volcanos EXPLANTS
explant <- seu[,rownames(seu@meta.data[seu@meta.data$tissue == "Slide_3" & seu@meta.data$fov %in% c(1,4,11,16,24,25,13,2,3,5,6,12,14,23),])]

png(filename = "~/SPATIAL/Maria_CosMx/Spatial_liver/figures/outs/volc_typeS.png",width = 10,height = 7,units = "in",res = 600)
volcano(anot = "new_anot", ct = "Hepatocytes", dif_col = "Type_S", seu = seu,id1 = "+",id2 = "-")
dev.off()


## Figure 2E -----------------------------------
genes <- c(
  "IGHG1", "CXCL9", "IGKC", "CD74", "IL32", "CXCL10", "STAT1", "SQLE", 
  "RARRES2", "FGFR2", "MYL12A", "SLPI", "IGHG2",
  "CCND1", "CLU", "H4C3", "SOD2", "IFITM3", "SOX2", "IRF3", "CTSD", "CRP", 
  "NPR2", "ACKR4", "ADGRV1", "CD14", "COL16A1", "MAP2K1", "NLRP1", "MYH11", 
  "CD3G", "FHIT", "CSF1R", "KRT1", "TLR2", "KRT86", "MTOR", "PLAC8", "CLOCK", 
  "TCF7", "ADGRE2", "CUZD1", "TGFBR1", "CFD", "FLT1", "IL7R"
)

genes <- c("IL32", "IRF3", "IFITM3", "SOD2", "CXCL9", 
               "HLA-DPB1", "CRP", "ACKR4", "ACTA2", "NLRP1")

genes <- c("IL32", "CXCL9")
### Density plot
unique(seu$Type_D)

type_d <- seu[,seu$Type_D %in% c("+","-") & seu$new_anot == "Hepatocytes" ]
raw_counts_mat <- t(as.matrix(GetAssayData(type_d, layer = "counts")[genes, ]))
df_counts <- as.data.frame(raw_counts_mat)
df_counts$Condition <- type_d$Type_D
plot_data <- df_counts %>%
  group_by(Condition) %>%
  summarise(across(all_of(genes), mean)) %>% 
  pivot_longer(cols = -Condition, names_to = "Gene", values_to = "Mean_Expression")
plot_data$Gene <- factor(plot_data$Gene, levels = genes)

ggplot(plot_data, aes(x = Gene, y = Mean_Expression, fill = Condition, group = Condition, color = Condition)) +
  geom_area(alpha = 0.4, position = "identity") +
  geom_line(size = 1) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(title = "Mean Raw Count Expression Profile",
       x = "Genes",
       y = "Mean Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# S2F10 (for D + )

g <- genes[1]
seu <- NormalizeData(seu,assay = "RNA")

for (g in genes) {
  
  message("Plotting ", g)
  
  png(
    filename = paste0("~/SPATIAL/Maria_CosMx/Spatial_liver/figures/outs/", g, ".png"),
    width = 5, height = 7, units = "in", res = 100
  )
  
  p <- plot.gene.spatial(
    meta = meta,
    expression_matrix = as.matrix(seu@assays$RNA$data),
    genes = g,
    polygon_file = read.csv(
      paste0(
        "/home/mmoro/SPATIAL/Maria_CosMx/Maria_v1/Polygons/Slide_2.csv"
      )
    ),
    x = "x_local_px",
    y = "y_local_px",
    col_fov = "fov",
    fov = 10,
    norm = FALSE, viridis = "D"
  )
  
  print(p)   # ← THIS is required inside a PNG device
  
  dev.off()
}


# S2F11 (for D -)

for (g in genes) {
  
  message("Plotting ", g)
  
  png(
    filename = paste0("~/SPATIAL/Maria_CosMx/Spatial_liver/figures/outs/", g, ".png"),
    width = 5, height = 7, units = "in", res = 100
  )
  
  p <- plot.gene.spatial(
    meta = meta,
    expression_matrix = as.matrix(seu@assays$RNA$data),
    genes = g,
    polygon_file = read.csv(
      paste0(
        "/home/mmoro/SPATIAL/Maria_CosMx/Maria_v1/Polygons/Slide_2.csv"
      )
    ),
    x = "x_local_px",
    y = "y_local_px",
    col_fov = "fov",
    fov = 11,
    norm = FALSE
  )
  
  print(p)   # ← THIS is required inside a PNG device
  
  dev.off()
}


## Genes
# S3F3 (for S + )


genes2 <- c(
  "CXCL8", "CHI3L1", "CRP", "RGS2", "SOD2", "IGFBP7", "TNFRSF12A", "VIM", 
  "BCL2L1", "TPM1", "LYZ", "GLUL", "SAT1", "CDKN1A", "NCAM1", 
  "PTGES", "KRT23", "IL32", "MIF", "DCN", "FN1", "CD24", "ATR", "TNFSF12", 
  "TTR", "ARG1", "IFI27", "IL17RA", "ENG", "HMGCS1", 
  "ADGRA3", "MFAP5", "CASP3", "HSP90AA1", "HBB", "CLEC1A", "CXCL14", 
  "DLL1", "OASL", "FGF2", "WNT3", "CCR1", "DUSP4", "FZD8", "DHRS2", "FASN"
)
genes2 <- c("CXCL8", "GLUL", "SOD2", "IGHG2", "SAT1", 
            "CHI3L1", "IFI27", "TTR", "IL17RA", "ACTA2")
type_s <- seu[,seu$Type_S %in% c("+","-") & seu$new_anot == "Hepatocytes" ]
raw_counts_mat <- t(as.matrix(GetAssayData(type_s, layer = "counts")[genes2, ]))
df_counts <- as.data.frame(raw_counts_mat)
df_counts$Condition <- type_s$Type_S
plot_data <- df_counts %>%
  group_by(Condition) %>%
  summarise(across(all_of(genes2), mean)) %>% 
  pivot_longer(cols = -Condition, names_to = "Gene", values_to = "Mean_Expression")
plot_data$Gene <- factor(plot_data$Gene, levels = genes2)

ggplot(plot_data, aes(x = Gene, y = Mean_Expression, fill = Condition, group = Condition, color = Condition)) +
  geom_area(alpha = 0.4, position = "identity") +
  geom_line(size = 1) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(title = "Mean Raw Count Expression Profile",
       x = "Genes",
       y = "Mean Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



for (g in genes2) {
  
  message("Plotting ", g)
  
  png(
    filename = paste0("~/SPATIAL/Maria_CosMx/Spatial_liver/figures/outs/", g, ".png"),
    width = 5, height = 7, units = "in", res = 100
  )
  
  p <- plot.gene.spatial(
    meta = meta,
    expression_matrix = as.matrix(seu@assays$RNA$data),
    genes = g,
    polygon_file = read.csv(
      paste0(
        "/home/mmoro/SPATIAL/Maria_CosMx/Maria_v1/Polygons/Slide_3.csv"
      )
    ),
    x = "x_local_px",
    y = "y_local_px",
    col_fov = "fov",
    fov = 3,
    norm = FALSE
  )
  
  print(p)   # ← THIS is required inside a PNG device
  
  dev.off()
}


# S3F1 (for S - )

for (g in genes2) {
  
  message("Plotting ", g)
  
  png(
    filename = paste0("~/SPATIAL/Maria_CosMx/Spatial_liver/figures/outs/", g, ".png"),
    width = 5, height = 7, units = "in", res = 100
  )
  
  p <- plot.gene.spatial(
    meta = meta,
    expression_matrix = as.matrix(seu@assays$RNA$data),
    genes = g,
    polygon_file = read.csv(
      paste0(
        "/home/mmoro/SPATIAL/Maria_CosMx/Maria_v1/Polygons/Slide_3.csv"
      )
    ),
    x = "x_local_px",
    y = "y_local_px",
    col_fov = "fov",
    fov = 1,
    norm = FALSE
  )
  
  print(p)   # ← THIS is required inside a PNG device
  
  dev.off()
}


## Final REAL PLOTS 22/01/26 Figure 3D 
library(ggplot2)
library(cowplot) # Required to extract the legend separate from the plot
library(grid)

# --- 1. SETUP & THEMES ---

genes <- c("IL32", "CXCL9")
out_dir <- "~/SPATIAL/Maria_CosMx/Spatial_liver/figures/outs/"

# Load polygon file once
poly_data <- read.csv("/home/mmoro/SPATIAL/Maria_CosMx/Maria_v1/Polygons/Slide_2.csv")

# Define a consistent Dark Theme with NO AXIS
theme_dark_clean <- theme(
  plot.background = element_rect(fill = "black", color = NA),
  panel.background = element_rect(fill = "black", color = NA),
  legend.background = element_rect(fill = "black", color = NA),
  legend.text = element_text(color = "white"),
  legend.title = element_text(color = "white"),
  
  # REMOVE ALL AXIS ELEMENTS
  axis.text = element_blank(),
  axis.title = element_blank(),
  axis.ticks = element_blank(),
  axis.line = element_blank(),
  
  title = element_text(color = "white"),
  panel.grid = element_blank()
)

# --- 2. LOOP THROUGH GENES ---

for (g in genes) {
  
  message("Processing ", g)
  
  # A. CALCULATE COMMON LIMITS
  # Identify cells belonging to FOV 10 and 11 to set shared color scale
  cells_fov10 <- rownames(meta)[meta$fov == 10]
  cells_fov11 <- rownames(meta)[meta$fov == 11]
  all_cells <- c(cells_fov10, cells_fov11)
  
  # Extract expression data
  expr_values <- seu@assays$RNA$data[g, all_cells]
  common_limits <- c(min(expr_values), max(expr_values))
  
  
  # B. CREATE & SAVE LEGEND (Separate Figure)
  
  # Create temp plot to extract legend
  p_temp <- plot.gene.spatial(
    meta = meta,
    expression_matrix = as.matrix(seu@assays$RNA$data),
    genes = g,
    polygon_file = poly_data,
    x = "x_local_px", y = "y_local_px",
    col_fov = "fov", fov = 10, 
    norm = FALSE, viridis = "D"
  ) +
    scale_fill_viridis_c(option = "D", limits = common_limits) +
    scale_color_viridis_c(option = "D", limits = common_limits) +
    theme_dark_clean +
    theme(legend.position = "right")
  
  # Extract and Save Legend
  legend_plot <- get_legend(p_temp)
  
  png(filename = paste0(out_dir, g, "_Legend.png"), width = 2, height = 5, units = "in", res = 400)
  grid.newpage()
  grid.draw(legend_plot)
  dev.off()
  
  
  # C. PLOT FOV 10 (Clean, Dark, No Legend)
  
  p10 <- plot.gene.spatial(
    meta = meta,
    expression_matrix = as.matrix(seu@assays$RNA$data),
    genes = g,
    polygon_file = poly_data,
    x = "x_local_px", y = "y_local_px",
    col_fov = "fov", fov = 10,
    norm = FALSE, viridis = "D"
  ) +
    scale_fill_viridis_c(option = "D", limits = common_limits) +
    scale_color_viridis_c(option = "D", limits = common_limits) +
    theme_dark_clean +
    theme(legend.position = "none") +
    labs(x = NULL, y = NULL) # Explicitly remove labels just in case
  
  png(filename = paste0(out_dir, g, "_FOV10.png"), width = 7, height = 7, units = "in", res = 400)
  print(p10)
  dev.off()
  
  
  # D. PLOT FOV 11 (Clean, Dark, No Legend)
  
  p11 <- plot.gene.spatial(
    meta = meta,
    expression_matrix = as.matrix(seu@assays$RNA$data),
    genes = g,
    polygon_file = poly_data,
    x = "x_local_px", y = "y_local_px",
    col_fov = "fov", fov = 11,
    norm = FALSE, viridis = "D"
  ) +
    scale_fill_viridis_c(option = "D", limits = common_limits) +
    scale_color_viridis_c(option = "D", limits = common_limits) +
    theme_dark_clean +
    theme(legend.position = "none") +
    labs(x = NULL, y = NULL)
  
  png(filename = paste0(out_dir, g, "_FOV11.png"), width = 7, height = 7, units = "in", res = 400)
  print(p11)
  dev.off()
}



## IFI27 & APOA1 slide3 fov 1, fov2


## Final REAL PLOTS 22/01/26 Figure 3D 
library(ggplot2)
library(cowplot) # Required to extract the legend separate from the plot
library(grid)

# --- 1. SETUP & THEMES ---
meta <- seu@meta.data[seu@meta.data$tissue == "Slide_3",]
genes <- c("IFI27", "SOD2")
out_dir <- "~/SPATIAL/Maria_CosMx/Spatial_liver/figures/outs/"

# Load polygon file once
poly_data <- read.csv("/home/mmoro/SPATIAL/Maria_CosMx/Maria_v1/Polygons/Slide_3.csv")

# Define a consistent Dark Theme with NO AXIS
theme_dark_clean <- theme(
  plot.background = element_rect(fill = "black", color = NA),
  panel.background = element_rect(fill = "black", color = NA),
  legend.background = element_rect(fill = "black", color = NA),
  legend.text = element_text(color = "white"),
  legend.title = element_text(color = "white"),
  
  # REMOVE ALL AXIS ELEMENTS
  axis.text = element_blank(),
  axis.title = element_blank(),
  axis.ticks = element_blank(),
  axis.line = element_blank(),
  
  title = element_text(color = "white"),
  panel.grid = element_blank()
)

# --- 2. LOOP THROUGH GENES ---

for (g in genes) {
  
  message("Processing ", g)
  
  # A. CALCULATE COMMON LIMITS
  # Identify cells belonging to FOV 1 and 3 to set shared color scale
  cells_fov1 <- rownames(meta)[meta$fov == 1]
  cells_fov3 <- rownames(meta)[meta$fov == 3]
  all_cells <- c(cells_fov1, cells_fov3)
  
  # Extract expression data
  expr_values <- seu@assays$RNA$data[g, all_cells]
  common_limits <- c(min(expr_values), max(expr_values))
  
  
  # B. CREATE & SAVE LEGEND (Separate Figure)
  
  # Create temp plot to extract legend
  p_temp <- plot.gene.spatial(
    meta = meta,
    expression_matrix = as.matrix(seu@assays$RNA$data),
    genes = g,
    polygon_file = poly_data,
    x = "x_local_px", y = "y_local_px",
    col_fov = "fov", fov = 1, 
    norm = FALSE, viridis = "D"
  ) +
    scale_fill_viridis_c(option = "D", limits = common_limits) +
    scale_color_viridis_c(option = "D", limits = common_limits) +
    theme_dark_clean +
    theme(legend.position = "right")
  
  # Extract and Save Legend
  legend_plot <- get_legend(p_temp)
  
  png(filename = paste0(out_dir, g, "_Legend.png"), width = 2, height = 5, units = "in", res = 400)
  grid.newpage()
  grid.draw(legend_plot)
  dev.off()
  
  
  # C. PLOT FOV 1 (Clean, Dark, No Legend)
  
  p1 <- plot.gene.spatial(
    meta = meta,
    expression_matrix = as.matrix(seu@assays$RNA$data),
    genes = g,
    polygon_file = poly_data,
    x = "x_local_px", y = "y_local_px",
    col_fov = "fov", fov = 1,
    norm = FALSE, viridis = "D"
  ) +
    scale_fill_viridis_c(option = "D", limits = common_limits) +
    scale_color_viridis_c(option = "D", limits = common_limits) +
    theme_dark_clean +
    theme(legend.position = "none") +
    labs(x = NULL, y = NULL) # Explicitly remove labels just in case
  
  png(filename = paste0(out_dir, g, "_FOV1.png"), width = 7, height = 7, units = "in", res = 400)
  print(p1)
  dev.off()
  
  
  # D. PLOT FOV 3 (Clean, Dark, No Legend)
  
  p3 <- plot.gene.spatial(
    meta = meta,
    expression_matrix = as.matrix(seu@assays$RNA$data),
    genes = g,
    polygon_file = poly_data,
    x = "x_local_px", y = "y_local_px",
    col_fov = "fov", fov = 3,
    norm = FALSE, viridis = "D"
  ) +
    scale_fill_viridis_c(option = "D", limits = common_limits) +
    scale_color_viridis_c(option = "D", limits = common_limits) +
    theme_dark_clean +
    theme(legend.position = "none") +
    labs(x = NULL, y = NULL)
  
  png(filename = paste0(out_dir, g, "_FOV3.png"), width = 7, height = 7, units = "in", res = 400)
  print(p3)
  dev.off()
}

