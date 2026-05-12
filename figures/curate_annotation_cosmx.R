# ==============================================================================
# Script: 0.curate_annotation_cosmx.R
# Description: Master metadata curation, FOV classification, and cell typing
# ==============================================================================

# Load necessary libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

# ==============================================================================
# Configuration & Paths
# Update these paths according to your local environment
# ==============================================================================
IN_SEURAT_PATH  <- "/path/to/seu.RDS"
OUT_SEURAT_PATH <- "/path/to/seurats_annotated.RDS"

# Ensure output directory exists
out_dir <- dirname(OUT_SEURAT_PATH)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# ==============================================================================
# 1. Annotation Dictionaries
# ==============================================================================
subset_type <- c(
  "hepatocytes"     = "Hepatocytes",
  "myeloids"        = "Myeloid_cells",
  "plasmas"         = "B_cell_lineage",
  "tcells"          = "T_cell_lineage",
  "non-parenquimal" = "Non_parenchymal",
  "IM-hepatocytes"  = "Im_hepatocytes"
)

grouped_anot <- c(
  "Non_parenchymal" = "Hepatocytes_Non_Parenchymal",
  "Hepatocytes"     = "Hepatocytes_Non_Parenchymal",
  "Myeloid_cells"   = "Immune_cells",
  "B_cell_lineage"  = "Immune_cells",
  "T_cell_lineage"  = "Immune_cells",
  "Im_hepatocytes"  = "Hepatocytes_Non_Parenchymal"
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

# ==============================================================================
# 2. Load Data
# ==============================================================================
message("Loading Seurat Object...")
seu <- readRDS(IN_SEURAT_PATH)

if(exists("JoinLayers", where="package:Seurat")) {
  seu <- JoinLayers(seu)
} else {
  seu <- Spatial_Package::JoinLayers(seu)
}

# ==============================================================================
# 3. Process Metadata
# ==============================================================================
message("Curating metadata fields...")

meta <- seu@meta.data

meta <- meta %>%
  mutate(
    # --- FOV Spatial Classifications ---
    red = (tissue == "Slide_2" & fov %in% c(21, 22, 23, 24)),
    
    orange = case_when(
      tissue == "Slide_1" & fov == 25                            ~ TRUE,
      tissue == "Slide_2" & fov %in% c(2, 3, 4, 6, 7, 10)        ~ TRUE,
      tissue == "Slide_3" & fov %in% c(6, 12, 14, 23)            ~ TRUE,
      TRUE                                                       ~ FALSE
    ),
    
    # --- Clinical / Antigen Levels ---
    antigen_level = case_when(
      tissue == "Slide_1" & fov %in% c(1:8, 25)                  ~ "neg",
      tissue == "Slide_3" & fov %in% c(1, 4, 11, 16, 24, 25)     ~ "neg",
      tissue == "Slide_2" & fov %in% c(13, 14, 15, 17, 19, 20)   ~ "S_low",
      tissue == "Slide_3" & fov == 13                            ~ "S_low",
      tissue == "Slide_2" & fov %in% c(16, 18)                   ~ "S_high",
      tissue == "Slide_3" & fov %in% c(2, 3, 5, 6, 12, 14, 23)   ~ "S_high",
      tissue == "Slide_3" & fov %in% c(7, 8, 9, 10, 15, 17:22)   ~ "S+/-",
      tissue == "Slide_2" & fov %in% c(1, 5, 6, 7, 9, 11, 12)    ~ "S&D_low",
      tissue == "Slide_2" & fov %in% c(2, 3, 4, 8, 10)           ~ "S&D",
      TRUE                                                       ~ NA_character_
    ),
    
    # --- Etiology ---
    etiology = case_when(
      tissue == "Slide_healthy"                                  ~ "HC",
      tissue == "Slide_2" & fov %in% c(13:21)                    ~ "HBV",
      tissue == "Slide_3" & fov %in% c(1:16, 23:25)              ~ "HBV",
      tissue == "Slide_2" & fov %in% c(1:12, 22:24)              ~ "HDV RNA+",
      tissue == "Slide_1" & fov %in% c(1:8, 25)                  ~ "HDV RNA-",
      tissue == "Slide_3" & fov %in% c(17:22)                    ~ "HDV RNA-",
      TRUE                                                       ~ NA_character_
    ),
    
    # --- Patient Type ---
    patient_type = case_when(
      tissue == "Slide_2" & fov %in% c(13:21)                    ~ "B07",
      tissue == "Slide_3" & fov %in% c(1:16, 23:25)              ~ "BH129",
      tissue == "Slide_2" & fov %in% c(1:5)                      ~ "D03",
      tissue == "Slide_2" & fov %in% c(6:12, 22:24)              ~ "D13",
      tissue == "Slide_1" & fov %in% c(1:8, 25)                  ~ "N10",
      tissue == "Slide_3" & fov %in% c(17:22)                    ~ "N02",
      TRUE                                                       ~ NA_character_
    ),
    
    # --- Antigen Etiology Combinations ---
    Ag_etiology = case_when(
      tissue == "Slide_1" & fov %in% c(1:8, 25)                  ~ "Neg_HDVn",
      tissue == "Slide_3" & fov %in% c(1, 4, 11, 16, 24, 25)     ~ "Neg_HBV",
      tissue == "Slide_2" & fov %in% c(13, 14, 15, 17, 19, 20)   ~ "S_low_HBV",
      tissue == "Slide_3" & fov == 13                            ~ "S_low_HBV",
      tissue == "Slide_2" & fov %in% c(16, 18)                   ~ "S_high_HBV",
      tissue == "Slide_3" & fov %in% c(2, 3, 5, 6, 12, 14, 23)   ~ "S_high_HBV",
      tissue == "Slide_3" & fov %in% c(7:10, 15)                 ~ "S+/-_HBV",
      tissue == "Slide_3" & fov %in% c(17:22)                    ~ "S+/-_HDVn",
      tissue == "Slide_2" & fov %in% c(1, 5, 6, 7, 9, 11, 12)    ~ "S&D_low_HDVp",
      tissue == "Slide_2" & fov %in% c(2, 3, 4, 8, 10)           ~ "S&D_HDVp",
      TRUE                                                       ~ NA_character_
    ),
    
    Ag_general_etiology = case_when(
      tissue == "Slide_1" & fov %in% c(1:8, 25)                  ~ "Neg_HDVn",
      tissue == "Slide_3" & fov %in% c(1, 4, 11, 16, 24, 25)     ~ "Neg_HBV",
      tissue == "Slide_2" & fov %in% c(13:20)                    ~ "S+_HBV",
      tissue == "Slide_3" & fov %in% c(2, 3, 5, 6, 12, 13, 14, 23) ~ "S+_HBV",
      tissue == "Slide_3" & fov %in% c(7:10, 15)                 ~ "S+/-_HBV",
      tissue == "Slide_3" & fov %in% c(17:22)                    ~ "S+/-_HDVn",
      tissue == "Slide_2" & fov %in% c(1:12)                     ~ "S&D_HDVp",
      TRUE                                                       ~ NA_character_
    ),
    
    Ag_general = case_when(
      tissue == "Slide_1" & fov %in% c(1:8, 25)                  ~ "Neg",
      tissue == "Slide_3" & fov %in% c(1, 4, 11, 16, 24, 25)     ~ "Neg",
      tissue == "Slide_2" & fov %in% c(13:20)                    ~ "S+",
      tissue == "Slide_3" & fov %in% c(2, 3, 5, 6, 12, 13, 14, 23) ~ "S+",
      tissue == "Slide_3" & fov %in% c(7:10, 15, 17:22)          ~ "S+/-",
      tissue == "Slide_2" & fov %in% c(1:12)                     ~ "S&D",
      TRUE                                                       ~ NA_character_
    ),
    
    # --- Map Cell Types ---
    subset       = unname(subset_type[subset]),
    refined      = unname(cell_types[refined]),
    new_anot     = unname(new_anot[refined]),
    grouped_anot = unname(grouped_anot[subset])
  )

# Assign curated metadata back to Seurat object
seu@meta.data <- meta

# ==============================================================================
# 4. Save Curated Object
# ==============================================================================
message(paste("Saving annotated Seurat object to:", OUT_SEURAT_PATH))
saveRDS(seu, OUT_SEURAT_PATH)

message("Annotation workflow complete.")