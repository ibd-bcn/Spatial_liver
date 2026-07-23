# ==============================================================================
# Script: 3.cellinteraction.R
# Description: Prepare data and process SCOTIA cell-cell interactions
# ==============================================================================

# Load dependencies
suppressPackageStartupMessages({
  library(Seurat)
  library(readr)
  library(ggplot2)
  library(plyr)
  library(dplyr)
})

# ==============================================================================
# Configuration & Paths
# Update these paths according to your local environment
# ==============================================================================
IN_OBJ_PATH    <- "/path/to/output/Objects/seurats_all.RDS"
LR_DB_PATH     <- "/path/to/database/human_lr_pair.rds"
SCOTIA_DIR     <- "/path/to/output/SCOTIA/Files/"
SCOTIA_RES_DIR <- "/path/to/output/SCOTIA/Results/"
PY_SCRIPT_PATH <- "/path/to/scripts/scotia_run.py"

# Ensure output directories exist
invisible(lapply(c(SCOTIA_DIR, SCOTIA_RES_DIR), dir.create, recursive = TRUE, showWarnings = FALSE))

# ==============================================================================
# Prepare Seurat Object & Broad Annotations
# ==============================================================================
seu <- readRDS(IN_OBJ_PATH)
seu <- JoinLayers(seu)
seu <- SCTransform(seu)

# Direct mapping from refined annotation to broad subset
broad_celltype_map <- c(
  "Gamma-delta T cells"       = "T_cells",
  "NK cells"                  = "T_cells",
  "NKT cells"                 = "T_cells",
  "Regulatory T cells"        = "T_cells",
  "Tem/Trm cytotoxic T cells" = "T_cells",
  "Effector helper T cells"   = "T_cells",
  "Naive T cells"             = "T_cells",
  "Memory B cells"            = "B_cells",
  "Plasma cells"              = "B_cells",
  "Naive B cells"             = "B_cells",
  "Monocytes"                 = "Myeloid_cells",
  "M2-LYVE1"                  = "Myeloid_cells",
  "DCs CD1C"                  = "Myeloid_cells",
  "i-Mac"                        = "Myeloid_cells",
  "i-KC"                       = "Kupffer_cells",
  "h-KC"                       = "Kupffer_cells",
  "Endothelial cells 1"       = "Non_parenchymal",
  "Endothelial cells 2"       = "Non_parenchymal",
  "Endothelial cells 3"       = "Non_parenchymal",
  "Endothelial cells 4"       = "Non_parenchymal",
  "Fibroblasts"               = "Non_parenchymal",
  "Smooth muscle cells"       = "Non_parenchymal",
  "Cholangiocytes"            = "Non_parenchymal",
  "Hepatocyte 1"              = "Hepatocytes",
  "Hepatocyte 2"              = "Hepatocytes",
  "Hepatocyte 3"              = "Hepatocytes",
  "Hepatocyte 5"              = "Hepatocytes",
  "Hepatocyte 6"              = "Hepatocytes"
)

# Apply mapping
seu@meta.data$subset <- unname(broad_celltype_map[seu@meta.data$refined])
genes <- rownames(seu)

# Save main metadata
meta_seu <- as.data.frame(seu@meta.data)
write_csv(meta_seu, file.path(SCOTIA_DIR, "meta_seu.csv"))

# ==============================================================================
# Export Data per Patient for SCOTIA
# ==============================================================================
patients <- unique(seu$tissue)

for (patient in patients) {
  # Extract metadata
  cell_id <- seu@meta.data$cell_names[seu@meta.data$tissue == patient]
  meta <- meta_seu[cell_id, ]
  
  df_meta <- data.frame(
    cell_id     = cell_id,
    fov         = meta$fov,
    annotation  = meta$subset,
    x_positions = meta$CenterX_global_px,
    y_positions = meta$CenterY_global_px
  )
  
  write_csv(df_meta, file.path(SCOTIA_DIR, paste0(patient, "_meta_def.csv")))
  
  # Extract expression (SCT counts)
  data_ex <- as.data.frame(t(seu@assays$SCT$counts[, cell_id]))
  data_ex$cell_id <- rownames(data_ex)
  data_ex$fov <- df_meta$fov
  data_ex <- data_ex[, c("cell_id", "fov", colnames(data_ex)[!colnames(data_ex) %in% c("cell_id", "fov")])]
  
  write_csv(data_ex, file.path(SCOTIA_DIR, paste0(patient, "_exp_def.csv")))
}

# ==============================================================================
# Ligand-Receptor Pairs Formatting
# ==============================================================================
# Database source: https://github.com/ZJUFanLab/CellTalkDB
lr_pair <- readRDS(LR_DB_PATH)
lr_pair <- lr_pair[lr_pair$ligand_gene_symbol %in% genes & lr_pair$receptor_gene_symbol %in% genes, ]
lr_pair <- lr_pair[, c("ligand_gene_symbol", "receptor_gene_symbol")]
colnames(lr_pair) <- c("l_gene", "r_gene")

write_csv(lr_pair, file.path(SCOTIA_DIR, "lr_pair.csv"))

# ==============================================================================
# Run Python SCOTIA Script
# ==============================================================================
# Note: Adjust 'taskset' arguments based on your cluster/HPC configuration
command <- paste("taskset -c 0,20 python3", PY_SCRIPT_PATH)
system(command)

# ==============================================================================
# Aggregate SCOTIA Results
# ==============================================================================
meta_seu <- read_csv(file.path(SCOTIA_DIR, "meta_seu.csv"), show_col_types = FALSE)

# Store results in a list instead of iterative rbind for better memory efficiency
ot_results_list <- list()

for (patient in patients) {
  fovs <- unique(meta_seu$fov[meta_seu$tissue == patient])
  
  for (fov in fovs) {
    ot_file <- file.path(SCOTIA_DIR, paste0(patient, "_fov_", fov, ".ot.csv"))
    meta_file <- file.path(SCOTIA_DIR, paste0(patient, "_fov_", fov, ".csv"))
    
    if (file.exists(ot_file) && file.exists(meta_file)) {
      meta <- read_delim(meta_file, delim = "\t", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)
      ot   <- read_delim(ot_file, delim = "\t", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)
      
      # Map IDs and coordinates
      ot$id_source <- mapvalues(ot$source_cell_idx, from = meta$index, to = meta$cell_id, warn_missing = FALSE)
      ot$id_receptor <- mapvalues(ot$receptor_cell_idx, from = meta$index, to = meta$cell_id, warn_missing = FALSE)
      
      ot$refined_source <- mapvalues(ot$id_source, from = meta_seu$cell_names, to = meta_seu$refined, warn_missing = FALSE)
      ot$refined_receptor <- mapvalues(ot$id_receptor, from = meta_seu$cell_names, to = meta_seu$refined, warn_missing = FALSE)
      
      ot$x_source <- as.numeric(mapvalues(ot$id_source, from = meta$cell_id, to = meta$x_positions, warn_missing = FALSE))
      ot$y_source <- as.numeric(mapvalues(ot$id_source, from = meta$cell_id, to = meta$y_positions, warn_missing = FALSE))
      ot$x_receptor <- as.numeric(mapvalues(ot$id_receptor, from = meta$cell_id, to = meta$x_positions, warn_missing = FALSE))
      ot$y_receptor <- as.numeric(mapvalues(ot$id_receptor, from = meta$cell_id, to = meta$y_positions, warn_missing = FALSE))
      
      ot$fov <- as.numeric(mapvalues(ot$id_source, from = meta_seu$cell_names, to = meta_seu$fov, warn_missing = FALSE))
      ot$tissue <- mapvalues(ot$id_source, from = meta_seu$cell_names, to = meta_seu$tissue, warn_missing = FALSE)
      
      ot_results_list[[length(ot_results_list) + 1]] <- ot
    }
  }
}

# Bind all rows into one final dataframe
df_final <- bind_rows(ot_results_list)

write_csv(df_final, file.path(SCOTIA_RES_DIR, "all_int.csv"))