# ==============================================================================
# Script: 1.qc.R
# Description: Spatial transcriptomics raw data processing and Quality Control
# ==============================================================================

# Load dependencies
suppressPackageStartupMessages({
  library(Seurat)
  library(readr)
  library(BiocParallel)
  library(InSituType)
  library(plyr)
  library(paletteer)
  library(harmony)
  library(stringr)
  library(ggplot2)
  library(outliers)
  library(purrr)
  library(data.table)
})

# ==============================================================================
# Configuration & Paths
# Update these paths according to your local environment
# ==============================================================================
RAW_DATA_DIR    <- "/path/to/raw_data/5_Raw_data/"
HEALTHY_REF_DIR <- "/path/to/reference_data/Healthy_sample/"
OUT_POLY_DIR    <- "/path/to/output/Polygons/"
OUT_MOL_DIR     <- "/path/to/output/Molecules/"
OUT_META_DIR    <- "/path/to/output/Metadata/"
OUT_OBJ_DIR     <- "/path/to/output/Objects/"

# Ensure output directories exist
invisible(lapply(c(OUT_POLY_DIR, OUT_MOL_DIR, OUT_META_DIR, OUT_OBJ_DIR), dir.create, recursive = TRUE, showWarnings = FALSE))

# ==============================================================================
# Read RAW FILES & Process
# ==============================================================================
bp <- MulticoreParam(workers = 3, progressbar = TRUE)
dirs <- list.dirs(path = RAW_DATA_DIR, recursive = FALSE)
shared.genes <- readRDS(file.path(HEALTHY_REF_DIR, "shared.genes.RDS"))

# Process each sample directory
seurats_r <- BiocParallel::bplapply(dirs, BPPARAM = bp, \(dir_path) {
  setwd(dir_path)
  tissue_name <- basename(dir_path)
  
  # Read Raw Data
  raw  <- read_csv(list.files(pattern = "exprMat_file.csv"), show_col_types = FALSE)
  meta <- read_csv(list.files(pattern = "metadata_file.csv"), show_col_types = FALSE)
  pols <- read_csv(list.files(pattern = "polygons.csv"), show_col_types = FALSE)
  mols <- read_csv(list.files(pattern = "tx_file.csv"), show_col_types = FALSE)
  
  # Remove empty cells
  raw  <- raw[raw$cell_ID != 0, ]
  pols <- pols[pols$cellID != 0, ]
  mols <- mols[mols$cell_ID != 0, ]
  
  # Format cell names
  pols$cell_names <- paste0(tissue_name, "_fov", pols$fov, "_", pols$cellID)
  mols$cell_names <- paste0(tissue_name, "_fov", mols$fov, "_", mols$cell_ID)
  meta$cell_names <- paste0(tissue_name, "_fov", meta$fov, "_", meta$cell_ID)
  
  # Add metadata info
  meta$tissue <- tissue_name
  
  # Separate counts and negative probes
  raw <- raw[, 3:ncol(raw)]
  neg <- raw[, grep(pattern = "NegPrb", x = colnames(raw))]
  count <- raw[, setdiff(colnames(raw), colnames(neg))]
  count <- count[, shared.genes]
  
  # Save processed dataframes
  write.csv(pols, file.path(OUT_POLY_DIR, paste0(tissue_name, ".csv")), row.names = FALSE)
  write.csv(mols, file.path(OUT_MOL_DIR, paste0(tissue_name, ".csv")), row.names = FALSE)
  write.csv(meta, file.path(OUT_META_DIR, paste0(tissue_name, ".csv")), row.names = FALSE)
  
  # Create Seurat object
  seu <- CreateSeuratObject(counts = t(count), meta.data = meta, assay = "RNA")
  
  # Add Negprob assay
  Negprob_assay <- CreateSeuratObject(counts = t(neg), assay = "Negprob")
  seu[["Negprob"]] <- CreateAssayObject(counts = Negprob_assay@assays$Negprob$counts)
  
  # Rename cells to standard format
  seu <- RenameCells(seu, new.names = seu@meta.data$cell_names)
  
  return(seu)
})

# Merge all sample Seurat objects
seurats <- reduce(seurats_r, merge)

# ==============================================================================
# Process Healthy Reference Sample
# ==============================================================================
meta_h <- as.data.frame(read_csv(file.path(HEALTHY_REF_DIR, "meta.csv"), show_col_types = FALSE))
rownames(meta_h) <- meta_h$cell_names

# Standardize coordinate metadata column names
colnames(meta_h) <- gsub("CenterX", "CenterX_global_px", colnames(meta_h))
colnames(meta_h) <- gsub("CenterY", "CenterY_global_px", colnames(meta_h))
colnames(meta_h) <- gsub("-", ".", colnames(meta_h))
colnames(meta_h) <- gsub("Mean.Membrane", "Mean.MembraneStain_B2M", colnames(meta_h))
colnames(meta_h) <- gsub("Max.Membrane", "Max.MembraneStain_B2M", colnames(meta_h))

# Filter metadata by FOVs
meta_h <- meta_h[meta_h$fov %in% c(2,3,4,7,26,28,34,45,46,48,78,96), ]

# Read Counts
counts_h <- fread(file.path(HEALTHY_REF_DIR, "counts.csv"), data.table = FALSE)
genes_h <- counts_h$target
rownames(counts_h) <- genes_h

common_cells <- intersect(colnames(counts_h), meta_h$cell_names)
counts_h <- counts_h[, common_cells]
meta_h <- meta_h[common_cells, ]
meta_h$tissue <- "Slide_healthy"

# Extract negative probes and bad codes
neg_h <- counts_h[grep(pattern = "NegPrb", x = rownames(counts_h)), ]
falsecode_h <- counts_h[grep(pattern = "FalseCode", x = rownames(counts_h)), ]
bad_cells <- c(rownames(neg_h), rownames(falsecode_h))

count_h <- counts_h[!rownames(counts_h) %in% bad_cells, ]
count_h <- count_h[shared.genes, ]

# Create healthy Seurat object
seu_h <- CreateSeuratObject(counts = count_h, meta.data = meta_h, assay = "RNA")
seu_h[["Negprob"]] <- CreateAssayObject(counts = neg_h)

# Merge main objects with healthy reference and join layers
seurats <- reduce(list(seurats, seu_h), merge)
seurats <- JoinLayers(seurats)

# ==============================================================================
# Quality Control (QC)
# ==============================================================================

# FLAG 1: Minimal counts per cell
# Exclude cells with fewer than 25 total RNA counts
seurats$flag1 <- seurats$nCount_RNA > 25

# FLAG 2: Proportion of negative counts
# Flag cells where > 5% of the counts per cell are negative probes
seurats$flag2 <- (seurats$nCount_Negprob / (seurats$nCount_RNA + seurats$nCount_Negprob)) < 0.05

# FLAG 3: Cell Complexity
# Total counts must exceed the number of detected genes in the cell (ratio >= 1)
seurats$flag3 <- (seurats$nCount_RNA / seurats$nFeature_RNA) >= 1

# FLAG 4: Area Outliers (Grubb's Test)
# Remove overly large or small cells based on cell area
data_vector <- seurats$Area

# Detect and remove large outliers
repeat {
  grubbs_test <- grubbs.test(data_vector, type = 10)
  if (grubbs_test$p.value > 0.01) break
  
  outliers <- regmatches(grubbs_test$alternative, gregexpr("\\d+", grubbs_test$alternative))
  data_vector <- data_vector[!(data_vector %in% as.numeric(unlist(outliers)))]
}

# Detect and remove small outliers
repeat {
  grubbs_test <- grubbs.test(data_vector, type = 10, opposite = TRUE)
  if (grubbs_test$p.value > 0.01) break
  
  outliers <- regmatches(grubbs_test$alternative, gregexpr("\\d+", grubbs_test$alternative))
  data_vector <- data_vector[!(data_vector %in% as.numeric(unlist(outliers)))]
}

seurats$flag4 <- (seurats$Area >= min(data_vector)) & (seurats$Area <= max(data_vector))

# Define final QC pass status
seurats$pass <- seurats$flag1 * seurats$flag2 * seurats$flag3 * seurats$flag4 

# ==============================================================================
# Save Outputs
# ==============================================================================

# Save complete Seurat object (Pre-QC filtering)
saveRDS(seurats, file.path(OUT_OBJ_DIR, "qc_seurats.RDS"))

# Save filtered Seurat object (QC passed only)
seurats_pass <- seurats[, seurats@meta.data$pass == 1]
saveRDS(seurats_pass, file.path(OUT_OBJ_DIR, "qc_seurats_pass.RDS"))