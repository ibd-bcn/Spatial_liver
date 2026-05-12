# ==============================================================================
# Script: 2.celltype.R
# Description: Cell typing and refined classification using InSituType
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
})

# ==============================================================================
# Configuration & Paths
# Update these paths according to your local environment
# ==============================================================================
# Input Data
QC_SEURAT_PATH  <- "path/to/output/Objects/qc_seurats_pass.RDS"
REF_TODAS_PATH  <- "path/to/reference_data/todas.RDS"

# Output Directories
OUT_OBJ_DIR     <- "path/to/output/Objects/"
OUT_MARKER_DIR  <- "path/to/output/Markers/"

# Ensure output directories exist
invisible(lapply(c(OUT_OBJ_DIR, OUT_MARKER_DIR), dir.create, recursive = TRUE, showWarnings = FALSE))

# Setup Colors
cols <- c(
  paletteer_d("ggsci::default_igv"),
  paletteer_d("ggsci::category20_d3"),
  paletteer_d("ggsci::default_ucscgb")
)

# ==============================================================================
# Global Coarse Classification
# ==============================================================================
seurats <- readRDS(QC_SEURAT_PATH)

# Load and prepare reference object
todas <- readRDS(REF_TODAS_PATH) 
todas@meta.data$subset <- ifelse(
  test = todas@meta.data$annotation %in% c("Hepatocyte 2", "Hepatocyte 4", "Hepatocyte 5", "Hepatocyte 6"),
  yes = "IM_hepatocytes",
  no = todas@meta.data$subset
)
todas@active.ident <- as.factor(todas$subset)
ref_data <- AverageExpression(object = todas)

# Prepare InSituType input
meta <- seurats@meta.data
counts_c <- t(as.matrix(seurats[["RNA"]]$counts))
negpb_c <- Matrix::rowMeans(t(as.matrix(seurats[["Negprob"]]$counts)))
ifdata <- as.matrix(meta[, c("Mean.PanCK", "Mean.CD45", "Mean.DAPI", "Mean.MembraneStain_B2M")])

cohort <- fastCohorting(mat = ifdata, gaussian_transform = TRUE)

# Annotate CELLS
sup <- insitutype(
  x = counts_c,
  neg = negpb_c,
  cohort = cohort,
  reference_profiles = as.matrix(ref_data$RNA),
  n_clusts = 0,
  assay_type = "RNA"
) 

# Map classifications back to Seurat metadata
clust_names <- as.data.frame(sup$clust)
clust_names$cell_names <- rownames(clust_names)
clust_names$prob <- sup$prob

seurats@meta.data$subset <- mapvalues(x = seurats@meta.data$cell_names, from = clust_names$cell_names, to = clust_names$`sup$clust`)
seurats@meta.data$subset_prob <- mapvalues(x = seurats@meta.data$cell_names, from = clust_names$cell_names, to = clust_names$prob)

# Threshold of probability >0.90
seurats.90 <- seurats[, seurats@meta.data$subset_prob > 0.90]
saveRDS(seurats.90, file.path(OUT_OBJ_DIR, "seu90.RDS"))

# Flightpath Plot
cols_main <- cols[seq_along(unique(sup$clust))]
names(cols_main) <- unique(sup$clust)
fp <- flightpath_plot(flightpath_result = NULL, insitutype_result = sup, col = cols_main[sup$clust])
print(fp)

# ==============================================================================
# REFINED CLASSIFICATION BLOCKS
# ==============================================================================
# We use the >0.90 thresholded object (seurats.90) for all refined classifications

# ------------------------------------------------------------------------------
# 1. Myeloids
# ------------------------------------------------------------------------------
todas_myeloids <- NormalizeData(todas[, todas@meta.data$subset == "myeloids"])
todas_myeloids@active.ident <- as.factor(todas_myeloids$annotation)
data_myeloids <- AverageExpression(object = todas_myeloids)

myeloids <- seurats.90[, seurats.90@meta.data$subset == "myeloids"]
ifdata_myeloids <- as.matrix(myeloids@meta.data[, c("Mean.PanCK", "Mean.CD45", "Mean.DAPI", "Mean.MembraneStain_B2M")])

sup_myeloids <- insitutype(
  x = t(myeloids[["RNA"]]$counts),
  neg = Matrix::rowMeans(t(myeloids[["Negprob"]]$counts)),
  cohort = fastCohorting(mat = ifdata_myeloids, gaussian_transform = TRUE),
  reference_profiles = as.matrix(data_myeloids$RNA),
  n_clusts = 0, assay_type = "RNA"
)

clust_names_my <- as.data.frame(sup_myeloids$clust)
myeloids@meta.data$refined <- mapvalues(x = rownames(myeloids@meta.data), from = rownames(clust_names_my), to = clust_names_my$`sup_myeloids$clust`)
myeloids@meta.data$refined_prob <- mapvalues(x = rownames(myeloids@meta.data), from = rownames(clust_names_my), to = sup_myeloids$prob)
saveRDS(myeloids, file.path(OUT_OBJ_DIR, "myeloids.RDS"))

# Markers
markers_my_ref <- FindAllMarkers(todas_myeloids, features = intersect(rownames(myeloids), rownames(todas_myeloids)), only.pos = TRUE)
myeloids <- NormalizeData(myeloids)
markers_my_cosmx <- FindAllMarkers(myeloids, group.by = "refined", only.pos = TRUE)
write.csv(markers_my_ref, file.path(OUT_MARKER_DIR, "myeloids.csv"))
write.csv(markers_my_cosmx, file.path(OUT_MARKER_DIR, "myeloids_cosmx.csv"))

# ------------------------------------------------------------------------------
# 2. T Cells
# ------------------------------------------------------------------------------
todas_tcells <- todas[, todas@meta.data$subset == "tcells" & !(todas@meta.data$annotation %in% c("Rb high", "Cycling T cells"))]
todas_tcells$annotation <- gsub("Tem cytotoxic T cells|Trm cytotoxic T cells", "Tem/Trm cytotoxic T cells", todas_tcells$annotation)
todas_tcells <- NormalizeData(todas_tcells)
todas_tcells@active.ident <- as.factor(todas_tcells$annotation)
data_tcells <- AverageExpression(object = todas_tcells)

tcells <- seurats.90[, seurats.90@meta.data$subset == "tcells"]
ifdata_tcells <- as.matrix(tcells@meta.data[, c("Mean.PanCK", "Mean.CD45", "Mean.DAPI", "Mean.MembraneStain_B2M")])

sup_tcells <- insitutype(
  x = t(tcells[["RNA"]]$counts),
  neg = Matrix::rowMeans(t(tcells[["Negprob"]]$counts)),
  cohort = fastCohorting(mat = ifdata_tcells, gaussian_transform = TRUE),
  reference_profiles = as.matrix(data_tcells$RNA),
  n_clusts = 0, assay_type = "RNA"
)

clust_names_tc <- as.data.frame(sup_tcells$clust)
tcells@meta.data$refined <- mapvalues(x = rownames(tcells@meta.data), from = rownames(clust_names_tc), to = clust_names_tc$`sup_tcells$clust`)
tcells@meta.data$refined_prob <- mapvalues(x = rownames(tcells@meta.data), from = rownames(clust_names_tc), to = sup_tcells$prob)
saveRDS(tcells, file.path(OUT_OBJ_DIR, "tcells.RDS"))

# Markers
markers_tc_ref <- FindAllMarkers(todas_tcells, features = intersect(rownames(tcells), rownames(todas_tcells)), only.pos = TRUE)
tcells <- NormalizeData(tcells)
markers_tc_cosmx <- FindAllMarkers(tcells, group.by = "refined", only.pos = TRUE)
write.csv(markers_tc_ref, file.path(OUT_MARKER_DIR, "tcells.csv"))
write.csv(markers_tc_cosmx, file.path(OUT_MARKER_DIR, "tcells_cosmx.csv"))

# ------------------------------------------------------------------------------
# 3. Plasmas
# ------------------------------------------------------------------------------
todas_plasmas <- NormalizeData(todas[, todas@meta.data$subset == "plasmas"])
todas_plasmas@active.ident <- as.factor(todas_plasmas$annotation)
data_plasmas <- AverageExpression(object = todas_plasmas)

plasmas <- seurats.90[, seurats.90@meta.data$subset == "plasmas"]
ifdata_plasmas <- as.matrix(plasmas@meta.data[, c("Mean.PanCK", "Mean.CD45", "Mean.DAPI", "Mean.MembraneStain_B2M")])

sup_plasmas <- insitutype(
  x = t(plasmas[["RNA"]]$counts),
  neg = Matrix::rowMeans(t(plasmas[["Negprob"]]$counts)),
  cohort = fastCohorting(mat = ifdata_plasmas, gaussian_transform = TRUE),
  reference_profiles = as.matrix(data_plasmas$RNA),
  n_clusts = 0, assay_type = "RNA"
)

clust_names_pl <- as.data.frame(sup_plasmas$clust)
plasmas@meta.data$refined <- mapvalues(x = rownames(plasmas@meta.data), from = rownames(clust_names_pl), to = clust_names_pl$`sup_plasmas$clust`)
plasmas@meta.data$refined_prob <- mapvalues(x = rownames(plasmas@meta.data), from = rownames(clust_names_pl), to = sup_plasmas$prob)
saveRDS(plasmas, file.path(OUT_OBJ_DIR, "plasmas.RDS"))

# Markers
markers_pl_ref <- FindAllMarkers(todas_plasmas, features = intersect(rownames(plasmas), rownames(todas_plasmas)), only.pos = TRUE)
plasmas <- NormalizeData(plasmas)
markers_pl_cosmx <- FindAllMarkers(plasmas, group.by = "refined", only.pos = TRUE)
write.csv(markers_pl_ref, file.path(OUT_MARKER_DIR, "plasmas.csv"))
write.csv(markers_pl_cosmx, file.path(OUT_MARKER_DIR, "plasmas_cosmx.csv"))

# ------------------------------------------------------------------------------
# 4. Fibroblasts (Non-parenchymal)
# ------------------------------------------------------------------------------
todas_fibroblasts <- NormalizeData(todas[, todas@meta.data$subset == "non_parenquimal" & todas@meta.data$annotation != "Schwann cells"])
todas_fibroblasts@active.ident <- as.factor(todas_fibroblasts$annotation)
data_fibroblasts <- AverageExpression(object = todas_fibroblasts)

fibroblasts <- seurats.90[, seurats.90@meta.data$subset == "non-parenquimal"]
ifdata_fibroblasts <- as.matrix(fibroblasts@meta.data[, c("Mean.PanCK", "Mean.CD45", "Mean.DAPI", "Mean.MembraneStain_B2M")])

sup_fibro <- insitutype(
  x = t(fibroblasts[["RNA"]]$counts),
  neg = Matrix::rowMeans(t(fibroblasts[["Negprob"]]$counts)),
  cohort = fastCohorting(mat = ifdata_fibroblasts, gaussian_transform = TRUE),
  reference_profiles = as.matrix(data_fibroblasts$RNA),
  n_clusts = 0, assay_type = "RNA"
)

clust_names_fb <- as.data.frame(sup_fibro$clust)
fibroblasts@meta.data$refined <- mapvalues(x = rownames(fibroblasts@meta.data), from = rownames(clust_names_fb), to = clust_names_fb$`sup_fibro$clust`)
fibroblasts@meta.data$refined_prob <- mapvalues(x = rownames(fibroblasts@meta.data), from = rownames(clust_names_fb), to = sup_fibro$prob)
saveRDS(fibroblasts, file.path(OUT_OBJ_DIR, "fibroblasts.RDS"))

# Markers
markers_fb_ref <- FindAllMarkers(todas_fibroblasts, features = intersect(rownames(fibroblasts), rownames(todas_fibroblasts)), only.pos = TRUE)
fibroblasts <- NormalizeData(fibroblasts)
markers_fb_cosmx <- FindAllMarkers(fibroblasts, group.by = "refined", only.pos = TRUE)
write.csv(markers_fb_ref, file.path(OUT_MARKER_DIR, "parench.csv"))
write.csv(markers_fb_cosmx, file.path(OUT_MARKER_DIR, "parench_cosmx.csv"))

# ------------------------------------------------------------------------------
# 5. Immune Hepatocytes
# ------------------------------------------------------------------------------
todas_im_hepatocytes <- todas[, todas@meta.data$subset == "IM_hepatocytes"]
todas_im_hepatocytes@meta.data$annotation <- gsub("Hepatocyte 4", "Hepatocyte 5", todas_im_hepatocytes@meta.data$annotation)
todas_im_hepatocytes <- NormalizeData(todas_im_hepatocytes)
todas_im_hepatocytes@active.ident <- as.factor(todas_im_hepatocytes$annotation)
data_im_hepatocytes <- AverageExpression(object = todas_im_hepatocytes)

im_hepatocytes <- seurats.90[, seurats.90@meta.data$subset == "IM-hepatocytes"]
ifdata_im_hepa <- as.matrix(im_hepatocytes@meta.data[, c("Mean.PanCK", "Mean.CD45", "Mean.DAPI", "Mean.MembraneStain_B2M")])

sup_im_hepa <- insitutype(
  x = t(im_hepatocytes[["RNA"]]$counts),
  neg = Matrix::rowMeans(t(im_hepatocytes[["Negprob"]]$counts)),
  cohort = fastCohorting(mat = ifdata_im_hepa, gaussian_transform = TRUE),
  reference_profiles = as.matrix(data_im_hepatocytes$RNA),
  n_clusts = 0, assay_type = "RNA"
)

clust_names_ih <- as.data.frame(sup_im_hepa$clust)
im_hepatocytes@meta.data$refined <- mapvalues(x = rownames(im_hepatocytes@meta.data), from = rownames(clust_names_ih), to = clust_names_ih$`sup_im_hepa$clust`)
im_hepatocytes@meta.data$refined_prob <- mapvalues(x = rownames(im_hepatocytes@meta.data), from = rownames(clust_names_ih), to = sup_im_hepa$prob)
saveRDS(im_hepatocytes, file.path(OUT_OBJ_DIR, "im_hepatocytes.RDS"))

# Markers
markers_ih_ref <- FindAllMarkers(todas_im_hepatocytes, features = intersect(rownames(im_hepatocytes), rownames(todas_im_hepatocytes)), only.pos = TRUE)
im_hepatocytes <- NormalizeData(im_hepatocytes)
markers_ih_cosmx <- FindAllMarkers(im_hepatocytes, group.by = "refined", only.pos = TRUE)
write.csv(markers_ih_ref, file.path(OUT_MARKER_DIR, "im_hepatocytes.csv"))
write.csv(markers_ih_cosmx, file.path(OUT_MARKER_DIR, "im_hepatocytes_cosmx.csv"))

# ------------------------------------------------------------------------------
# 6. Hepatocytes
# ------------------------------------------------------------------------------
todas_hepatocytes <- NormalizeData(todas[, todas@meta.data$subset == "hepatocytes"])
todas_hepatocytes@active.ident <- as.factor(todas_hepatocytes$annotation)
data_hepatocytes <- AverageExpression(object = todas_hepatocytes)

hepatocytes <- seurats.90[, seurats.90@meta.data$subset == "hepatocytes"]
ifdata_hepatocytes <- as.matrix(hepatocytes@meta.data[, c("Mean.PanCK", "Mean.CD45", "Mean.DAPI", "Mean.MembraneStain_B2M")])

sup_hepa <- insitutype(
  x = t(hepatocytes[["RNA"]]$counts),
  neg = Matrix::rowMeans(t(hepatocytes[["Negprob"]]$counts)),
  cohort = fastCohorting(mat = ifdata_hepatocytes, gaussian_transform = TRUE),
  reference_profiles = as.matrix(data_hepatocytes$RNA),
  n_clusts = 0, assay_type = "RNA"
)

clust_names_h <- as.data.frame(sup_hepa$clust)
hepatocytes@meta.data$refined <- mapvalues(x = rownames(hepatocytes@meta.data), from = rownames(clust_names_h), to = clust_names_h$`sup_hepa$clust`)
hepatocytes@meta.data$refined_prob <- mapvalues(x = rownames(hepatocytes@meta.data), from = rownames(clust_names_h), to = sup_hepa$prob)
saveRDS(hepatocytes, file.path(OUT_OBJ_DIR, "hepatocytes.RDS"))

# Markers
markers_h_ref <- FindAllMarkers(todas_hepatocytes, features = intersect(rownames(hepatocytes), rownames(todas_hepatocytes)), only.pos = TRUE)
hepatocytes <- NormalizeData(hepatocytes)
markers_h_cosmx <- FindAllMarkers(hepatocytes, group.by = "refined", only.pos = TRUE)
write.csv(markers_h_ref, file.path(OUT_MARKER_DIR, "hepatocytes.csv"))
write.csv(markers_h_cosmx, file.path(OUT_MARKER_DIR, "hepato_cosmx.csv"))

# ==============================================================================
# Merge All Refined Objects
# ==============================================================================
# We merge the existing objects from memory directly to save time
seurats_all <- reduce(list(tcells, plasmas, myeloids, fibroblasts, hepatocytes, im_hepatocytes), merge)

saveRDS(seurats_all, file.path(OUT_OBJ_DIR, "seurats_all.RDS"))
write.csv(seurats_all@meta.data, file.path(OUT_OBJ_DIR, "meta.csv"), row.names = FALSE)