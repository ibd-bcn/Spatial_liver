#
# Libraries ----------------------------------------------------------------
#
message('Loading libraries')

library(Seurat)
library(plyr)
library(ggplot2)
library(viridis)

#
# Extra functions -------------------------------------------------------------
#
message('Loading functions')
source('~/function.R')
#
# Data loading -----------------------------------------------------------------
#
message('Loading final subset objects')
plasmas <- readRDS("~/subsets/bcells/harmony3/plasmas_harmony.RDS")
hepatocytes <- readRDS("~/subsets/hepatocytes/harmony_new2/hepatocytes_harmony.RDS")
myeloids <- readRDS("~/subsets/myeloids/harmony4/myeloids.RDS")
fibroblasts <- readRDS("~/subsets/parenquimal/harmony2/fibroblasts_harmony.RDS")
tcells <- readRDS("~/subsets/tcells/harmony2/tcells_harmony.RDS")

#
# Manual annotation ------------------------------------------------------------
#
myeloid_annotation<- data.frame(
  Cluster = c(0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L),
  Myeloids = c("Monocytes", "KC1", "Neutrophils", "DCs CD1C", "KC2", "M2-LYVE1", "M1", "M1", "Mast cells", "Monocytes")
)


plasma_annotation <- data.frame(
  Cluster = c(0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L,8L,9L,10L),
  Plasmas = c("Plasma cells", "Memory B cells", "Plasma cells", "Naive B cells",
              "Plasma cells", "pDC", "Memory B cells", "Plasma cells",
              "Cycling B lineage cells", "Naive B cells", "Plasma cells")
)


tcell_annotation <- data.frame(
  Cluster = c(0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 11L, 12L),
  Tcells =  c("Trm cytotoxic T cells", "Naive T cells", "Effector helper T cells",
                        "Effector helper T cells", "Rb high", "Tem cytotoxic T cells",
                        "NKT cells", "Regulatory T cells", "Gamma-delta T cells",
                        "Naive T cells", "NK cells", "Effector helper T cells",
                        "Cycling T cells")

)

stroma_annotation <- data.frame(
  Cluster = c(0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 11L, 12L, 13L),
  Stromal = c("Endothelial cells 2", "Endothelial cells 1", "Smooth muscle cells",
              "Fibroblasts", "Endothelial cells 2", "Fibroblasts",
              "Endothelial cells 3", "Endothelial cells 4", "Endothelial cells 1",
              "Endothelial cells 1", "Endothelial cells 1", "Fibroblasts",
              "Schwann cells", "Fibroblasts")
)

hepatocytes_annotation <- data.frame(
  Cluster = c(0L, 1L, 2L, 3L, 4L, 5L, 6L),
  Epithelium = c("Hepatocyte 1", "Hepatocyte 2", "Hepatocyte 3",
                 "Cholangiocytes", "Hepatocyte 4", "Hepatocyte 5", "Hepatocyte 6")
)


hepatocytes$annotation <- mapvalues(x = hepatocytes$RNA_snn_res.0.1,
                                  from = hepatocytes_annotation$Cluster,
                                  to = hepatocytes_annotation$Epithelium)

tcells$annotation <- mapvalues(x = tcells$RNA_snn_res.0.5,
                                     from = tcell_annotation$Cluster,
                                     to = tcell_annotation$Tcells)

plasmas$annotation <- mapvalues(x = plasmas$RNA_snn_res.0.5,
                                      from = plasma_annotation$Cluster,
                                      to = plasma_annotation$Plasmas)


myeloids$annotation<- mapvalues(x = myeloids$RNA_snn_res.0.5,
                                       from = myeloid_annotation$Cluster,
                                       to = myeloid_annotation$Myeloids)


fibroblasts$annotation <- mapvalues(x = fibroblasts$RNA_snn_res.0.5,
                                     from = stroma_annotation$Cluster,
                                     to = stroma_annotation$Stromal)

dir.create("~/final_pieces")
saveRDS(fibroblasts, "~/final_pieces/parenquimal.RDS")
saveRDS(myeloids, "~/final_pieces/myeloids.RDS")
saveRDS(tcells, "~/final_pieces/tcells.RDS")
saveRDS(plasmas, "~/final_pieces/plasmas.RDS")
saveRDS(hepatocytes, "~/final_pieces/hepatocytes.RDS")
