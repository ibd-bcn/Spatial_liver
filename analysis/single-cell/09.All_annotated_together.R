#
# Libraries ----------------------------------------------------------------
#
message('Loading libraries')

library(Seurat)
library(plyr)
library(ggplot2)
library(viridis)
library(harmony)

#
# Extra functions -------------------------------------------------------------
#
message('Loading functions')
source('~/function.R')
#
# Data loading -----------------------------------------------------------------
#
message('Loading annotated subset objects')
plasmas <- readRDS('~/final_pieces/plasmas.RDS')
tcells <- readRDS('~/final_pieces/tcells.RDS')
hepatocytes <- readRDS('~/final_pieces/hepatocytes.RDS')
myeloids <- readRDS('~/final_pieces/myeloids.RDS')
fibroblasts <- readRDS('~/final_pieces/parenquimal.RDS')

# Together


plasmas$subset <- 'plasmas'
tcells$subset <- 'tcells'
hepatocytes$subset <- 'hepatocytes'
myeloids$subset <- 'myeloids'
fibroblasts$subset <- 'parenquimal'


#
# Data joining -----------------------------------------------------------------
#
todas <- merge(plasmas, tcells)
todas <- merge(todas, hepatocytes)
todas <- merge(todas, myeloids)
todas <- merge(todas, fibroblasts)
todas <- JoinLayers(todas)

#
# UMAP Re-generation
#
todas <- seurat_to_pca(todas)
PCS <- select_pcs(todas, 2)
PCS2 <- select_pcs(todas, 1.6)

ElbowPlot(todas, ndims = 100) +
  geom_vline(xintercept = PCS) +
  geom_vline(xintercept = PCS2, colour="#BB0000")+
  annotate(geom="text", x=PCS-5, y=3, label= paste("sdev > 2; PCs =", PCS),
           color="black")+
  annotate(geom="text", x=PCS2-5, y=4, label= paste("sdev > 1.6; PCs =", PCS2),
           color="#BB0000")+
  labs(title = paste0('merge all'))

todas <- FindNeighbors(todas, dims = 1:45)
todas <- RunUMAP(todas, dims = 1:45)

DimPlot(todas, group.by = "annotation")
#
# Harmony
#
library(harmony)
todas <- RunHarmony(todas, group.by = c('sample', 'sequencing'), dims.use = 1:45)
ElbowPlot(todas, ndims = 100, reduction = 'harmony') +
  geom_vline(xintercept = 35, linetype = 2) +
  geom_vline(xintercept = 25, linetype = 2) +
  labs(title = paste0('Harmony - 34PCS'))
todas <- FindNeighbors(todas, reduction = "harmony", dims = 1:30)
todas <- RunUMAP(todas, dims=1:30, reduction= "harmony")

DimPlot(todas, group.by = 'annotation', label=T)+
  labs(title='All  cells - Harmony')

#hnf4a, ttr, tf y apoa1

FeaturePlot(todas, c("HNF4A", "TTR", "TF", "APOA1"), order = T)
FeaturePlot(todas, c("CD3E", "CD8A"), order = T)

saveRDS(todas, "~/final_pieces/todas.RDS")
