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
message('Loading B cell subset')
plasmas <- readRDS(file = '~/subsets/bcells.RDS')





FeaturePlot(plasmas, features= c('MS4A1', 'DERL3', 'EPCAM'), order = T)
FeaturePlot(plasmas, features= c('CD3E', 'CD3D', 'CD3G', 'EPCAM'), order = T)
#
# #
# # cell doublets from other subsets
# #

plasmas <- bcells
counts <- plasmas@assays$RNA$counts
p <- grep("CD3E$|CD3D$|CD3G$|EPCAM$",rownames(plasmas))
pp <- which(Matrix::colSums(counts[p,])>0)
length(pp)
# 311
xx <-setdiff(colnames(plasmas), names(pp))
#
plasmas <- subset(plasmas,cells = xx)
plasmas
# 33906 features across 6410 samples within 1 assay
#
# #
# # %MT filtering
# #
plasmas <- plasmas[,plasmas$percent.mt < 25]
plasmas
# 33906 features across 6320 samples within 1 assay
#
# #
# # Genes with no counts filtering
# #
counts <- plasmas@assays$RNA$counts
pp <- which(Matrix::rowSums(counts)==0)
length(pp) #6726
xx <-setdiff(rownames(plasmas), names(pp))
plasmas <- subset(plasmas, features = xx)
plasmas
# 27180 features across 6320 samples within 1 assay
#
# #
# # Normalization, scaling and UMAP generation -----------------------------------
# #
plasmas <- seurat_to_pca(plasmas)

#
PCS <- select_pcs(plasmas, 2)
PCS2 <- select_pcs(plasmas, 1.6) #14
ElbowPlot(plasmas, ndims = 100) +
  geom_vline(xintercept = 22, linetype = 2)+
  geom_vline(xintercept = PCS) +
  geom_vline(xintercept = PCS2, colour="#BB0000")+
  annotate(geom="text", x=PCS-5, y=3, label= paste("sdev > 2; PCs =", PCS),
           color="black")+
  annotate(geom="text", x=PCS2-5, y=4, label= paste("sdev > 1.6; PCs =", PCS2),
           color="#BB0000")+
  labs(title = paste0('92 single-cell samples from both IBD and HC plasmas'))
#
#
plasmas<- FindNeighbors(plasmas,  dims = 1:36, reduction = 'pca')
plasmas<-RunUMAP(plasmas, dims=1:36, reduction = 'pca')
#
DimPlot(plasmas, group.by = 'sample')
#
#
# #
# # Louvain clustering with batch correction -------------------------------------
# #
plasmas <- RunHarmony(plasmas, group.by = c('sample', 'sequencing'), dims.use = 1:36)
ElbowPlot(plasmas, ndims = 100, reduction = 'harmony') +
  geom_vline(xintercept =14, linetype = 2) +
  labs(title = paste0('Harmony - 45PCS'))
plasmas<- FindNeighbors(plasmas, reduction = "harmony", dims = 1:20)
plasmas<-RunUMAP(plasmas, dims=1:20, reduction= "harmony")
#
dir.create('~/subsets/bcells/harmony')

plasmas <- resolutions(plasmas, resolutions = c(0.5,0.7,0.9,1.1,1.3,1.5),
                       workingdir = '~/subsets/bcells/harmony',
                       title = 'plasmas harmony cosmx')
saveRDS(plasmas, file = '~/subsets/bcells/harmony/plasmas_harmony.RDS')


# curating plasmas
plasmas <- plasmas[,plasmas$RNA_snn_res.0.5 != 6]
plasmas <- plasmas[,plasmas$RNA_snn_res.0.5 != 14]

plasmas <- seurat_to_pca(plasmas)

#
PCS <- select_pcs(plasmas, 2)
PCS2 <- select_pcs(plasmas, 1.6) #14
ElbowPlot(plasmas, ndims = 100) +
  geom_vline(xintercept = 22, linetype = 2)+
  geom_vline(xintercept = PCS) +
  geom_vline(xintercept = PCS2, colour="#BB0000")+
  annotate(geom="text", x=PCS-5, y=3, label= paste("sdev > 2; PCs =", PCS),
           color="black")+
  annotate(geom="text", x=PCS2-5, y=4, label= paste("sdev > 1.6; PCs =", PCS2),
           color="#BB0000")+
  labs(title = paste0('92 single-cell samples from both IBD and HC plasmas'))
#
#
plasmas<- FindNeighbors(plasmas,  dims = 1:40, reduction = 'pca')
plasmas<-RunUMAP(plasmas, dims=1:40, reduction = 'pca')
#
DimPlot(plasmas, group.by = 'sample')
#
#
# #
# # Louvain clustering with batch correction -------------------------------------
# #
plasmas <- RunHarmony(plasmas, group.by = c('sample', 'sequencing'), dims.use = 1:40)
ElbowPlot(plasmas, ndims = 100, reduction = 'harmony') +
  geom_vline(xintercept =14, linetype = 2) +
  labs(title = paste0('Harmony - 45PCS'))
plasmas<- FindNeighbors(plasmas, reduction = "harmony", dims = 1:20)
plasmas<-RunUMAP(plasmas, dims=1:20, reduction= "harmony")
#
dir.create('~/subsets/bcells/harmony2')

plasmas <- resolutions(plasmas, resolutions = c(0.5,0.7,0.9,1.1,1.3,1.5),
                       workingdir = '~/subsets/bcells/harmony2',
                       title = 'plasmas harmony cosmx')
saveRDS(plasmas, file = '~/subsets/bcells/harmony2/plasmas_harmony.RDS')


# second round of curation

plasmas <- plasmas[,plasmas$RNA_snn_res.0.5 != 5]
plasmas <- plasmas[,plasmas$RNA_snn_res.0.5 != 9]
plasmas <- seurat_to_pca(plasmas)

#
PCS <- select_pcs(plasmas, 2)
PCS2 <- select_pcs(plasmas, 1.6) #14
ElbowPlot(plasmas, ndims = 100) +
  geom_vline(xintercept = 22, linetype = 2)+
  geom_vline(xintercept = PCS) +
  geom_vline(xintercept = PCS2, colour="#BB0000")+
  annotate(geom="text", x=PCS-5, y=3, label= paste("sdev > 2; PCs =", PCS),
           color="black")+
  annotate(geom="text", x=PCS2-5, y=4, label= paste("sdev > 1.6; PCs =", PCS2),
           color="#BB0000")
#
#
plasmas<- FindNeighbors(plasmas,  dims = 1:40, reduction = 'pca')
plasmas<-RunUMAP(plasmas, dims=1:40, reduction = 'pca')
#
DimPlot(plasmas, group.by = 'sample')
#
#
# #
# # Louvain clustering with batch correction -------------------------------------
# #
plasmas <- RunHarmony(plasmas, group.by = c('sample', 'sequencing'), dims.use = 1:40)
ElbowPlot(plasmas, ndims = 100, reduction = 'harmony') +
  geom_vline(xintercept =14, linetype = 2) +
  labs(title = paste0('Harmony - 45PCS'))
plasmas<- FindNeighbors(plasmas, reduction = "harmony", dims = 1:20)
plasmas<-RunUMAP(plasmas, dims=1:20, reduction= "harmony")
#
dir.create('~/subsets/bcells/harmony3')

plasmas <- resolutions(plasmas, resolutions = c(0.5,0.7,0.9,1.1,1.3,1.5),
                       workingdir = '/subsets/bcells/harmony3',
                       title = 'plasmas harmony cosmx')
saveRDS(plasmas, file = '~/subsets/bcells/harmony3/plasmas_harmony.RDS')
