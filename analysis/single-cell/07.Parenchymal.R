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
message('Loading parenchymal/endothelial subset')
parenquimal <- readRDS(file = '~/subsets/parenquimal.RDS')

# parenquimal


parenquimal
# 33906 features across 15388 samples within 1 assay

#  remove cells with genes from other subsets
#
FeaturePlot(parenquimal, features= c('CD3E', 'CD3G', 'MS4A1', 'DERL3', 'GNLY'), order = T)
FeaturePlot(parenquimal, features= c('IGHA1', 'IGHG1', 'EPCAM', 'JCHAIN', 'ALB'), order = T)
counts <- parenquimal@assays$RNA$counts
p <- grep("CD3E$|CD3D$|CD3G$|MS4A1$|^DERL3$|^EPCAM$|^JCHAIN$|^GNLY$",rownames(parenquimal))
pp <- which(Matrix::colSums(counts[p,])>0)
length(pp)
# 1247
xx <-setdiff(colnames(parenquimal), names(pp))
parenquimal <- subset(parenquimal,cells = xx)
parenquimal
# 33514 features across 15738 samples within 1 assay

# # remove genes from IGs

gg <- rownames(parenquimal)[c(grep("^IGH",rownames(parenquimal)),
                              grep("^IGK", rownames(parenquimal)),
                              grep("^IGL", rownames(parenquimal)))]
genes <- setdiff(rownames(parenquimal),gg)
parenquimal <- subset(parenquimal,features = genes)
parenquimal
#33637 features across 14141 samples within 1 assay

#  %MT filtering
#
VlnPlot(parenquimal, feature = 'percent.mt')
parenquimal <- parenquimal[,parenquimal$percent.mt < 25]
# 33249 features across 14259 samples within 1 assay
#
counts <- parenquimal@assays$RNA$counts
pp <- which(Matrix::rowSums(counts)==0)
length(pp)
# 4160
xx <-setdiff(rownames(parenquimal), names(pp))
parenquimal <- subset(parenquimal, features = xx)
parenquimal
#29477 features across 12945 samples within 1 assay

#
#  Normalization, scaling and UMAP generation -----------------------------------
#
parenquimal <- seurat_to_pca(parenquimal)


PCS <- select_pcs(parenquimal, 2)
PCS2 <- select_pcs(parenquimal, 1.6) # 95
ElbowPlot(parenquimal, ndims = 100) +
  geom_vline(xintercept = PCS) +
  geom_vline(xintercept = PCS2, colour="#BB0000")+
  geom_vline(xintercept = 35, colour="gray")+
  annotate(geom="text", x=PCS-5, y=3, label= paste("sdev > 2; PCs =", PCS),
           color="black")+
  annotate(geom="text", x=PCS2-5, y=4, label= paste("sdev > 1.6; PCs =", PCS2),
           color="#BB0000")+
  labs(title = paste0(''))

parenquimal<- FindNeighbors(parenquimal,  dims = 1:41, reduction = 'pca')
parenquimal<-RunUMAP(parenquimal, dims=1:41, reduction = 'pca')
#
DimPlot(parenquimal, group.by = 'sample')
#

parenquimal <- RunHarmony(parenquimal, group.by = c('sample', 'sequencing'), dims.use = 1:41)
ElbowPlot(parenquimal, ndims = 100, reduction = 'harmony') +
  geom_vline(xintercept = 35, linetype = 2) +
  labs(title = paste0('Harmony - 27PCS'))
parenquimal<- FindNeighbors(parenquimal, reduction = "harmony", dims = 1:30)
parenquimal<-RunUMAP(parenquimal, dims=1:30, reduction= "harmony")

DimPlot(parenquimal, group.by = c('sequencing'))

dir.create('~/subsets/parenquimal')
dir.create('~/subsets/parenquimal/harmony')
parenquimal <- resolutions(parenquimal, resolutions = c(0.3,0.5,0.7,0.9,1.1,1.3,1.5),
                           workingdir = '~/subsets/parenquimal/harmony',
                           title = 'parenquimal_harmony cosmx')
saveRDS(parenquimal, '~/subsets/parenquimal/harmony/parenquimal_harmony.RDS')


# curating parenquimal
parenquimal <- parenquimal[,parenquimal$RNA_snn_res.0.7 != 3]
parenquimal <- parenquimal[,parenquimal$RNA_snn_res.0.7 != 9]
parenquimal <- parenquimal[,parenquimal$RNA_snn_res.0.7 != 12]
parenquimal <- parenquimal[,parenquimal$RNA_snn_res.0.7 != 14]
parenquimal <- parenquimal[,parenquimal$RNA_snn_res.0.7 != 15]
parenquimal <- parenquimal[,parenquimal$RNA_snn_res.0.7 != 19]
parenquimal <- parenquimal[,parenquimal$RNA_snn_res.0.7 != 20]
parenquimal <- seurat_to_pca(parenquimal)


PCS <- select_pcs(parenquimal, 2)
PCS2 <- select_pcs(parenquimal, 1.6) # 95
ElbowPlot(parenquimal, ndims = 100) +
  geom_vline(xintercept = PCS) +
  geom_vline(xintercept = PCS2, colour="#BB0000")+
  geom_vline(xintercept = 35, colour="gray")+
  annotate(geom="text", x=PCS-5, y=3, label= paste("sdev > 2; PCs =", PCS),
           color="black")+
  annotate(geom="text", x=PCS2-5, y=4, label= paste("sdev > 1.6; PCs =", PCS2),
           color="#BB0000")+
  labs(title = paste0(''))

parenquimal<- FindNeighbors(parenquimal,  dims = 1:41, reduction = 'pca')
parenquimal<-RunUMAP(parenquimal, dims=1:41, reduction = 'pca')
#
DimPlot(parenquimal, group.by = 'sample')
#

parenquimal <- RunHarmony(parenquimal, group.by = c('sample', 'sequencing'), dims.use = 1:41)
ElbowPlot(parenquimal, ndims = 100, reduction = 'harmony') +
  geom_vline(xintercept = 35, linetype = 2) +
  labs(title = paste0('Harmony - 27PCS'))
parenquimal<- FindNeighbors(parenquimal, reduction = "harmony", dims = 1:30)
parenquimal<-RunUMAP(parenquimal, dims=1:30, reduction= "harmony")

DimPlot(parenquimal, group.by = c('sequencing'))

dir.create('~/subsets/parenquimal/harmony2')

parenquimal <- resolutions(parenquimal, resolutions = c(0.3,0.5,0.7,0.9,1.1,1.3,1.5),
                           workingdir = '~/subsets/parenquimal/harmony2',
                           title = 'parenquimal_harmony cosmx')
saveRDS(parenquimal, '~/subsets/parenquimal/harmony2/parenquimal_harmony.RDS')
