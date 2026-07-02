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
message('Loading hepatocyte subset')
hepatocytes <- readRDS(file = '~/subsets/hepatocytes.RDS')

# hepatocytes

hepatocytes

# 33906 features across 34082 samples within 1 assay

#  remove cells with genes from other subsets
#
FeaturePlot(hepatocytes, features= c('CD3E', 'CD3G', 'MS4A1', 'DERL3', 'GNLY'), order = T)
FeaturePlot(hepatocytes, features= c('IGHA1', 'IGHG1', 'EPCAM', 'JCHAIN', 'ALB'), order = T)
counts <- hepatocytes@assays$RNA$counts
p <- grep("CD3E$|CD3D$|CD3G$|MS4A1$|^DERL3$|^JCHAIN$|^GNLY$",rownames(hepatocytes))
pp <- which(Matrix::colSums(counts[p,])>0)
length(pp)
# 3005
xx <-setdiff(colnames(hepatocytes), names(pp))
hepatocytes <- subset(hepatocytes,cells = xx)
hepatocytes
# 33906 features across 31077 samples within 1 assay

# # remove genes from IGs

gg <- rownames(hepatocytes)[c(grep("^IGH",rownames(hepatocytes)),
                              grep("^IGK", rownames(hepatocytes)),
                              grep("^IGL", rownames(hepatocytes)))]
genes <- setdiff(rownames(hepatocytes),gg)
hepatocytes <- subset(hepatocytes,features = genes)
hepatocytes
#33637 features across 31077 samples within 1 assay

#  %MT filtering
#
VlnPlot(hepatocytes, feature = 'percent.mt')
hepatocytes <- hepatocytes[,hepatocytes$percent.mt < 30]
# 33637 features across 31077 samples within 1 assay
#
counts <- hepatocytes@assays$RNA$counts
pp <- which(Matrix::rowSums(counts)==0)
length(pp)
# 2888
xx <-setdiff(rownames(hepatocytes), names(pp))
hepatocytes <- subset(hepatocytes, features = xx)
hepatocytes
#30749 features across 31077 samples within 1 assay
hepatocytes <- seurat_to_pca(hepatocytes)


PCS <- select_pcs(hepatocytes, 2)
PCS2 <- select_pcs(hepatocytes, 1.6) # 26
ElbowPlot(hepatocytes, ndims = 100) +
  geom_vline(xintercept = PCS) +
  geom_vline(xintercept = PCS2, colour="#BB0000")+
  geom_vline(xintercept = 35, colour="gray")+
  annotate(geom="text", x=PCS-5, y=3, label= paste("sdev > 2; PCs =", PCS),
           color="black")+
  annotate(geom="text", x=PCS2-5, y=4, label= paste("sdev > 1.6; PCs =", PCS2),
           color="#BB0000")+
  labs(title = paste0(''))

hepatocytes<- FindNeighbors(hepatocytes,  dims = 1:26, reduction = 'pca')
hepatocytes<-RunUMAP(hepatocytes, dims=1:26, reduction = 'pca')
#
DimPlot(hepatocytes, group.by = 'sample')
#

hepatocytes <- RunHarmony(hepatocytes, group.by = c('sample', 'sequencing'), dims.use = 1:26)
ElbowPlot(hepatocytes, ndims = 100, reduction = 'harmony') +
  geom_vline(xintercept = 35, linetype = 2) +
  labs(title = paste0('Harmony - 27PCS'))
hepatocytes<- FindNeighbors(hepatocytes, reduction = "harmony", dims = 1:20)
hepatocytes<-RunUMAP(hepatocytes, dims=1:20, reduction= "harmony")

DimPlot(hepatocytes, group.by = c('sequencing'))

dir.create('~/subsets/hepatocytes/')
dir.create('~/subsets/hepatocytes/harmony')
hepatocytes <- resolutions(hepatocytes, resolutions = c(0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5),
                           workingdir = '~/subsets/hepatocytes/harmony',
                           title = 'hepatocytes cosmx harmony')
saveRDS(hepatocytes, '~/subsets/hepatocytes/harmony/hepatocytes_harmony.RDS')

# cleaning

hepatocytes <- hepatocytes[, hepatocytes$RNA_snn_res.0.1 != 10]
hepatocytes <- hepatocytes[, hepatocytes$RNA_snn_res.0.1 != 12]
hepatocytes <- hepatocytes[, hepatocytes$RNA_snn_res.0.1 != 14]
hepatocytes <- seurat_to_pca(hepatocytes)


PCS <- select_pcs(hepatocytes, 2)
PCS2 <- select_pcs(hepatocytes, 1.6) # 26
ElbowPlot(hepatocytes, ndims = 100) +
  geom_vline(xintercept = PCS) +
  geom_vline(xintercept = PCS2, colour="#BB0000")+
  geom_vline(xintercept = 35, colour="gray")+
  annotate(geom="text", x=PCS-5, y=3, label= paste("sdev > 2; PCs =", PCS),
           color="black")+
  annotate(geom="text", x=PCS2-5, y=4, label= paste("sdev > 1.6; PCs =", PCS2),
           color="#BB0000")+
  labs(title = paste0(''))

hepatocytes<- FindNeighbors(hepatocytes,  dims = 1:26, reduction = 'pca')
hepatocytes<-RunUMAP(hepatocytes, dims=1:26, reduction = 'pca')
#
DimPlot(hepatocytes, group.by = 'sample')
#

hepatocytes <- RunHarmony(hepatocytes, group.by = c('sample', 'sequencing'), dims.use = 1:26)
ElbowPlot(hepatocytes, ndims = 100, reduction = 'harmony') +
  geom_vline(xintercept = 35, linetype = 2) +
  labs(title = paste0('Harmony - 27PCS'))
hepatocytes<- FindNeighbors(hepatocytes, reduction = "harmony", dims = 1:20)
hepatocytes<-RunUMAP(hepatocytes, dims=1:20, reduction= "harmony")

DimPlot(hepatocytes, group.by = c('sequencing'))

dir.create('~/subsets/hepatocytes/harmony')

hepatocytes <- resolutions(hepatocytes, resolutions = c(0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5),
                           workingdir = '~/subsets/hepatocytes/harmony',
                           title = 'hepatocytes cosmx harmony')
saveRDS(hepatocytes, '~/subsets/hepatocytes/harmony/hepatocytes_harmony.RDS')
##adding cells from other subset that were hepatocytes
# cluster is cluster 3 from first 0.7 of parenquimal that were hepatocytes
hepatocytes <- merge(hepatocytes, cluster)
hepatocytes <- JoinLayers(hepatocytes)
hepatocytes <- seurat_to_pca(hepatocytes)


PCS <- select_pcs(hepatocytes, 2)
PCS2 <- select_pcs(hepatocytes, 1.6) # 26
ElbowPlot(hepatocytes, ndims = 100) +
  geom_vline(xintercept = PCS) +
  geom_vline(xintercept = PCS2, colour="#BB0000")+
  geom_vline(xintercept = 35, colour="gray")+
  annotate(geom="text", x=PCS-5, y=3, label= paste("sdev > 2; PCs =", PCS),
           color="black")+
  annotate(geom="text", x=PCS2-5, y=4, label= paste("sdev > 1.6; PCs =", PCS2),
           color="#BB0000")+
  labs(title = paste0(''))

hepatocytes<- FindNeighbors(hepatocytes,  dims = 1:26, reduction = 'pca')
hepatocytes<-RunUMAP(hepatocytes, dims=1:26, reduction = 'pca')
#
DimPlot(hepatocytes, group.by = 'sample')
#

hepatocytes <- RunHarmony(hepatocytes, group.by = c('sample', 'sequencing'), dims.use = 1:26)
ElbowPlot(hepatocytes, ndims = 100, reduction = 'harmony') +
  geom_vline(xintercept = 35, linetype = 2) +
  labs(title = paste0('Harmony - 27PCS'))
hepatocytes<- FindNeighbors(hepatocytes, reduction = "harmony", dims = 1:20)
hepatocytes<-RunUMAP(hepatocytes, dims=1:20, reduction= "harmony")

DimPlot(hepatocytes, group.by = c('sequencing'))

dir.create('~/subsets/hepatocytes/harmony_new')

hepatocytes <- resolutions(hepatocytes, resolutions = c(0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5),
                           workingdir = '~/subsets/hepatocytes/harmony_new',
                           title = 'hepatocytes cosmx harmony')
saveRDS(hepatocytes, '~/subsets/hepatocytes/harmony_new/hepatocytes_harmony.RDS')

# cleaning round 2

hepatocytes <- hepatocytes[,hepatocytes$RNA_snn_res.0.5 != 10]
hepatocytes <- hepatocytes[,hepatocytes$RNA_snn_res.0.5 != 12]
hepatocytes <- hepatocytes[,hepatocytes$RNA_snn_res.0.5 != 14]
hepatocytes <- seurat_to_pca(hepatocytes)


PCS <- select_pcs(hepatocytes, 2)
PCS2 <- select_pcs(hepatocytes, 1.6) # 26
ElbowPlot(hepatocytes, ndims = 100) +
  geom_vline(xintercept = PCS) +
  geom_vline(xintercept = PCS2, colour="#BB0000")+
  geom_vline(xintercept = 35, colour="gray")+
  annotate(geom="text", x=PCS-5, y=3, label= paste("sdev > 2; PCs =", PCS),
           color="black")+
  annotate(geom="text", x=PCS2-5, y=4, label= paste("sdev > 1.6; PCs =", PCS2),
           color="#BB0000")+
  labs(title = paste0(''))

hepatocytes<- FindNeighbors(hepatocytes,  dims = 1:26, reduction = 'pca')
hepatocytes<-RunUMAP(hepatocytes, dims=1:26, reduction = 'pca')
#
DimPlot(hepatocytes, group.by = 'sample')
#

hepatocytes <- RunHarmony(hepatocytes, group.by = c('sample', 'sequencing'), dims.use = 1:26)
ElbowPlot(hepatocytes, ndims = 100, reduction = 'harmony') +
  geom_vline(xintercept = 35, linetype = 2) +
  labs(title = paste0('Harmony - 27PCS'))
hepatocytes<- FindNeighbors(hepatocytes, reduction = "harmony", dims = 1:20)
hepatocytes<-RunUMAP(hepatocytes, dims=1:20, reduction= "harmony")

DimPlot(hepatocytes, group.by = c('sequencing'))

dir.create('~/subsets/hepatocytes/harmony_new2')

hepatocytes <- resolutions(hepatocytes, resolutions = c(0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5),
                           workingdir = '~/subsets/hepatocytes/harmony_new2',
                           title = 'hepatocytes cosmx harmony')
saveRDS(hepatocytes, '~/subsets/hepatocytes/harmony_new2/hepatocytes_harmony.RDS')
