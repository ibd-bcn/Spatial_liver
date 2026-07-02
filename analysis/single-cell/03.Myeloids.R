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
message('Loading myeloid subset')
myeloids <- readRDS(file = '~/subsets/myeloids.RDS')

# subsets ----------------------------------------------------------------------------

#
# Data filtering  --------------------------------------------------------------
#
myeloids
# 33906 features across 13045 samples within 1 assay

#  remove cells with genes from other subsets
#


FeaturePlot(myeloids, features= c('CD3E', 'CD3G', 'MS4A1', 'DERL3', 'GNLY'), order = T)
FeaturePlot(myeloids, features= c('IGHA1', 'IGHG1', 'EPCAM', 'JCHAIN', 'ALB'), order = T)
counts <- myeloids@assays$RNA$counts
p <- grep("CD3E$|CD3D$|CD3G$|MS4A1$|^HNF4A$|^DCN$|^JCHAIN$|^GNLY$",rownames(myeloids))
pp <- which(Matrix::colSums(counts[p,])>0)
length(pp)
# 1338
xx <-setdiff(colnames(myeloids), names(pp))
myeloids <- subset(myeloids,cells = xx)
myeloids
# 33906 features across 11707 samples within 1 assay

# # remove genes from IGs

gg <- rownames(myeloids)[c(grep("^IGH",rownames(myeloids)),
                           grep("^IGK", rownames(myeloids)),
                           grep("^IGL", rownames(myeloids)))]
genes <- setdiff(rownames(myeloids),gg)


myeloids <- subset(myeloids, features = genes)

myeloids # 33637 features across 11707 samples within 1 assay


#  %MT filtering
#
VlnPlot(myeloids, feature = 'percent.mt')
myeloids <- myeloids[,myeloids$percent.mt < 25] #33637 features across 10889 samples within 1 assay

#
counts <- myeloids@assays$RNA$counts
pp <- which(Matrix::rowSums(counts)==0)
length(pp)
#  5551
xx <-setdiff(rownames(myeloids), names(pp))


myeloids <- subset(myeloids, features = xx)
myeloids
# 28086 features across 10889 samples within 1 assay

#
#  Normalization, scaling and UMAP generation -----------------------------------
#
options(future.globals.maxSize = 8000 * 1024^2)
myeloids <- seurat_to_pca(myeloids)


PCS <- select_pcs(myeloids, 2)
PCS2 <- select_pcs(myeloids, 1.6) # 95
ElbowPlot(myeloids, ndims = 100) +
  geom_vline(xintercept = PCS) +
  geom_vline(xintercept = PCS2, colour="#BB0000")+
  geom_vline(xintercept = 35, colour="gray")+
  annotate(geom="text", x=PCS-5, y=3, label= paste("sdev > 2; PCs =", PCS),
           color="black")+
  annotate(geom="text", x=PCS2-5, y=4, label= paste("sdev > 1.6; PCs =", PCS2),
           color="#BB0000")+
  labs(title = paste0(''))

myeloids<- FindNeighbors(myeloids,  dims = 1:30, reduction = 'pca')
myeloids<-RunUMAP(myeloids, dims=1:30, reduction = 'pca')
#
DimPlot(myeloids, group.by = 'sequencing')
DimPlot(myeloids, group.by = "pathology")

#



head(seurat_obj@meta.data)

# batch correction

library(harmony)
myeloids <- RunHarmony(myeloids, group.by.vars = c('sequencing', 'sample'), dims.use = 1:30)
ElbowPlot(myeloids, ndims = 100, reduction = 'harmony') +
  geom_vline(xintercept = 20, linetype = 2) +
  labs(title = paste0('Harmony - 20PCS'))
myeloids<- FindNeighbors(myeloids, reduction = "harmony", dims = 1:20)
myeloids<-RunUMAP(myeloids, dims=1:20, reduction= "harmony")

DimPlot(myeloids, group.by = c('sequencing'))

dir.create('~/subsets/myeloids/')
dir.create('~/subsets/myeloids/harmony')

myeloids <- resolutions(myeloids, resolutions = c(0.5,0.7,0.9,1.1,1.3,1.5),
                        workingdir = '~/subsets/myeloids/harmony',
                        title = 'myeloids_harmony cosmx')
saveRDS(myeloids, '~/subsets/harmony/myeloids.RDS')

# we eliminate several clusters that are not myeloid 

myeloids <- myeloids[,myeloids$RNA_snn_res.0.7 != 11]
myeloids <- myeloids[,myeloids$RNA_snn_res.0.7 != 12]
myeloids <- myeloids[,myeloids$RNA_snn_res.0.7 != 16]
myeloids <- myeloids[,myeloids$RNA_snn_res.0.7 != 17]
myeloids <- myeloids[,myeloids$RNA_snn_res.0.7 != 20]



# rerun normalization and batch correction

myeloids <- seurat_to_pca(myeloids)


PCS <- select_pcs(myeloids, 2)
PCS2 <- select_pcs(myeloids, 1.6) # 95
ElbowPlot(myeloids, ndims = 100) +
  geom_vline(xintercept = PCS) +
  geom_vline(xintercept = PCS2, colour="#BB0000")+
  geom_vline(xintercept = 35, colour="gray")+
  annotate(geom="text", x=PCS-5, y=3, label= paste("sdev > 2; PCs =", PCS),
           color="black")+
  annotate(geom="text", x=PCS2-5, y=4, label= paste("sdev > 1.6; PCs =", PCS2),
           color="#BB0000")+
  labs(title = paste0(''))

myeloids<- FindNeighbors(myeloids,  dims = 1:30, reduction = 'pca')
myeloids<-RunUMAP(myeloids, dims=1:30, reduction = 'pca')
#
DimPlot(myeloids, group.by = 'sequencing')
DimPlot(myeloids, group.by = "pathology")

#



library(harmony)
myeloids <- RunHarmony(myeloids, group.by.vars = c('sequencing', 'sample'), dims.use = 1:30)
ElbowPlot(myeloids, ndims = 100, reduction = 'harmony') +
  geom_vline(xintercept = 20, linetype = 2) +
  labs(title = paste0('Harmony - 20PCS'))
myeloids<- FindNeighbors(myeloids, reduction = "harmony", dims = 1:20)
myeloids<-RunUMAP(myeloids, dims=1:20, reduction= "harmony")

DimPlot(myeloids, group.by = c('sequencing'))

dir.create('~/subsets/myeloids/harmony2')


myeloids <- resolutions(myeloids, resolutions = c(0.5,0.7,0.9,1.1,1.3,1.5),
                        workingdir = '~/subsets/myeloids/harmony2',
                        title = 'myeloids_harmony cosmx')
saveRDS(myeloids, '~/myeloids/harmony2/myeloids.RDS')

# second run of curating the myeloid compartment at a higher resolution 

myeloids <- myeloids[,myeloids$RNA_snn_res.1.5 != 2]
myeloids <- myeloids[,myeloids$RNA_snn_res.1.5 != 3]
myeloids <- myeloids[,myeloids$RNA_snn_res.1.5 != 11]
myeloids <- myeloids[,myeloids$RNA_snn_res.1.5 != 15]
myeloids <- myeloids[,myeloids$RNA_snn_res.1.5 != 17]
myeloids <- myeloids[,myeloids$RNA_snn_res.1.5 != 18]
myeloids <- seurat_to_pca(myeloids)


PCS <- select_pcs(myeloids, 2)
PCS2 <- select_pcs(myeloids, 1.6) # 95
ElbowPlot(myeloids, ndims = 100) +
  geom_vline(xintercept = PCS) +
  geom_vline(xintercept = PCS2, colour="#BB0000")+
  geom_vline(xintercept = 35, colour="gray")+
  annotate(geom="text", x=PCS-5, y=3, label= paste("sdev > 2; PCs =", PCS),
           color="black")+
  annotate(geom="text", x=PCS2-5, y=4, label= paste("sdev > 1.6; PCs =", PCS2),
           color="#BB0000")+
  labs(title = paste0(''))

myeloids<- FindNeighbors(myeloids,  dims = 1:30, reduction = 'pca')
myeloids<-RunUMAP(myeloids, dims=1:30, reduction = 'pca')
#
DimPlot(myeloids, group.by = 'sequencing')
DimPlot(myeloids, group.by = "pathology")

#


library(harmony)
myeloids <- RunHarmony(myeloids, group.by.vars = c('sequencing', 'sample'), dims.use = 1:30)
ElbowPlot(myeloids, ndims = 100, reduction = 'harmony') +
  geom_vline(xintercept = 20, linetype = 2) +
  labs(title = paste0('Harmony - 20PCS'))
myeloids<- FindNeighbors(myeloids, reduction = "harmony", dims = 1:20)
myeloids<-RunUMAP(myeloids, dims=1:20, reduction= "harmony")

DimPlot(myeloids, group.by = c('sequencing'))

dir.create('~/subsets/myeloids/harmony3')


myeloids <- resolutions(myeloids, resolutions = c(0.5,0.7,0.9,1.1,1.3,1.5),
                        workingdir = '~/subsets/myeloids/harmony3',
                        title = 'myeloids_harmony cosmx')
saveRDS(myeloids, '~/subsets/myeloids/harmony3/myeloids.RDS')


# Third run of curating the myeloid compartment 

myeloids <- readRDS("~/subsets/myeloids/harmony3/myeloids.RDS")

myeloids <- myeloids[,myeloids$RNA_snn_res.0.7 != 10]
myeloids <- myeloids[,myeloids$RNA_snn_res.0.7 != 13]

myeloids <- seurat_to_pca(myeloids)


PCS <- select_pcs(myeloids, 2)
PCS2 <- select_pcs(myeloids, 1.6) # 95
ElbowPlot(myeloids, ndims = 100) +
  geom_vline(xintercept = PCS) +
  geom_vline(xintercept = PCS2, colour="#BB0000")+
  geom_vline(xintercept = 35, colour="gray")+
  annotate(geom="text", x=PCS-5, y=3, label= paste("sdev > 2; PCs =", PCS),
           color="black")+
  annotate(geom="text", x=PCS2-5, y=4, label= paste("sdev > 1.6; PCs =", PCS2),
           color="#BB0000")+
  labs(title = paste0(''))

myeloids<- FindNeighbors(myeloids,  dims = 1:30, reduction = 'pca')
myeloids<-RunUMAP(myeloids, dims=1:30, reduction = 'pca')
#
DimPlot(myeloids, group.by = 'sequencing')
DimPlot(myeloids, group.by = "pathology")

#

library(harmony)
myeloids <- RunHarmony(myeloids, group.by.vars = c('sequencing', 'sample'), dims.use = 1:30)
ElbowPlot(myeloids, ndims = 100, reduction = 'harmony') +
  geom_vline(xintercept = 20, linetype = 2) +
  labs(title = paste0('Harmony - 20PCS'))
myeloids<- FindNeighbors(myeloids, reduction = "harmony", dims = 1:20)
myeloids<-RunUMAP(myeloids, dims=1:20, reduction= "harmony")

DimPlot(myeloids, group.by = c('sequencing'))

dir.create('~/subsets/myeloids/harmony4')


myeloids <- resolutions(myeloids, resolutions = c(0.5,0.7,0.9,1.1,1.3,1.5),
                        workingdir = '~/subsets/myeloids/harmony4',
                        title = 'myeloids_harmony cosmx')
saveRDS(myeloids, '~/subsets/myeloids/harmony4/myeloids.RDS')
