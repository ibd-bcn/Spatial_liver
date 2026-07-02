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
message('Loading T cell subset')
tcells <- readRDS(file = '~/subsets/tcells.RDS')

# Tcells
#
# remove cells with genes from other subsets
# #
counts <- tcells@assays$RNA$counts

p <- grep("MS4A1$|DERL3$|EPCAM$|^HBB$|^FCGR3B$|^S100A12$|^CD14$|^ALB$|^JCHAIN$|^HNF4A$|^DCN$|^TTR$",rownames(tcells))
pp <- which(Matrix::colSums(counts[p,])>0)
length(pp)
# 25243
xx <-setdiff(colnames(tcells), names(pp))
tcells <- subset(tcells,cells = xx)
tcells
# 33906 features across 11175 samples within 1 assay
#
#
# #
# # remove genes from IGs
# #
gg <- rownames(tcells)[c(grep("^IGH",rownames(tcells)),
                         grep("^IGK", rownames(tcells)),
                         grep("^IGL", rownames(tcells)))]
genes <- setdiff(rownames(tcells),gg)


tcells <- subset(tcells, features = genes)
tcells


# tcells
# 33637 features across 11175 samples within 1 assay

#
# #
# # Filter by %MT genes
# #
VlnPlot(tcells, features = 'percent.mt')
tcells <- tcells[,tcells$percent.mt < 25]
# 33249 features across 11462 samples within 1 assay
#
# #
# # genes with no expression out
# #
counts <- tcells@assays$RNA$counts
pp <- which(Matrix::rowSums(counts)==0)
length(pp) #7615
xx <-setdiff(rownames(tcells), names(pp))
#
tcells <- subset(tcells, features = xx)
tcells
# 24444 features across 11050 samples within 1 assay
#
#
# #
# # Normalization, scaling and UMAP generation -----------------------------------
# #
options(future.globals.maxSize = 8000 * 1024^2)
tcells <- seurat_to_pca(tcells)


PCS <- select_pcs(tcells, 2)
PCS2 <- select_pcs(tcells, 1.6) #22
ElbowPlot(tcells, ndims = 100) +
  geom_vline(xintercept = PCS) +
  geom_vline(xintercept = 23, linetype = 2) +
  geom_vline(xintercept = PCS2, colour="#BB0000")+
  annotate(geom="text", x=PCS-5, y=3, label= paste("sdev > 2; PCs =", PCS),
           color="black")+
  annotate(geom="text", x=PCS2-5, y=4, label= paste("sdev > 1.6; PCs =", PCS2),
           color="#BB0000")
tcells<- FindNeighbors(tcells,  dims = 1:27, reduction = 'pca')
tcells<-RunUMAP(tcells, dims=1:27, reduction = 'pca')

DimPlot(tcells, group.by = 'sequencing')
# #
# # Louvain clustering with batch correction -------------------------------------
# #
tcells <- RunHarmony(tcells, group.by = c('sample', 'sequencing'), dims.use = 1:27)
ElbowPlot(tcells, ndims = 100, reduction = 'harmony') +
  geom_vline(xintercept = 15, linetype = 2) +
  labs(title = paste0('Harmony - 23PCS'))
tcells<- FindNeighbors(tcells, reduction = "harmony", dims = 1:20)
tcells<-RunUMAP(tcells, dims=1:20, reduction= "harmony")

DimPlot(tcells, group.by = 'pathology')
#
#
FeaturePlot(tcells, features= c('CD3E','CD3D', 'GNLY', 'TTR'), order = T)


dir.create('~/subsets/tcells/')
dir.create('~/subsets/tcells/harmony')


tcells <- resolutions(tcells,resolutions = c(0.5,0.7,0.9,1.1,1.3,1.5),
                      workingdir = '~/subsets/tcells/harmony',
                      title = 'tcells_harmony cosmx')
saveRDS(tcells,'~/subsets/tcells/harmony/tcells_harmony.RDS')

# curating t cell subset

tcells <- tcells[,tcells$RNA_snn_res.0.7 != 5]
tcells <- tcells[,tcells$RNA_snn_res.0.7 != 15]
tcells <- seurat_to_pca(tcells)


PCS <- select_pcs(tcells, 2)
PCS2 <- select_pcs(tcells, 1.6) #22
ElbowPlot(tcells, ndims = 100) +
  geom_vline(xintercept = PCS) +
  geom_vline(xintercept = 23, linetype = 2) +
  geom_vline(xintercept = PCS2, colour="#BB0000")+
  annotate(geom="text", x=PCS-5, y=3, label= paste("sdev > 2; PCs =", PCS),
           color="black")+
  annotate(geom="text", x=PCS2-5, y=4, label= paste("sdev > 1.6; PCs =", PCS2),
           color="#BB0000")
tcells<- FindNeighbors(tcells,  dims = 1:28, reduction = 'pca')
tcells<-RunUMAP(tcells, dims=1:28, reduction = 'pca')

DimPlot(tcells, group.by = 'sequencing')
# #
# # Louvain clustering with batch correction -------------------------------------
# #
tcells <- RunHarmony(tcells, group.by = c('sample', 'sequencing'), dims.use = 1:28)
ElbowPlot(tcells, ndims = 100, reduction = 'harmony') +
  geom_vline(xintercept = 15, linetype = 2) +
  labs(title = paste0('Harmony - 23PCS'))
tcells<- FindNeighbors(tcells, reduction = "harmony", dims = 1:20)
tcells<-RunUMAP(tcells, dims=1:20, reduction= "harmony")

DimPlot(tcells, group.by = 'pathology')
#
#
FeaturePlot(tcells, features= c('CD3E','CD3D', 'GNLY', 'TTR'), order = T)

dir.create('~/subsets/tcells/harmony2')


tcells <- resolutions(tcells,resolutions = c(0.5,0.7,0.9,1.1,1.3,1.5),
                      workingdir = '~/subsets/tcells/harmony2',
                      title = 'tcells_harmony cosmx')
saveRDS(tcells,'~/subsets/tcells/harmony2/tcells_harmony.RDS')
