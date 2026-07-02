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
source('~/function.R') # https://github.com/ibd-bcn/tofacitinib_ibd/blob/main/Analysis/extra_functions/functions_rnaseq.R
#

#
# Data loading -----------------------------------------------------------------
#
message('Loading filtered merged object')
seudata_f <- readRDS('~/data/seudata_f.RDS')

#
# Data processing --------------------------------------------------------------
#
seudata_f <- seurat_to_pca(seudata_f)

PCS <- select_pcs(seudata_f, 2) # 27
PCS2 <- select_pcs(seudata_f, 1.6) # 41
p <- ElbowPlot(seudata_f, ndims = 100) +
  geom_vline(xintercept = PCS) +
  geom_vline(xintercept = PCS2, colour="#BB0000")+
  annotate(geom="text", x=PCS-5, y=3, label= paste("sdev > 2; PCs =", PCS),
           color="black")+
  annotate(geom="text", x=PCS2-5, y=4, label= paste("sdev > 1.6; PCs =", PCS2),
           color="#BB0000")+
  labs(title = paste0('altogether'))

seudata_f <- FindNeighbors(seudata_f,  dims = 1:42, reduction = 'pca')
seudata_f <- RunUMAP(seudata_f, dims=1:42, reduction = 'pca')

# plot the output UMAP
DimPlot(seudata_f, group.by = 'sequencing') + labs(title = 'altogether')

seudata_f <- resolutions(seudata_f, resolutions = c(0.3,0.5),
                         workingdir = '~/data',
                         title = 'wide_todas')

saveRDS(seudata_f,'~/data/seudata_f.RDS')
# 33514 features across 95982 samples within 1 assay


#
# Splitting the cells into subsets ---------------------------------------------
#
## getting subsets!


message('Cluster annotation and splitting')
subsets <- tibble::tribble(
  ~cluster,    ~subset,
  0L,   "tcells",
  1L,   "hepatocytes",
  2L,   "tcells",
  3L,   "hepatocytes",
  4L,   "endothelial",
  5L,   "myeloids",
  6L,   "hepatocytes",
  7L,   "endothelial",
  8L,   "endothelial",
  9L,   "myeloids",
  10L,  "bcells",
  11L,  "myeloids",
  12L,  "tcells",
  13L,  "hepatocytes",
  14L,  "hepatocytes",
  15L,  "eliminate",
  16L,  "bcells",
  17L,  "tcells",
  18L,  "hepatocytes",
  19L, "endothelial",
  20L, "hepatocytes",
  21L, "hepatocytes",
  22L, "hepatocytes",
  23L, "myeloids",
  24L, "eliminate",
  25L, "eliminate",
  26L, "endothelial",
  27L, "hepatocytes",
  28L, "hepatocytes",
  29L, "bcells",
  30L, "myeloids",
  31L, "hepatocytes"


)

tcells <- seudata_f[,seudata_f$RNA_snn_res.0.3 %in% subsets$cluster[subsets$subset == 'tcells']]
bcells <- seudata_f[,seudata_f$RNA_snn_res.0.3 %in% subsets$cluster[subsets$subset == 'bcells']]
myeloids <- seudata_f[,seudata_f$RNA_snn_res.0.3 %in% subsets$cluster[subsets$subset == 'myeloids']]
parenquimal <- seudata_f[,seudata_f$RNA_snn_res.0.3 %in% subsets$cluster[subsets$subset == 'endothelial']]
hepatocytes <- seudata_f[,seudata_f$RNA_snn_res.0.3 %in% subsets$cluster[subsets$subset == 'hepatocytes']]

dir.create("~/subsets/")
saveRDS(tcells, file = '~/subsets/tcells.RDS')
saveRDS(bcells, file = '~/subsets/bcells.RDS')
saveRDS(myeloids, file = '~/subsets/myeloids.RDS')
saveRDS(hepatocytes, file = '~/subsets/hepatocytes.RDS')
saveRDS(parenquimal, file = '~/subsets/parenquimal.RDS')
