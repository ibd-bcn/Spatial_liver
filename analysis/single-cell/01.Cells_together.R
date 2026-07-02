#
# Libraries
#
message('Loading libraries')

library(Seurat)
library(plyr)
library(ggplot2)
library(DropletUtils)
library(SingleCellExperiment)
library(scater)
library(scran)
library(scDblFinder)
library(viridis)
library(readr)
library(fs)
library(dplyr)
#
# Extra functions
#
message('Loading functions')
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
source('~/function.R') # extra functions can be found at https://github.com/ibd-bcn/tofacitinib_ibd/blob/main/Analysis/extra_functions/functions_rnaseq.R


samples <- c(
             "L743", "L758", "L896", "L902", "L917", "L933", "46", "48", "49", "50", "51", "52", "CHB9", "CHB7", "FC6", "CHB3", "CHB2", "CHB8")

patient_id <- c(
                "HDV01", "HDV02", "HDV05", "HDV06", "HDV07", "HDV08",
                "HC1", "HC2", "HC3", "HC4", "HC5", "HC6" ,"HBV1", "HBV2", "HBV3","HBV4","HBV5", "HBV6")

sequencing <- c(
                "GEX 5'", "GEM-X 5'", "GEM-X 5'", "GEM-X 5'", "GEM-X 5'", "GEM-X 5'",
                "GEX 3' v2", "GEX 3'", "GEX 3'", "GEX 3'", "GEX 3'", "GEX 3'", "GEX 3'", "GEX 3'", "GEX 3'", "GEX 3'", "GEX 3'", "GEX 3'")

pathology <- c(
               "HDV", "HDV", "HDV", "HDV", "HDV", "HDV",
               "HC", "HC", "HC", "HC", "HC", "HC", "HBV", "HBV", "HBV", "HBV", "HBV", "HBV" )

hbsag <- c(
           "+", "+", "+", "+", "+", "+",
           "n.a.", "n.a.", "n.a.", "n.a.", "n.a.", "n.a.", "+", "+", "-", "+", "+", "+")

hbeag <- c(
           "-", "-", "-", "-", "+", "-",
           "n.a.", "n.a.", "n.a.", "n.a.", "n.a.", "n.a.", "-", "-","-","+","+","+")

# Create the data frame
metadata <- data.frame(
  Samples = samples,
  Patient_id = patient_id,
  Sequencing = sequencing,
  Pathology = pathology,
  HBsAg = hbsag,
  HBeAg = hbeag,
  stringsAsFactors = FALSE
)
metadata


# Load necessary libraries
library(Seurat)
library(dplyr)
library(scDblFinder)
library(fs)

# Base directory containing the raw data
base_dir <- "~/raw_data"

# Initialize a list to store Seurat objects
list_data <- list()


# obtain filtered matrices
find_filtered_matrix <- function(base_dir, sample_name) {
  all_dirs <- dir_ls(base_dir, recurse = TRUE, type = "directory")
  filtered_dir <- all_dirs[grepl(sample_name, all_dirs) & grepl("filtered_feature_bc_matrix$", all_dirs)]
  if (length(filtered_dir) > 0) {
    return(filtered_dir[1]) # Return the first match
  } else {
    return(NA) # Return NA if no match is found
  }
}


for (sample_name in metadata$Samples) {
  # Locate the folder
  filtered_dir <- find_filtered_matrix(base_dir, sample_name)
  if (is.na(filtered_dir)) {
    cat("No directory found for sample:", sample_name, "\n")
    next
  }
  cat("Processing sample:", sample_name, "Directory:", filtered_dir, "\n")

  # Load data
  sce2 <- tryCatch(Read10X(filtered_dir), error = function(e) {
    cat("Error loading data for sample:", sample_name, "\n")
    return(NULL)
  })
  if (is.null(sce2)) next
  if (ncol(sce2) == 0 || nrow(sce2) == 0) {
    cat("Empty data for sample:", sample_name, "\n")
    next
  }

  # Filter features
  keep_feature <- rowSums(sce2 > 0) > 0
  sce2 <- sce2[keep_feature, , drop = FALSE]
  cat("Filtered dimensions of sce2:", dim(sce2), "\n")

  if (ncol(sce2) < 2 || nrow(sce2) < 2) {
    cat("Insufficient data for sample:", sample_name, "\n")
    next
  }

  # Detect doublets
  sce2 <- scDblFinder(sce2, verbose = FALSE)

  # Create metadata
  current_metadata <- metadata %>% filter(Samples == sample_name)
  if (nrow(current_metadata) == 0) {
    cat("No metadata found for sample:", sample_name, "\n")
    next
  }
  m4 <- data.frame(
    'sample' = rep(sample_name, ncol(sce2)),
    'doublet' = sce2$scDblFinder.class,
    'patient_id' = current_metadata$Patient_id,
    'sequencing' = current_metadata$Sequencing,
    'pathology' = current_metadata$Pathology,
    'hbsag' = current_metadata$HBsAg,
    'hbeag' = current_metadata$HBeAg
  )
  rownames(m4) <- colnames(sce2)

  # Create Seurat object
  data <- tryCatch(CreateSeuratObject(
    counts = counts(sce2),
    min.features = 100,
    project = sample_name,
    assay = "RNA",
    meta.data = m4
  ), error = function(e) {
    cat("Error creating Seurat object for sample:", sample_name, "\nError details:", e$message, "\n")
    return(NULL)
  })
  if (is.null(data)) next

  # Add mitochondrial percentage
  data[["percent.mt"]] <- PercentageFeatureSet(object = data, pattern = "^MT-")

  # Store in list
  list_data[[sample_name]] <- data
}


list_data <- unlist(list_data)
seudata <- list_data[[1]]

#
# Merging raw data
#
seudata <- merge(seudata, list_data[2:length(list_data)])
seudata <- JoinLayers(seudata)

dir.create("~/data")
saveRDS(seudata, file = '~/data/seudata.RDS')

#
# Raw data exploration and filtering
#
VlnPlot(seudata, features = c('percent.mt', 'nFeature_RNA'), group.by = 'sample')
VlnPlot(seudata, features = c( 'nCount_RNA'), group.by = 'sample')

# Singlet

seudata <- seudata[,seudata$doublet == 'singlet']
counts <- seudata@assays$RNA$counts
pp <- which(Matrix::rowSums(counts)==0)
length(pp)
# 198

xx <-setdiff(rownames(seudata), names(pp))
seudata <- subset(seudata, features = xx)
seudata # 33906 features across 164977 samples within 1 assay


meta <- seudata@meta.data
meta$density <- get_density(meta$percent.mt, meta$nFeature_RNA, n = 100)
ggplot(meta) +
  geom_point(aes(percent.mt, nFeature_RNA, color = density)) +
  geom_hline(yintercept = 200)+
  scale_color_viridis() +
  theme_classic() +
  theme(plot.title = element_text(size = 25))+
  labs(title = '')
dev.off()

seudata_f <- seudata[, seudata$percent.mt < 30 & seudata$nFeature_RNA > 200]

# 33906 features across 110388 samples within 1 assay
# Saving raw and filtered data
#
message('Saving raw and filtered data')

saveRDS(seudata_f,'~/data/seudata_f.RDS')
