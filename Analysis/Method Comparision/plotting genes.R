#Library
library(SpatialExperiment)
library(ggspavis)
library(STexampleData)
library(zellkonverter)
library(RColorBrewer)
#library(sceasy)
#library(SeuratDisk)

#1. Importing the .h5ad objects and setting them up as "SpatialExperiment" for plotting. 
#Set all the file paths
fp_chuvio_L1_1 <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/astear/data/rcd_annot/chuvio_L1_1_RCTD_lvl3.h5ad"
fp_chuvio_L1_2 <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/astear/data/rcd_annot/chuvio_L1_2_RCTD_lvl3.h5ad"
fp_chuvio_L2_1 <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/astear/data/rcd_annot/chuvio_L2_1_RCTD_lvl3.h5ad"
fp_chuvio_L3_1 <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/astear/data/rcd_annot/chuvio_L3_1_RCTD_lvl3.h5ad"
fp_chuvio_L4_1 <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/astear/data/rcd_annot/chuvio_L4_1_RCTD_lvl3.h5ad"


# Use the readH5AD function with the specified file path
chuvio_L1_1 <- readH5AD(fp_chuvio_L1_1)
chuvio_L1_2 <- readH5AD(fp_chuvio_L1_2)
chuvio_L2_1 <- readH5AD(fp_chuvio_L2_1)
chuvio_L3_1 <- readH5AD(fp_chuvio_L3_1)
chuvio_L4_1 <- readH5AD(fp_chuvio_L4_1)

#2. Plotting the Sum of gene expression
#colData(adata_Lung_sub_chuvio)$sum_CD80_CDH1_CTLA4_EGFR <- colSums(counts(adata_Lung_sub_chuvio))




# Getting a sum of counts of genes: Filter rows based on row names
sample <- chuvio_L1_1
sample <- sample2
#filtering
#counts_matrix <- assay(sample, "logcounts")
#rows_to_keep <- rowSums(counts_matrix >= 5) > 0
#sample_filtered <- sample[rows_to_keep, ]

rows_to_sum <- c("ERBB2", "EGFR")
filtered_counts <- counts(sample)[rows_to_sum, , drop = FALSE]
# Sum the counts
sum_counts <- colSums(filtered_counts)
# Add the sum to the colData
colData(sample_filtered)$sum_CCL5_CCR1 <- sum_counts

colors <- brewer.pal(n = 9, name = "Blues")
#not working -> Need a SpatialExperiment object not a SingleCellExperiment
plotSpots(sample_filtered, x_coord = "x_centroid", y_coord = "y_centroid",
          annotate = "sum_CCL5_CCR1", pal = colors, in_tissue = NULL, y_reverse = TRUE, pt.size = 0.01)


#testing
colnames(colData(adata_chuvio_sub_chromium))

#3. looking at the direction of cell signaling 
plotSpots(adata_chuvio_set_78, x_coord = "x_centroid", y_coord = "y_centroid",annotate = "sum_CD80_CDH1_CTLA4_EGFR", in_tissue = NULL, y_reverse = TRUE, pt.size = 0.01)

# Check the resulting SpatialExperiment object
print(spatial_exp)

#Import the anndata objects



#Code from Estella for plotting the cells
#gspavis::plotSpots(spe2_sub, annotate = x, in_tissue = NULL, y_reverse = TRUE, pt.size = 0.05))
# ggspavis::plotSpots(spe2_sub, annotate = "KRT7",
#                     title = paste0(slice_name, ": spatial distribution of log total gene expression"))



#Testing

sample2 <- SingleCellExperiment(assays = list(counts = counts_matrix))

# Transfer relevant metadata and colData from 'sample' to 'sample2'
metadata(sample2) <- metadata(sample)
colData(sample2) <- colData(sample)
rowData(sample2) <- rowData(sample)

#Chacking
counts_matrix <- assay(sample2, "counts")



