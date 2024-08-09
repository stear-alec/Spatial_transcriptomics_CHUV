#Library
library(SpatialExperiment)
library(ggspavis)
library(STexampleData)
library(zellkonverter)
library(RColorBrewer)
#library(sceasy)
#library(SeuratDisk)

#1. Importing the .h5ad objects and setting them up as "SpatialExperiment" for plotting. 
##25 genes
fp_Original_20g <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/astear/Programming/synthetic_data/sub20genes/chuvio_L1_1_20genes.h5ad"
Original_20g <- readH5AD(fp_Original_20g)
fp_synthetic_20g <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/astear/Programming/synthetic_data/sub20genes/V1_simu_sce_20genes.h5ad"
synthetic_20g <- readH5AD(fp_synthetic_20g)

## 20K random cells
fp_Original_20k <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/astear/Programming/synthetic_data/sub20K/L1_1_sub20K.h5ad"
Original_20k <- readH5AD(fp_Original_20k)
fp_synthetic_20k <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/astear/Programming/synthetic_data/sub20K/V1_simu_sce_20K.h5ad"
synthetic_20k <- readH5AD(fp_synthetic_20k)

## Tu-Lb cells
fp_Original_20g <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/astear/Programming/synthetic_data/sub20genes/chuvio_L1_1_20genes.h5ad"
Original_20g <- readH5AD(fp_Original_20g)
fp_synthetic_20g <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/astear/Programming/synthetic_data/sub20genes/V1_simu_sce_20genes.h5ad"
synthetic_20g <- readH5AD(fp_synthetic_20g)



#NON NORMALISED - original
# Filter out rows with total gene expression less than 20
file = Original_20k

total_gene_expr <- colSums(assay(file, "counts"))
cols_to_keep <- total_gene_expr >= 20
sample_filtered <- file[, cols_to_keep]

# Plotting the Sum of gene expression
##colors <- brewer.pal(n = 1, name = "Blues")
sample_filtered$library_size <- colSums(assay(sample_filtered, "counts"))
#counts_matrix <- assay(sample_filtered, "logcounts")
plotSpots(sample_filtered, x_coord = "x_centroid", y_coord = "y_centroid", annotate = "library_size", 
          in_tissue = NULL, y_reverse = TRUE, pt.size = 0.01)

#NON NORMALISED - synthetic
# Filter out rows with total gene expression less than 20
file = synthetic_20k

total_gene_expr <- colSums(assay(file, "counts"))
cols_to_keep <- total_gene_expr >= 20
sample_filtered <- file[, cols_to_keep]

# Plotting the Sum of gene expression
##colors <- brewer.pal(n = 1, name = "Blues")
sample_filtered$library_size <- colSums(assay(sample_filtered, "counts"))
#counts_matrix <- assay(sample_filtered, "logcounts")
plotSpots(sample_filtered, x_coord = "spatial1", y_coord = "spatial2", annotate = "library_size", 
          in_tissue = NULL, y_reverse = TRUE, pt.size = 0.01)


#LOG-NORMALISED
# Filter out rows with total gene expression less than 20
total_gene_expr <- colSums(assay(file, "logcounts"))
cols_to_keep <- total_gene_expr >= log(20)
sample_filtered <- file[, cols_to_keep]

# Plotting the Sum of gene expression
##colors <- brewer.pal(n = 1, name = "Blues")
sample_filtered$library_size <- colSums(assay(sample_filtered, "logcounts"))
#counts_matrix <- assay(sample_filtered, "logcounts")
plotSpots(sample_filtered, x_coord = "x_centroid", y_coord = "y_centroid", annotate = "library_size", 
          in_tissue = NULL, y_reverse = TRUE, pt.size = 0.01)


