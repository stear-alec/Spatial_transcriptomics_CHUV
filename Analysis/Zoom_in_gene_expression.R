#Library
library(SpatialExperiment)
library(ggspavis)
library(STexampleData)
library(zellkonverter)

#Import plots of chuvio and lung
orig_chuvio <- readH5AD("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/xenium/xenium_objs/anndata/chuvio_L1_1.h5ad")
orig_lung <- readH5AD("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/xenium/xenium_objs/anndata/lung_L1_1.h5ad")

# Zoom in on their interactions 



# Can we see significant interactions (found in other samples) on the original one ? 
plotSpots(orig_chuvio, x_coord = "x_centroid", y_coord = "y_centroid",
          annotate = "sum_VEGFA_VIM_CD44", in_tissue = NULL, y_reverse = TRUE, pt.size = 0.01)




