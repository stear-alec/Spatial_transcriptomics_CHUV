library(scDesign3)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(viridis)
library(zellkonverter)
theme_set(theme_bw())

# Replace "path/to/file.h5ad" with the actual path to your H5AD file
chuvio_L1_1 <- readH5AD("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/astear/Analysis/Giotto/data/chuvio_L1_1_RCTD_lvl3_processed.h5ad")
