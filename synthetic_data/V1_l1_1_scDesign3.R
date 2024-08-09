


library(base)
#search()
.libPaths()
lib_path = "/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/astear/pkgs/miniforge3/envs/R_scdesign/lib/R/library"
library(scDesign3)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(viridis, lib.loc = lib_path)
library(zellkonverter, lib.loc = lib_path)
library(reticulate)
#setting up reticulate
use_condaenv(condaenv = "R_scdesign", conda = "auto", required = NULL)

theme_set(theme_bw())

# Path to your .h5ad file 
h5ad_file <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/astear/Programming/synthetic_data/subsets/L1_1_sub20K.h5ad"

# Read the .h5ad file and convert it to a SingleCellExperiment object
l1_1 <- zellkonverter::readH5AD(h5ad_file)

# View the SingleCellExperiment object
l1_1

#Change the name of the columns to get spatial1, 2
# Assuming example_sce is your SingleCellExperiment object
col_data <- colData(l1_1)

# Rename the columns
colnames(col_data)[colnames(col_data) == "x_centroid"] <- "spatial1"
colnames(col_data)[colnames(col_data) == "y_centroid"] <- "spatial2"

# Assign the modified colData back to the SingleCellExperiment object
colData(l1_1) <- col_data

# Verify the changes
colnames(colData(l1_1))

#Change the name for the cell type, 
set.seed(123)
example_simu <- scdesign3(
  sce = l1_1,
  assay_use = "counts",
  celltype = "singler_annotation",
  pseudotime = NULL,
  spatial = c("spatial1", "spatial2"),
  other_covariates = NULL,
  mu_formula = "s(spatial1, spatial2, bs = 'gp', k= 100)",
  sigma_formula = "1",
  family_use = "nb",
  n_cores = 2,
  usebam = TRUE,
  corr_formula = "1",
  copula = "gaussian",
  DT = TRUE,
  pseudo_obs = FALSE,
  return_model = FALSE,
  nonzerovar = FALSE,
  parallelization = "pbmcapply",
)



simu_sce <- SingleCellExperiment(list(counts = example_simu$new_count), colData = example_simu$new_covariate)
logcounts(simu_sce) <- log1p(counts(simu_sce))

# saving the synthesized dataset as .h5ad file
zellkonverter::writeH5AD(simu_sce, "V1_simu_sce.h5ad")

counts_dense <- as.matrix(t(log1p(counts(l1_1))))
VISIUM_dat_test <- data.frame(counts_dense) %>% as_tibble() %>% dplyr::mutate(X = colData(l1_1)$spatial1, Y = colData(l1_1)$spatial2) %>% tidyr::pivot_longer(-c("X", "Y"), names_to = "Gene", values_to = "Expression") %>% dplyr::mutate(Method = "Reference")

simu_sce_dense <- as.matrix(t(log1p(counts(simu_sce))))
VISIUM_dat_scDesign3 <- data.frame(simu_sce_dense) %>% as_tibble() %>% dplyr::mutate(X = colData(simu_sce)$spatial1, Y = colData(simu_sce)$spatial2) %>% tidyr::pivot_longer(-c("X", "Y"), names_to = "Gene", values_to = "Expression") %>% dplyr::mutate(Method = "scDesign3")

VISIUM_dat <- bind_rows(VISIUM_dat_test, VISIUM_dat_scDesign3) %>% dplyr::mutate(Method = factor(Method, levels = c("Reference", "scDesign3")))

# Your existing code here

plot <- VISIUM_dat %>% 
  filter(Gene %in% rownames(l1_1)[1:5]) %>% 
  ggplot(aes(x = X, y = Y, color = Expression)) + 
  geom_point(size = 0.5) + 
  scale_colour_gradientn(colors = viridis_pal(option = "magma")(10), limits=c(0, 8)) + 
  coord_fixed(ratio = 1) + 
  facet_grid(Method ~ Gene) + 
  theme_gray()

# Save plot in a bigger format
ggsave("V1_larger_plot_l1_1.png", plot, width = 15, height = 10, units = "in", dpi = 300)

# Display the plot
print(plot)