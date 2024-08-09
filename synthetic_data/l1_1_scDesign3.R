library(scDesign3)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(viridis)
library(zellkonverter)
library(reticulate)
theme_set(theme_bw())

h5ad_file <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/xenium/xenium_objs/anndata/chuvio_L1_1.h5ad"

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

l1_1 <- l1_1[1:10, ]

set.seed(123)
example_simu <- scdesign3(
  sce = l1_1,
  assay_use = "counts",
  celltype = "singler_annotation",
  pseudotime = NULL,
  spatial = c("spatial1", "spatial2"),
  other_covariates = NULL,
  mu_formula = "s(spatial1, spatial2, bs = 'gp', k= 400)",
  sigma_formula = "1",
  family_use = "nb",
  n_cores = 2,
  usebam = FALSE,
  corr_formula = "1",
  copula = "gaussian",
  DT = TRUE,
  pseudo_obs = FALSE,
  return_model = FALSE,
  nonzerovar = FALSE,
  parallelization = "pbmcapply"
)

simu_sce <- SingleCellExperiment(list(counts = example_simu$new_count), colData = example_simu$new_covariate)
logcounts(simu_sce) <- log1p(counts(simu_sce))

VISIUM_dat_test <- data.frame(t(log1p(counts(l1_1)))) %>% as_tibble() %>% dplyr::mutate(X = colData(l1_1)$spatial1, Y = colData(l1_1)$spatial2) %>% tidyr::pivot_longer(-c("X", "Y"), names_to = "Gene", values_to = "Expression") %>% dplyr::mutate(Method = "Reference")
VISIUM_dat_scDesign3 <- data.frame(t(log1p(counts(simu_sce)))) %>% as_tibble() %>% dplyr::mutate(X = colData(simu_sce)$spatial1, Y = colData(simu_sce)$spatial2) %>% tidyr::pivot_longer(-c("X", "Y"), names_to = "Gene", values_to = "Expression") %>% dplyr::mutate(Method = "scDesign3")
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
ggsave("larger_plot_test_vignette.png", plot, width = 15, height = 10, units = "in", dpi = 300)

# Display the plot
print(plot)