# Trying Giotto 
python_path = "/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/astear/pkgs/miniforge3/envs/giotto_env/bin/"
default_python_path = python_path
my_python_path = python_path

#installGiottoEnvironment(force_environment = TRUE, mini_install_path = "/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/astear/pkgs/miniforge3/")
# conda env create -n giotto_env -f ./python/configuration/genv_current.yml

#debugging
#rs.restartR()

#default_instrs <- createGiottoInstructions()
#default_python_path <- default_instrs$python_path

# Direct reticulate to use Python within the Giotto Environment
reticulate::use_python(default_python_path)
my_instructions = createGiottoInstructions(python_path = python_path)

setwd("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/astear/Analysis/Giotto")

#Start of normal use
library(Giotto)
library(GiottoData)
library(SpatialExperiment)
library(ggspavis)
library(STexampleData)
library(zellkonverter)
library(SingleCellExperiment)
library(data.table)

#Testing with h5ad , create directly the giotto object from this ? check the class of the resulting object? 
#orig_chuvio <- readH5AD("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/xenium/xenium_objs/anndata/chuvio_L1_1.h5ad")
#asGiotto(orig_chuvio, transfer_features = TRUE, verbose = TRUE)
#chuvio_L1_1 <- readH5AD("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/astear/Analysis/Giotto/data/chuvio_L1_1_RCTD_lvl3_processed.h5ad")


#rownames(orig_chuvio)

# trying to import a h5ad file
chuvio_L1_1 <- anndataToGiotto(anndata_path ="/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/astear/Analysis/Giotto/data/chuvio_L1_1_RCTD_lvl3_processed.h5ad",
                                python_path = my_python_path,
                               metadata_cols = c("total_counts", "pct_counts_mt"))

#Trying to make my own giotto object from scratch

### Extract the logcounts assay and save as data.table


path_to_locations = "/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/astear/Analysis/Giotto/data/spatial_chuvio_L1_1.txt"
cell_locations = data.table::fread(path_to_locations)

expression_matrix <- fread(file = "log_counts.csv", drop = 1)
cell_locations <- fread(file = "spatial_chuvo_L1_1.csv", drop = 1)
my_giotto_object = createGiottoObject(raw_exprs = expression_matrix,
                                      spatial_locs = cell_locations,
                                      instructions = my_instructions,)

#Following vignette 
my_giotto_object <- filterGiotto(gobject = seqfish_mini, 
                                 expression_threshold = 0.5, 
                                 gene_det_in_min_cells = 20, 
                                 min_det_genes_per_cell = 0)
my_giotto_object <- normalizeGiotto(gobject = my_giotto_object)
## create network (required for binSpect methods)
my_giotto_object = createSpatialNetwork(gobject = my_giotto_object, minimum_k = 2)


# binSpect kmeans method
km_spatialgenes = binSpect(my_giotto_object, bin_method = 'kmeans')

spatGenePlot(my_giotto_object, expression_values = 'scaled', 
             genes = km_spatialgenes[1:2]$genes, point_size = 3,
             point_shape = 'border', point_border_stroke = 0.1, cow_n_col = 2)

# binSpect rank method
rnk_spatialgenes = binSpect(my_giotto_object, bin_method = 'rank')

spatGenePlot(my_giotto_object, expression_values = 'scaled', 
             genes = rnk_spatialgenes[1:2]$genes, point_size = 3,
             point_shape = 'border', point_border_stroke = 0.1, cow_n_col = 2)

# silhouetteRank method
silh_spatialgenes = silhouetteRank(my_giotto_object)

spatGenePlot(my_giotto_object, expression_values = 'scaled', 
             genes = silh_spatialgenes[1:2]$genes,  point_size = 3,
             point_shape = 'border', point_border_stroke = 0.1, cow_n_col = 2)

# spatialDE method
spatDE_spatialgenes = spatialDE(my_giotto_object)
results = data.table::as.data.table(spatDE_spatialgenes$results)
setorder(results, -LLR)

spatGenePlot(my_giotto_object, expression_values = 'scaled', 
             genes = results$g[1:2],  point_size = 3,
             point_shape = 'border', point_border_stroke = 0.1, cow_n_col = 2)

# spark method
spark_spatialgenes = spark(my_giotto_object)
setorder(spark_spatialgenes, adjusted_pvalue, combined_pvalue)

spatGenePlot(my_giotto_object, expression_values = 'scaled', 
             genes = spark_spatialgenes[1:2]$genes,  point_size = 3,
             point_shape = 'border', point_border_stroke = 0.1, cow_n_col = 2)

# trendsceek method
trendsc_spatialgenes = trendSceek(my_giotto_object)
trendsc_spatialgenes = data.table::as.data.table(trendsc_spatialgenes)

spatGenePlot(my_giotto_object, expression_values = 'scaled', 
             genes = trendsc_spatialgenes[1:2]$gene,  point_size = 3,
             point_shape = 'border', point_border_stroke = 0.1, cow_n_col = 2)
