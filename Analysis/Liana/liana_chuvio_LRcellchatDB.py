# import liana
import liana as li
# needed for visualization and toy data
import scanpy as sc
import pandas as pd
import numpy as np
import os

path = "/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/astear/data/rcd_annot/chuvio_L1_1_RCTD_lvl3.h5ad"
chuvio_L1_1 = sc.read_h5ad(path)
path = "/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/astear/data/rcd_annot/chuvio_L1_2_RCTD_lvl3.h5ad"
chuvio_L1_2 = sc.read_h5ad(path)
path = "/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/astear/data/rcd_annot/chuvio_L2_1_RCTD_lvl3.h5ad"
chuvio_L2_1 = sc.read_h5ad(path)
path = "/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/astear/data/rcd_annot/chuvio_L3_1_RCTD_lvl3.h5ad"
chuvio_L3_1 = sc.read_h5ad(path)
path = "/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/astear/data/rcd_annot/chuvio_L4_1_RCTD_lvl3.h5ad"
chuvio_L4_1 = sc.read_h5ad(path)


li.method.rank_aggregate(chuvio_L1_2,
                           groupby='singler_annotation',
                           resource_name = 'cellchatdb',
                           expr_prop=0.1, # must be expressed in expr_prop fraction of cells
                           min_cells = 5,
                           n_perms = 1000,
                           use_raw = False, # run on log- and library-normalized counts
                           verbose = True,
                           inplace = True
                          )

chuvio_L1_2_liana_cellchatdb = chuvio_L1_2.uns['liana_res']
chuvio_L1_2_liana_cellchatdb.to_csv('chuvio_L1_2_liana_cellchatdb.csv', index=False)

li.method.rank_aggregate(chuvio_L2_1,
                           groupby='singler_annotation',
                           resource_name = 'cellchatdb',
                           expr_prop=0.1, # must be expressed in expr_prop fraction of cells
                           min_cells = 5,
                           n_perms = 1000,
                           use_raw = False, # run on log- and library-normalized counts
                           verbose = True,
                           inplace = True
                          )

chuvio_L2_1_liana_cellchatdb = chuvio_L2_1.uns['liana_res']
chuvio_L2_1_liana_cellchatdb.to_csv('chuvio_L2_1_liana_cellchatdb.csv', index=False)

li.method.rank_aggregate(chuvio_L3_1,
                           groupby='singler_annotation',
                           resource_name = 'cellchatdb',
                           expr_prop=0.1, # must be expressed in expr_prop fraction of cells
                           min_cells = 5,
                           n_perms = 1000,
                           use_raw = False, # run on log- and library-normalized counts
                           verbose = True,
                           inplace = True
                          )

chuvio_L3_1_liana_cellchatdb = chuvio_L3_1.uns['liana_res']
chuvio_L3_1_liana_cellchatdb.to_csv('chuvio_L3_1_liana_cellchatdb.csv', index=False)


