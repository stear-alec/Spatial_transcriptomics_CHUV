# import liana
import liana as li
# needed for visualization and toy data
import scanpy as sc
import pandas as pd
import numpy as np
import os


file_name = "breast_B1_1"

path = f"/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/xenium/xenium_objs/anndata/{file_name}.h5ad"
file = sc.read_h5ad(path)
li.method.rank_aggregate(file,
                           groupby='rctd_first_type',
                           resource_name = 'cellchatdb',
                           expr_prop=0.1, # must be expressed in expr_prop fraction of cells
                           min_cells = 5,
                           n_perms = 1000,
                           use_raw = False, # run on log- and library-normalized counts
                           verbose = True,
                           inplace = True
                          )
file_cellchatdb = file.uns['liana_res']
file_cellchatdb.to_csv(f'{file_name}_liana_cellchatdb.csv', index=False)

file_name = "breast_B1_2"

path = f"/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/xenium/xenium_objs/anndata/{file_name}.h5ad"
file = sc.read_h5ad(path)
li.method.rank_aggregate(file,
                           groupby='rctd_first_type',
                           resource_name = 'cellchatdb',
                           expr_prop=0.1, # must be expressed in expr_prop fraction of cells
                           min_cells = 5,
                           n_perms = 1000,
                           use_raw = False, # run on log- and library-normalized counts
                           verbose = True,
                           inplace = True
                          )
file_cellchatdb = file.uns['liana_res']
file_cellchatdb.to_csv(f'{file_name}_liana_cellchatdb.csv', index=False)

file_name = "breast_B2_1"

path = f"/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/xenium/xenium_objs/anndata/{file_name}.h5ad"
file = sc.read_h5ad(path)
li.method.rank_aggregate(file,
                           groupby='rctd_first_type',
                           resource_name = 'cellchatdb',
                           expr_prop=0.1, # must be expressed in expr_prop fraction of cells
                           min_cells = 5,
                           n_perms = 1000,
                           use_raw = False, # run on log- and library-normalized counts
                           verbose = True,
                           inplace = True
                          )
file_cellchatdb = file.uns['liana_res']
file_cellchatdb.to_csv(f'{file_name}_liana_cellchatdb.csv', index=False)

file_name = "breast_B3_1"

path = f"/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/xenium/xenium_objs/anndata/{file_name}.h5ad"
file = sc.read_h5ad(path)
li.method.rank_aggregate(file,
                           groupby='rctd_first_type',
                           resource_name = 'cellchatdb',
                           expr_prop=0.1, # must be expressed in expr_prop fraction of cells
                           min_cells = 5,
                           n_perms = 1000,
                           use_raw = False, # run on log- and library-normalized counts
                           verbose = True,
                           inplace = True
                          )
file_cellchatdb = file.uns['liana_res']
file_cellchatdb.to_csv(f'{file_name}_liana_cellchatdb.csv', index=False)

file_name = "breast_B4_1"

path = f"/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/xenium/xenium_objs/anndata/{file_name}.h5ad"
file = sc.read_h5ad(path)
li.method.rank_aggregate(file,
                           groupby='rctd_first_type',
                           resource_name = 'cellchatdb',
                           expr_prop=0.1, # must be expressed in expr_prop fraction of cells
                           min_cells = 5,
                           n_perms = 1000,
                           use_raw = False, # run on log- and library-normalized counts
                           verbose = True,
                           inplace = True
                          )
file_cellchatdb = file.uns['liana_res']
file_cellchatdb.to_csv(f'{file_name}_liana_cellchatdb.csv', index=False)