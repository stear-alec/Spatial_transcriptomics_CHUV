import os
import gc
import ot
import pickle
import anndata
import scanpy as sc
import pandas as pd
import numpy as np
from scipy import sparse
from scipy.stats import spearmanr, pearsonr
from scipy.spatial import distance_matrix
import matplotlib.pyplot as plt

import commot as ct

path_L1_1 = "/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/astear/Analysis/TensorC2C/data/chuvio_L1_1_VISTA.h5ad"
L1_1 = sc.read_h5ad(path_L1_1)

path_L1_2 = "/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/astear/Analysis/TensorC2C/data/chuvio_L1_2_VISTA.h5ad"
L1_2 = sc.read_h5ad(path_L1_2)

path_chuvio_L4_1 = "/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/xenium/xenium_objs/anndata/chuvio_L4_1.h5ad"
L4_1 = sc.read_h5ad(path_chuvio_L4_1)

L1_1.obsm['spatial'] = L1_1.obs[['x_centroid', 'y_centroid']].to_numpy()
L1_2.obsm['spatial'] = L1_2.obs[['x_centroid', 'y_centroid']].to_numpy()
L4_1.obsm['spatial'] = L4_1.obs[['x_centroid', 'y_centroid']].to_numpy()

L1_1.raw = L1_1
L1_2.raw = L1_2
L4_1.raw = L4_1

L1_1_dis500 = L1_1.copy()
L1_2_dis500 = L1_2.copy()
L4_1_dis500 = L4_1.copy()

df_cellchat = ct.pp.ligand_receptor_database(species='human', signaling_type='Secreted Signaling', database='CellChat')
print(df_cellchat.shape)

L1_1_cellchat_filtered = ct.pp.filter_lr_database(df_cellchat, L1_1_dis500, min_cell_pct=0.00000001)
L1_2_cellchat_filtered = ct.pp.filter_lr_database(df_cellchat, L1_2_dis500, min_cell_pct=0.00000001)
L4_1_cellchat_filtered = ct.pp.filter_lr_database(df_cellchat, L4_1_dis500, min_cell_pct=0.00000001)


ct.tl.spatial_communication(L1_1_dis500,
    database_name='cellchat', df_ligrec=L1_1_cellchat_filtered, dis_thr=500, heteromeric=True, pathway_sum=True)

ct.tl.spatial_communication(L1_2_dis500,
    database_name='cellchat', df_ligrec=L1_2_cellchat_filtered, dis_thr=500, heteromeric=True, pathway_sum=True)

ct.tl.spatial_communication(L4_1_dis500,
    database_name='cellchat', df_ligrec=L4_1_cellchat_filtered, dis_thr=500, heteromeric=True, pathway_sum=True)

L1_1_dis500.write("./L1_1_dis500_commot.h5ad")
L1_2_dis500.write("./L1_2_dis500_commot.h5ad")
L4_1_dis500.write("./L4_1_dis500_commot.h5ad")
