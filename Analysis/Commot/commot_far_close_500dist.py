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

#Import the datasets


# importing datasets
path = '/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/astear/Liana+/rank_aggregate/Liana_RA_subsetTULB/anndata_subsetTULB/'

TU01plusLB_chuvio_L1_1 = sc.read_h5ad(f'{path}TU01plusLB_chuvio_L1_1.h5ad')
TU01plusLB_chuvio_L1_2 = sc.read_h5ad(f'{path}TU01plusLB_chuvio_L1_2.h5ad')
TU01plusLB_chuvio_L2_1 = sc.read_h5ad(f'{path}TU01plusLB_chuvio_L2_1.h5ad')
TU01plusLB_chuvio_L3_1 = sc.read_h5ad(f'{path}TU01plusLB_chuvio_L3_1.h5ad')
TU01plusLB_chuvio_L4_1 = sc.read_h5ad(f'{path}TU01plusLB_chuvio_L4_1.h5ad')

TU_MINUS01plusLB_chuvio_L1_1 = sc.read_h5ad(f'{path}TU_MINUS01plusLB_chuvio_L1_1.h5ad')
TU_MINUS01plusLB_chuvio_L1_2 = sc.read_h5ad(f'{path}TU_MINUS01plusLB_chuvio_L1_2.h5ad')
TU_MINUS01plusLB_chuvio_L2_1 = sc.read_h5ad(f'{path}TU_MINUS01plusLB_chuvio_L2_1.h5ad')
TU_MINUS01plusLB_chuvio_L3_1 = sc.read_h5ad(f'{path}TU_MINUS01plusLB_chuvio_L3_1.h5ad')
TU_MINUS01plusLB_chuvio_L4_1 = sc.read_h5ad(f'{path}TU_MINUS01plusLB_chuvio_L4_1.h5ad')

#Updating the spatial coordinates

TU01plusLB_chuvio_L1_1.obsm['spatial'] = TU01plusLB_chuvio_L1_1.obs[['x_centroid', 'y_centroid']].to_numpy()
TU01plusLB_chuvio_L1_2.obsm['spatial'] = TU01plusLB_chuvio_L1_2.obs[['x_centroid', 'y_centroid']].to_numpy()
TU01plusLB_chuvio_L2_1.obsm['spatial'] = TU01plusLB_chuvio_L2_1.obs[['x_centroid', 'y_centroid']].to_numpy()
TU01plusLB_chuvio_L3_1.obsm['spatial'] = TU01plusLB_chuvio_L3_1.obs[['x_centroid', 'y_centroid']].to_numpy()
TU01plusLB_chuvio_L4_1.obsm['spatial'] = TU01plusLB_chuvio_L4_1.obs[['x_centroid', 'y_centroid']].to_numpy()

TU_MINUS01plusLB_chuvio_L1_1.obsm['spatial'] = TU_MINUS01plusLB_chuvio_L1_1.obs[['x_centroid', 'y_centroid']].to_numpy()
TU_MINUS01plusLB_chuvio_L1_2.obsm['spatial'] = TU_MINUS01plusLB_chuvio_L1_2.obs[['x_centroid', 'y_centroid']].to_numpy()
TU_MINUS01plusLB_chuvio_L2_1.obsm['spatial'] = TU_MINUS01plusLB_chuvio_L2_1.obs[['x_centroid', 'y_centroid']].to_numpy()
TU_MINUS01plusLB_chuvio_L3_1.obsm['spatial'] = TU_MINUS01plusLB_chuvio_L3_1.obs[['x_centroid', 'y_centroid']].to_numpy()
TU_MINUS01plusLB_chuvio_L4_1.obsm['spatial'] = TU_MINUS01plusLB_chuvio_L4_1.obs[['x_centroid', 'y_centroid']].to_numpy()

#Add the raw data
TU01plusLB_chuvio_L1_1.raw = TU01plusLB_chuvio_L1_1
TU01plusLB_chuvio_L1_2.raw = TU01plusLB_chuvio_L1_2
TU01plusLB_chuvio_L2_1.raw = TU01plusLB_chuvio_L2_1
TU01plusLB_chuvio_L3_1.raw = TU01plusLB_chuvio_L3_1
TU01plusLB_chuvio_L4_1.raw = TU01plusLB_chuvio_L4_1

TU_MINUS01plusLB_chuvio_L1_1.raw = TU_MINUS01plusLB_chuvio_L1_1
TU_MINUS01plusLB_chuvio_L1_2.raw = TU_MINUS01plusLB_chuvio_L1_2
TU_MINUS01plusLB_chuvio_L2_1.raw = TU_MINUS01plusLB_chuvio_L2_1
TU_MINUS01plusLB_chuvio_L3_1.raw = TU_MINUS01plusLB_chuvio_L3_1
TU_MINUS01plusLB_chuvio_L4_1.raw = TU_MINUS01plusLB_chuvio_L4_1

#Copy to these samller datasets
TU01plusLB_chuvio_L1_1_dis500 = TU01plusLB_chuvio_L1_1.copy()
TU01plusLB_chuvio_L1_2_dis500 = TU01plusLB_chuvio_L1_2.copy()
TU01plusLB_chuvio_L2_1_dis500 = TU01plusLB_chuvio_L2_1.copy()
TU01plusLB_chuvio_L3_1_dis500 = TU01plusLB_chuvio_L3_1.copy()
TU01plusLB_chuvio_L4_1_dis500 = TU01plusLB_chuvio_L4_1.copy()

TU_MINUS01plusLB_chuvio_L1_1_dis500 = TU_MINUS01plusLB_chuvio_L1_1.copy()
TU_MINUS01plusLB_chuvio_L1_2_dis500 = TU_MINUS01plusLB_chuvio_L1_2.copy()
TU_MINUS01plusLB_chuvio_L2_1_dis500 = TU_MINUS01plusLB_chuvio_L2_1.copy()
TU_MINUS01plusLB_chuvio_L3_1_dis500 = TU_MINUS01plusLB_chuvio_L3_1.copy()
TU_MINUS01plusLB_chuvio_L4_1_dis500 = TU_MINUS01plusLB_chuvio_L4_1.copy()

df_cellchat = ct.pp.ligand_receptor_database(species='human', signaling_type='Secreted Signaling', database='CellChat')
print(df_cellchat.shape)

#Get the type of dataset I need
TU01plusLB_chuvio_L1_1_cellchat_filtered = ct.pp.filter_lr_database(df_cellchat, TU01plusLB_chuvio_L1_1_dis500, min_cell_pct=0.005)
TU01plusLB_chuvio_L1_2_cellchat_filtered = ct.pp.filter_lr_database(df_cellchat, TU01plusLB_chuvio_L1_2_dis500, min_cell_pct=0.005)
TU01plusLB_chuvio_L2_1_cellchat_filtered = ct.pp.filter_lr_database(df_cellchat, TU01plusLB_chuvio_L1_1_dis500, min_cell_pct=0.005)
TU01plusLB_chuvio_L3_1_cellchat_filtered = ct.pp.filter_lr_database(df_cellchat, TU01plusLB_chuvio_L1_2_dis500, min_cell_pct=0.005)
TU01plusLB_chuvio_L4_1_cellchat_filtered = ct.pp.filter_lr_database(df_cellchat, TU01plusLB_chuvio_L1_1_dis500, min_cell_pct=0.005)

TU_MINUS01plusLB_chuvio_L1_1_cellchat_filtered = ct.pp.filter_lr_database(df_cellchat, TU_MINUS01plusLB_chuvio_L1_1_dis500, min_cell_pct=0.005)
TU_MINUS01plusLB_chuvio_L1_2_cellchat_filtered = ct.pp.filter_lr_database(df_cellchat, TU_MINUS01plusLB_chuvio_L1_2_dis500, min_cell_pct=0.005)
TU_MINUS01plusLB_chuvio_L2_1_cellchat_filtered = ct.pp.filter_lr_database(df_cellchat, TU_MINUS01plusLB_chuvio_L2_1_dis500, min_cell_pct=0.005)
TU_MINUS01plusLB_chuvio_L3_1_cellchat_filtered = ct.pp.filter_lr_database(df_cellchat, TU_MINUS01plusLB_chuvio_L3_1_dis500, min_cell_pct=0.005)
TU_MINUS01plusLB_chuvio_L4_1_cellchat_filtered = ct.pp.filter_lr_database(df_cellchat, TU_MINUS01plusLB_chuvio_L4_1_dis500, min_cell_pct=0.005)

ct.tl.spatial_communication(TU01plusLB_chuvio_L1_1_dis500,
    database_name='cellchat', df_ligrec=TU01plusLB_chuvio_L1_1_cellchat_filtered, dis_thr=500, heteromeric=True, pathway_sum=True)

ct.tl.spatial_communication(TU01plusLB_chuvio_L1_2_dis500,
    database_name='cellchat', df_ligrec=TU01plusLB_chuvio_L1_2_cellchat_filtered, dis_thr=500, heteromeric=True, pathway_sum=True)

ct.tl.spatial_communication(TU01plusLB_chuvio_L2_1_dis500,
    database_name='cellchat', df_ligrec=TU01plusLB_chuvio_L2_1_cellchat_filtered, dis_thr=500, heteromeric=True, pathway_sum=True)

ct.tl.spatial_communication(TU01plusLB_chuvio_L3_1_dis500,
    database_name='cellchat', df_ligrec=TU01plusLB_chuvio_L3_1_cellchat_filtered, dis_thr=500, heteromeric=True, pathway_sum=True)

ct.tl.spatial_communication(TU01plusLB_chuvio_L4_1_dis500,
    database_name='cellchat', df_ligrec=TU01plusLB_chuvio_L4_1_cellchat_filtered, dis_thr=500, heteromeric=True, pathway_sum=True)

ct.tl.spatial_communication(TU_MINUS01plusLB_chuvio_L1_1_dis500,
    database_name='cellchat', df_ligrec=TU_MINUS01plusLB_chuvio_L1_1_cellchat_filtered, dis_thr=500, heteromeric=True, pathway_sum=True)

ct.tl.spatial_communication(TU_MINUS01plusLB_chuvio_L1_2_dis500,
    database_name='cellchat', df_ligrec=TU_MINUS01plusLB_chuvio_L1_2_cellchat_filtered, dis_thr=500, heteromeric=True, pathway_sum=True)

ct.tl.spatial_communication(TU_MINUS01plusLB_chuvio_L2_1_dis500,
    database_name='cellchat', df_ligrec=TU_MINUS01plusLB_chuvio_L2_1_cellchat_filtered, dis_thr=500, heteromeric=True, pathway_sum=True)

ct.tl.spatial_communication(TU_MINUS01plusLB_chuvio_L3_1_dis500,
    database_name='cellchat', df_ligrec=TU_MINUS01plusLB_chuvio_L3_1_cellchat_filtered, dis_thr=500, heteromeric=True, pathway_sum=True)

ct.tl.spatial_communication(TU_MINUS01plusLB_chuvio_L4_1_dis500,
    database_name='cellchat', df_ligrec=TU_MINUS01plusLB_chuvio_L4_1_cellchat_filtered, dis_thr=500, heteromeric=True, pathway_sum=True)

TU01plusLB_chuvio_L1_1_dis500.write("data/TU01plusLB_chuvio_L1_1_dis500.h5ad")
TU01plusLB_chuvio_L1_2_dis500.write("data/TU01plusLB_chuvio_L1_2_dis500.h5ad")
TU01plusLB_chuvio_L2_1_dis500.write("data/TU01plusLB_chuvio_L2_1_dis500.h5ad")
TU01plusLB_chuvio_L3_1_dis500.write("data/TU01plusLB_chuvio_L3_1_dis500.h5ad")
TU01plusLB_chuvio_L4_1_dis500.write("data/TU01plusLB_chuvio_L4_1_dis500.h5ad")
TU_MINUS01plusLB_chuvio_L1_1_dis500.write("data/TU_MINUS01plusLB_chuvio_L1_1_dis500.h5ad")
TU_MINUS01plusLB_chuvio_L1_2_dis500.write("data/TU_MINUS01plusLB_chuvio_L1_2_dis500.h5ad")
TU_MINUS01plusLB_chuvio_L2_1_dis500.write("data/TU_MINUS01plusLB_chuvio_L2_1_dis500.h5ad")
TU_MINUS01plusLB_chuvio_L3_1_dis500.write("data/TU_MINUS01plusLB_chuvio_L3_1_dis500.h5ad")
TU_MINUS01plusLB_chuvio_L4_1_dis500.write("data/TU_MINUS01plusLB_chuvio_L4_1_dis500.h5ad")
