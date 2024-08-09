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



file = "breast_B4_1"

file_name = file
path = f"/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/xenium/xenium_objs/anndata/{file}.h5ad"
file = sc.read_h5ad(path)
file.obsm['spatial'] = file.obs[['x_centroid', 'y_centroid']].to_numpy()
file.raw = file
#Copying, not sure this step is necessary
file_dis500 = file.copy()
#Using the CellPhoneDB database ot find Ligand Receptor genes
df_cellchat = ct.pp.ligand_receptor_database(species='human', signaling_type='Secreted Signaling', database='CellChat')
ct.tl.spatial_communication(file_dis500,
    database_name='cellchat', df_ligrec=df_cellchat, dis_thr=500, heteromeric=True, pathway_sum=True)
file_dis500.write(f"data/{file_name}_dis500.h5ad")






