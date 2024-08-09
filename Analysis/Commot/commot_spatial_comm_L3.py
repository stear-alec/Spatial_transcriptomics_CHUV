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


path_chuvio_L3_1 = "/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/xenium/xenium_objs/anndata/chuvio_L3_1.h5ad"
chuvio_L3_1 = sc.read_h5ad(path_chuvio_L3_1)


chuvio_L3_1.obsm['spatial'] = chuvio_L3_1.obs[['x_centroid', 'y_centroid']].to_numpy()



chuvio_L3_1.raw = chuvio_L3_1

#Copying, not sure this step is necessary
chuvio_L3_1_dis500 = chuvio_L3_1.copy()

#Using the CellPhoneDB database ot find Ligand Receptor genes
df_cellchat = ct.pp.ligand_receptor_database(species='human', signaling_type='Secreted Signaling', database='CellChat')

#Running Algorithm
#ct.tl.spatial_communication(chuvio_L1_1_dis500,
    #database_name='cellchat', df_ligrec=df_cellchat, dis_thr=500, heteromeric=True, pathway_sum=True)
#chuvio_L1_1_dis500.write("data/chuvio_L1_1_dis500.h5ad")

#ct.tl.spatial_communication(chuvio_L1_2_dis500,
    #database_name='cellchat', df_ligrec=df_cellchat, dis_thr=500, heteromeric=True, pathway_sum=True)
#chuvio_L1_2_dis500.write("data/chuvio_L1_2_dis500.h5ad")

#ct.tl.spatial_communication(chuvio_L2_1_dis500,
    #database_name='cellchat', df_ligrec=df_cellchat, dis_thr=500, heteromeric=True, pathway_sum=True)
#chuvio_L2_1_dis500.write("data/chuvio_L2_1_dis500.h5ad")

ct.tl.spatial_communication(chuvio_L3_1_dis500,
    database_name='cellchat', df_ligrec=df_cellchat, dis_thr=500, heteromeric=True, pathway_sum=True)
chuvio_L3_1_dis500.write("data/chuvio_L3_1_dis500.h5ad")







