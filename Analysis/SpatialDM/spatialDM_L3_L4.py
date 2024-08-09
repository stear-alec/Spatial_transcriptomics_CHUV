import os
import gc
import pickle
import anndata
import scanpy as sc
import pandas as pd
import numpy as np
from scipy import sparse
from scipy.stats import spearmanr, pearsonr
from scipy.spatial import distance_matrix
import matplotlib.pyplot as plt
import spatialdm as sdm

file = "breast_B1_1"

path = f"/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/xenium/xenium_objs/anndata/{file}.h5ad"
file = sc.read_h5ad(path)

file.obsm['spatial'] = chuvio_B1_1.obs[['x_centroid', 'y_centroid']].to_numpy()

file.raw = file

#running the method
sdm.weight_matrix(file, l=120, cutoff=0.2, single_cell=False) # weight_matrix by rbf kernel
sdm.extract_lr(file, 'human', min_cell=1)      # find overlapping LRs from CellChatDB
sdm.spatialdm_global(file, 100, specified_ind=None, method='both', nproc=1)     # global Moran selection (orig 1000 set to 100)
sdm.sig_pairs(file, method='permutation', fdr=True, threshold=0.1)     # select significant pairs
sdm.spatialdm_local(file, n_perm=100, method='both', specified_ind=None, nproc=1)     # local spot selection
sdm.sig_spots(file, method='permutation', fdr=False, threshold=0.1)     # significant local spots

breast_B1_1_sDM = file

file = file.uns['global_res'].sort_values(by='fdr')





sdm.weight_matrix(chuvio_L4_1, l=120, cutoff=0.2, single_cell=False) # weight_matrix by rbf kernel
sdm.extract_lr(chuvio_L4_1, 'human', min_cell=1)      # find overlapping LRs from CellChatDB
sdm.spatialdm_global(chuvio_L4_1, 100, specified_ind=None, method='both', nproc=1)     # global Moran selection (orig 1000 set to 100)
sdm.sig_pairs(chuvio_L4_1, method='permutation', fdr=True, threshold=0.1)     # select significant pairs
sdm.spatialdm_local(chuvio_L4_1, n_perm=100, method='both', specified_ind=None, nproc=1)     # local spot selection
sdm.sig_spots(chuvio_L4_1, method='permutation', fdr=False, threshold=0.1)     # significant local spots

chuvio_L4_1_sDM = chuvio_L4_1

#getting the top genes
chuvio_L3_1_sDM_top_G = chuvio_L3_1_sDM.uns['global_res'].sort_values(by='fdr')
chuvio_L4_1_sDM_top_G = chuvio_L4_1_sDM.uns['global_res'].sort_values(by='fdr')

chuvio_L3_1_sDM_top_G = chuvio_L3_1_sDM.uns['global_res'].sort_values(by='fdr')
chuvio_L4_1_sDM_top_G.to_csv('chuvio_L4_1_sDM_top_G.csv', index=False)

