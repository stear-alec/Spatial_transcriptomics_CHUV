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

file_name = file
path = f"/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/xenium/xenium_objs/anndata/{file}.h5ad"
file = sc.read_h5ad(path)
file.obsm['spatial'] = file.obs[['x_centroid', 'y_centroid']].to_numpy()
file.raw = file
#running the method
sdm.weight_matrix(file, l=120, cutoff=0.2, single_cell=False) # weight_matrix by rbf kernel
sdm.extract_lr(file, 'human', min_cell=1)      # find overlapping LRs from CellChatDB
sdm.spatialdm_global(file, 100, specified_ind=None, method='both', nproc=1)     # global Moran selection (orig 1000 set to 100)
sdm.sig_pairs(file, method='permutation', fdr=True, threshold=0.1)     # select significant pairs
sdm.spatialdm_local(file, n_perm=100, method='both', specified_ind=None, nproc=1)     # local spot selection
sdm.sig_spots(file, method='permutation', fdr=False, threshold=0.1)     # significant local spots
file = file.uns['global_res'].sort_values(by='fdr')
file.to_csv(f'{file_name}_sDM_top_G.csv', index=False)

