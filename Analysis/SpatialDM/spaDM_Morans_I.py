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

import spatialdm.plottings as pl

import matplotlib.pyplot as plt


path = "/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/xenium/xenium_objs/anndata/breast_B1_1.h5ad"
chuvio_B1_1 = sc.read_h5ad(path)
path = "/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/xenium/xenium_objs/anndata/breast_B1_2.h5ad"
chuvio_B1_2 = sc.read_h5ad(path)
path = "/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/xenium/xenium_objs/anndata/breast_B2_1.h5ad"
chuvio_B2_1 = sc.read_h5ad(path)
path = "/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/xenium/xenium_objs/anndata/breast_B3_1.h5ad"
chuvio_B3_1 = sc.read_h5ad(path)
path = "/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/xenium/xenium_objs/anndata/breast_B4_1.h5ad"
chuvio_B4_1 = sc.read_h5ad(path)

chuvio_B1_1.obsm['spatial'] = chuvio_B1_1.obs[['x_centroid', 'y_centroid']].to_numpy()
chuvio_B1_2.obsm['spatial'] = chuvio_B1_2.obs[['x_centroid', 'y_centroid']].to_numpy()
chuvio_B2_1.obsm['spatial'] = chuvio_B2_1.obs[['x_centroid', 'y_centroid']].to_numpy()
chuvio_B3_1.obsm['spatial'] = chuvio_B3_1.obs[['x_centroid', 'y_centroid']].to_numpy()
chuvio_B4_1.obsm['spatial'] = chuvio_B4_1.obs[['x_centroid', 'y_centroid']].to_numpy()

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

chuvio_L1_1.obsm['spatial'] = chuvio_L1_1.obs[['x_centroid', 'y_centroid']].to_numpy()
chuvio_L1_2.obsm['spatial'] = chuvio_L1_2.obs[['x_centroid', 'y_centroid']].to_numpy()
chuvio_L2_1.obsm['spatial'] = chuvio_L2_1.obs[['x_centroid', 'y_centroid']].to_numpy()
chuvio_L3_1.obsm['spatial'] = chuvio_L3_1.obs[['x_centroid', 'y_centroid']].to_numpy()
chuvio_L4_1.obsm['spatial'] = chuvio_L4_1.obs[['x_centroid', 'y_centroid']].to_numpy()



adata = chuvio_L4_1

adata.raw = adata
#running the method
sdm.weight_matrix(adata, l=120, cutoff=0.2, single_cell=False) # weight_matrix by rbf kernel
sdm.extract_lr(adata, 'human', min_cell=1)      # find overlapping LRs from CellChatDB
sdm.spatialdm_global(adata, 100, specified_ind=None, method='both', nproc=1)     # global Moran selection (orig 1000 set to 100)
sdm.sig_pairs(adata, method='permutation', fdr=True, threshold=0.1)     # select significant pairs
sdm.spatialdm_local(adata, n_perm=100, method='both', specified_ind=None, nproc=1)     # local spot selection
sdm.sig_spots(adata, method='permutation', fdr=False, threshold=0.1)     # significant local spots

plt.figure(figsize=(20,10))
for i in range(1):
    plt.subplot(2, 2, i + 2)
    plt.scatter(adata.obsm['spatial'][:,0], adata.obsm['spatial'][:,1], marker = 's', s=1);
    plt.axis('equal')
    pl.plot_pairs(adata, ['EREG_EGFR'], marker='s', s=0.3, figsize=(60, 10))
plt.savefig('chuvio_L4_1_EREG_EGFR.pdf')