# import liana
import liana as li
# needed for visualization and toy data
import scanpy as sc

#Import Dataset
path_file = "/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/astear/Analysis/Neutrophils_comparision/data_set/chuvio_L4_1_VISTA.h5ad"

from anndata import read_h5ad
import sys
adata = read_h5ad(path_file)
adata.raw = adata

from liana.method import singlecellsignalr, connectome, cellphonedb, natmi, logfc, cellchat, geometric_mean

def Liana_to_output(groupy):
    global_vars = globals()
    liana_output = adata.uns['liana_res']
    global_vars[f"liana_output_{groupy}"] = liana_output
    global_vars[f"significant_interactions_{groupy}"] = len(liana_output)
    global_vars[f"list_targets_{groupy}"] = liana_output['target'].unique()

# Run rank_aggregate
li.mt.rank_aggregate(adata, groupby='New_annotation', expr_prop=0.1, verbose=True)

Liana_to_output(groupy="New_annotation")
#Make sure the above groupby is used in indicating the object to be saved as csv in the line beneath
liana_output_New_annotation.to_csv("liana_output_chuvio_L4_1_VISTA.csv", index=False)