#library
import os

import pandas as pd
import numpy as np

import cell2cell as c2c

import warnings
warnings.filterwarnings('ignore')

#Data
data_folder = '/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/astear/Analysis/TensorC2C/data/liana-outputs/'
output_folder = '/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/astear/Analysis/TensorC2C/data/tc2c-outputs/'

#fetching the previous made tensors
tensor = c2c.io.read_data.load_tensor(os.path.join(output_folder, 'BALF-Tensor-5chuv_TU_LB_23_04_24.pkl'))
meta_tensor = c2c.io.load_variable_with_pickle(os.path.join(output_folder, 'BALF-Tensor-5chuv_TU_LB_Metadata_23_04_24.pkl'))

#Running the deconvolution
tensor2 = c2c.analysis.run_tensor_cell2cell_pipeline(tensor,
                                                     meta_tensor,
                                                     copy_tensor=True, # Whether to output a new tensor or modifying the original
                                                     rank=None, # Number of factors to perform the factorization. If None, it is automatically determined by an elbow analysis
                                                     tf_optimization='robust', # To define how robust we want the analysis to be.
                                                     random_state=0, # Random seed for reproducibility
                                                     device='cpu', # Device to use. If using GPU and PyTorch, use 'cuda'. For CPU use 'cpu'
                                                     elbow_metric='error', # Metric to use in the elbow analysis.
                                                     smooth_elbow=False, # Whether smoothing the metric of the elbow analysis.
                                                     upper_rank=25, # Max number of factors to try in the elbow analysis
                                                     tf_init='random', # Initialization method of the tensor factorization
                                                     tf_svd='numpy_svd', # Type of SVD to use if the initialization is 'svd'
                                                     cmaps=None, # Color palettes to use in color each of the dimensions. Must be a list of palettes.
                                                     sample_col='Element', # Columns containing the elements in the tensor metadata
                                                     group_col='Category', # Columns containing the major groups in the tensor metadata
                                                     fig_fontsize=14, # Fontsize of the figures generated
                                                     output_folder=output_folder, # Whether to save the figures in files. If so, a folder pathname must be passed
                                                     output_fig=True, # Whether to output the figures. If False, figures won't be saved a files if a folder was passed in output_folder.
                                                     fig_format='pdf', # File format of the figures.
                                                    )

#Save this tensor file
c2c.io.export_variable_with_pickle(tensor2, '/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/astear/Analysis/TensorC2C/data/tc2c-outputs/L1_TU_LB/BALF-Tensor-5chuv_TU_LBOut_23_04_24.pkl')