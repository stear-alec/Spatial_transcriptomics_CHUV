#!/bin/bash

###############################
#          HPC setup          #
###############################

##SBATCH --partition gpu
##SBATCH --gres gpu:1
##SBATCH --gres-flags enforce-binding

#SBATCH --partition cpu

#SBATCH --output out

#SBATCH --time=15:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=200gb

# Activate the virtual env
source ~/.bashrc
conda activate R_scdesign

# Change to the directory containing the Python script
cd /work/PRTNR/CHUV/DIR/rgottar1/spatial/env/astear/Programming/synthetic_data/subTulb

# Run the R script
Rscript LBLT_l1_1_scDesign3.R