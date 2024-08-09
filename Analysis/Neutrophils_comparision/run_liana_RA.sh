#!/bin/bash

###############################
#          HPC setup          #
###############################

##SBATCH --partition gpu
##SBATCH --gres gpu:1
##SBATCH --gres-flags enforce-binding

#SBATCH --partition cpu

#SBATCH --output out

#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=100gb

# Activate the virtual environment
source activate CellPhoneDB

# Change to the directory containing the Python script
cd /work/PRTNR/CHUV/DIR/rgottar1/spatial/env/astear/Analysis/Neutrophils_comparision

# Run the Python script
python Liana_chuvio_L4_1.py