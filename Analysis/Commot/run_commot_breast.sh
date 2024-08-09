#!/bin/bash

###############################
#          HPC setup          #
###############################

##SBATCH --partition gpu
##SBATCH --gres gpu:1
##SBATCH --gres-flags enforce-binding

#SBATCH --partition cpu

#SBATCH --output out

#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=500gb

# Activate the virtual environment
source activate commot2

# Change to the directory containing the Python script
cd /work/PRTNR/CHUV/DIR/rgottar1/spatial/env/astear/Analysis/Commot

# Run the Python script
python commot_breast.py