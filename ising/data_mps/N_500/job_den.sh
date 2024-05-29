#!/bin/bash -l
#SBATCH --job-name=den_end_500
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --time=48:00:00
#SBATCH --mail-user=amit.gangapuram@lrz.de
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --array=1-9
#SBATCH --output /dss/dsshome1/07/di93hog/quench_dyn/ising/data_mps/N_500/den_output/den_500_exp_%A_%a
#SBATCH --error /dss/dsshome1/07/di93hog/quench_dyn/ising/data_mps/N_500/den_error/den_500_exp_%A_%a

/dss/dsshome1/07/di93hog/julia_1pt6/bin/julia /dss/dsshome1/07/di93hog/quench_dyn/ising/data_mps/N_500/data_ising.jl
