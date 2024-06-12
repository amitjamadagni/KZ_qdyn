#!/bin/bash -l
#SBATCH --job-name=den_end_500
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --time=48:00:00
#SBATCH --array=1-9

./julia /quench_dyn/ising/data_mps/N_500/data_ising.jl
