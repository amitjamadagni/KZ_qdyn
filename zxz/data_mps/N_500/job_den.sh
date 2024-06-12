#!/bin/bash -l
#SBATCH --job-name=den_end_500
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --array=1

./julia quench_dyn/zxz/data_mps/N_500/data_zxz.jl
