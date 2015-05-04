#!/bin/bash

#SBATCH -J EQKSearch_results_SAMC_B1000_avg_Eikonal
#SBATCH -o results_SAMC_B1000_avg_Eikonal/log_EQKSearch_%j.out
#SBATCH -e results_SAMC_B1000_avg_Eikonal/log_EQKSearch_%j.err
#SBATCH -N 1
#SBATCH --ntasks-per-node=12
#SBATCH --time=15:00:00
#SBATCH --mem=1gb

ulimit -c unlimited
module load gcc/4.8.0

export OMP_NUM_THREADS=12; time /projects/yeti4009/eqkhyposolver/EQKSolver results_SAMC_B1000_avg_Eikonal/param.txt -c 1> results_SAMC_B1000_avg_Eikonal/log_EQKSearch.out 2> results_SAMC_B1000_avg_Eikonal/log_EQKSearch.err

