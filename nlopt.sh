#!/bin/bash

# module load mpi/openmpi-4.0
# SBATCH --gres=gpu:volta:1
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -o logs/sg98-opt-global-full-%j-%a.out
#SBATCH -a 1-5
# 
# module load anaconda/2020b
# conda init
# conda activate

# echo $SLURM_ARRAY_TASK_ID
# echo $LLSUB_RANK

# julia nlopt_gap_global.jl $SLURM_ARRAY_TASK_ID
# julia nlopt_degen.jl $SLURM_ARRAY_TASK_ID
# julia nlopt_hinge.jl $SLURM_ARRAY_TASK_ID 147
