#!/bin/bash

## slurm options

#SBATCH -p intel-sc3         
#SBATCH -q normal             # normal qos
#SBATCH -J scmethy              # job name
#SBATCH -e hcc73.log
#SBATCH -N 4
#SBATCH -n 32
#SBATCH --mem=200G
#SBATCH -a 1-6


hostname
module load R/4.0.5
R
Rscript --no-save '/storage/zhangyanxiaoLab/zhangliwen/Projects/sc-LIHC/analysis/scMehty/hcc7.R'

