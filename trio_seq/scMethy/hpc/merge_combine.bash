#!/bin/bash

## slurm options

#SBATCH -p intel-sc3         
#SBATCH -q normal             # normal qos
#SBATCH -J scmethy              # job name
#SBATCH -e merge.log
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=80G



hostname
module load R/4.0.5
R
Rscript --no-save '/storage/zhangyanxiaoLab/zhangliwen/Projects/sc-LIHC/analysis/scMehty/merge_combine.R'

