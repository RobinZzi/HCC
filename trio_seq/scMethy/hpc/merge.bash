#!/bin/bash

## slurm options

#SBATCH -p intel-sc3         
#SBATCH -q normal             # normal qos
#SBATCH -J scmethy              # job name
#SBATCH -e 2.log
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=40G
#SBATCH -a 1-6


hostname

R
Rscript --no-save '/storage/zhangyanxiaoLab/zhangliwen/Projects/sc-LIHC/analysis/scMehty/merge.R'

