#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
R CMD BATCH --no-save FA_src.R FA_src.out