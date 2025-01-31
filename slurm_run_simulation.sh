#!/usr/bin/env bash
#SBATCH --ntasks=1                      # number of processes to run
#SBATCH --nodes=1                       # default number of nodes
#SBATCH --partition=cpu_compute         # good enough for what I need
#SBATCH --cpus-per-task=5               # for a multithredded job
#SBATCH --mem=96g                       # memory
#SBATCH --job-name=clusterSim             # job name
#SBATCH --output=outfiles/simulation_%J.txt    # output file

module add R
Rscript workflow.R
