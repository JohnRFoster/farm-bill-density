#!/usr/bin/env bash
#SBATCH --ntasks=1                      # number of processes to run
#SBATCH --nodes=1                       # default number of nodes
#SBATCH --partition=cpu_compute         # good enough for what I need
#SBATCH --cpus-per-task=5               # for a multithredded job
#SBATCH --mem=24g                       # memory
#SBATCH --job-name=farmBillSim             # job name
#SBATCH --output=outfiles/farmBillSim_%a.txt    # output file

module add R
Rscript workflow.R $SLURM_ARRAY_TASK_ID
