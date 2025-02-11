#!/usr/bin/env bash
#SBATCH --ntasks=1                      # number of processes to run
#SBATCH --nodes=1                       # default number of nodes
#SBATCH --partition=cpu_compute         # good enough for what I need
#SBATCH --cpus-per-task=5               # for a multithredded job
#SBATCH --mem=48g                       # memory
#SBATCH --job-name=farmBillSimArray             # job name
#SBATCH --output=outfiles/farmBillSimArray_%a.txt    # output file
#SBATCH --array=1-10

module add R
Rscript workflow.R $SLURM_ARRAY_TASK_ID
