#!/usr/bin/env bash
#SBATCH --ntasks=1                      # number of processes to run
#SBATCH --nodes=1                       # default number of nodes
#SBATCH --partition=cpu_compute         # good enough for what I need
#SBATCH --job-name=traceplots            # job name
#SBATCH --output=outfiles/traceplots.txt    # output file

module add R
Rscript R/make_simulation_traceplots.R
