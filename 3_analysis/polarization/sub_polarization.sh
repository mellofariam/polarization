#!/bin/bash

#SBATCH -J polarization         # Job name
#SBATCH -o output/job.%j.out    # Name of stdout output file (%j expands to jobId)
#SBATCH -n 1                    # Number of tasks
#SBATCH -t 24:00:00             # Run time
#SBATCH -p commons              # Desired partition
#SBATCH --account commons       # Desired account
#SBATCH --mem=64GB
#SBATCH --array=61-90

source /home/mm146/.conda/envs/work/bin/activate
conda activate work
which python

replica=${SLURM_ARRAY_TASK_ID}

python compute-polarization.py $chr_number $condition $replica
