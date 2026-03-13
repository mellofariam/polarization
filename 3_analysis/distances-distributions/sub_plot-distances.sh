#!/bin/bash

#SBATCH -J plt.dist             # Job name
#SBATCH -o output/job.%j.out    # Name of stdout output file (%j expands to jobId)
#SBATCH -n 1                    # Number of tasks
#SBATCH -t 24:00:00             # Run time
#SBATCH -p commons              # Desired partition
#SBATCH --account commons       # Desired account
#SBATCH --mem=64GB

source /home/mm146/.conda/envs/work/bin/activate
conda activate work

python plot-distances-$target.py
