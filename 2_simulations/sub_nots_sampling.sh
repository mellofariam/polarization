#!/bin/bash

#SBATCH -J sim1chr              # Job name
#SBATCH -o output/job.%j.out    # Name of stdout output file (%j expands to jobId)
#SBATCH -n 1                    # Number of tasks
#SBATCH -t 24:00:00             # Run time
#SBATCH -p commons              # Desired partition
#SBATCH --account commons       # Desired account
#SBATCH --gres=gpu:1            # request GPU
#SBATCH --array=1-30

# Launch an MPI-based executable
source /home/mm146/.conda/envs/work/bin/activate
conda activate work

which python

i=${SLURM_ARRAY_TASK_ID}

python single-chr-simulation.py $condition $i
