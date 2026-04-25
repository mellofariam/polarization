#!/bin/bash

#SBATCH -J asbly                # Job name
#SBATCH -o output/job.%j.out    # Name of stdout output file (%j expands to jobId)
#SBATCH -n 1                    # Number of tasks
#SBATCH -t 24:00:00             # Run time
#SBATCH -p commons              # Desired partition
#SBATCH --account commons       # Desired account
#SBATCH --gres=gpu:1            # Request GPU
#SBATCH --array=1-32

source /home/mm146/.conda/envs/work/bin/activate
conda activate work

which python

condition=$1
i=${SLURM_ARRAY_TASK_ID}

SCRIPT_DIR="/home/mm146/Polarization/polarization/4_simulations/sim_assembly.py"
OUTPUT_BASE="/work/cms16/mm146/Polarization/4_IMR-90_hg38_chr17/2_assembly"
INPUT_BASE="/work/cms16/mm146/Polarization/4_IMR-90_hg38_chr17/1_chromosome-collapse"

mkdir -p $OUTPUT_BASE/$condition/output

python $SCRIPT_DIR ${OUTPUT_BASE} 17 $condition ${i} ${INPUT_BASE} CUDA \
    > $OUTPUT_BASE/$condition/output/job.${SLURM_JOB_ID}.${i}.out
