#!/bin/bash

#SBATCH -J polarization
#SBATCH -o /scratch/mm146/Polarization/output/job.polarization.%j.out
#SBATCH -n 4
#SBATCH -t 24:00:00
#SBATCH --account=commons
#SBATCH --partition=commons
#SBATCH --mem=8GB
#SBATCH --array=1-32

source /home/mm146/.conda/envs/work/bin/activate
conda activate work

which python

# remember to define condition and strength with --export=condition=$condition,strength=$strength when calling sbatch
replica=${SLURM_ARRAY_TASK_ID}

chr_folder="/scratch/mm146/Polarization/4_IMR-90_hg38_chr1"

python compute-polarization.py \
    $condition \
    $replica \
    $chr_folder/8_sampling/$strength \
    $chr_folder/9_analysis/$strength/polarization/data
