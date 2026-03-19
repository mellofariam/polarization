#!/bin/bash

#SBATCH -J dist_ij
#SBATCH -o /scratch/mm146/Polarization/3_IMR-90_hg38_fields-nb/4_analysis/output/job.dist_ij.%j.out
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH --account=ctbp-onuchic
#SBATCH --partition=ctbp-onuchic,ctbp-common,commons
#SBATCH --mem=16GB
#SBATCH --array=1-32

source /opt/apps/software/Mamba/23.1.0-4/bin/activate
conda activate "$HOME/work"

which python

replica=${SLURM_ARRAY_TASK_ID}

python compute-distance-histogram.py $condition $replica 159 3119