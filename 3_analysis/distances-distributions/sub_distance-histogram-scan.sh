#!/bin/bash

#SBATCH -J dist-hist
#SBATCH -o /scratch/mm146/Polarization/output/job.dist-hist.%j.out
#SBATCH -n 4
#SBATCH -t 24:00:00
#SBATCH --account=ctbp-onuchic
#SBATCH --partition=ctbp-onuchic,ctbp-common,commons
#SBATCH --mem=64GB
#SBATCH --array=1-32

source /home/mm146/.conda/envs/work/bin/activate
conda activate work

which python

replica=${SLURM_ARRAY_TASK_ID}

traj_folder="/scratch/mm146/Polarization/4_IMR-90_hg38_chr1/8_sampling/$strength"
output_folder="/scratch/mm146/Polarization/4_IMR-90_hg38_chr1/9_analysis/$strength/distances-distributions/data"

python compute-distance-histogram.py $condition $replica $traj_folder $output_folder
    
