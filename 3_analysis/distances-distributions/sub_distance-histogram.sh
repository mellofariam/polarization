#!/bin/bash

#SBATCH -J dist-hist
#SBATCH -o /scratch/mm146/Polarization/5_IMR-90_hg38_chr17/4_analysis/output/job.dist-hist.%j.out
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

python compute-distance-histogram.py \\
    $condition \\
    $replica \\
    /scratch/mm146/Polarization/5_IMR-90_hg38_chr17/3_sampling \\
    /scratch/mm146/Polarization/5_IMR-90_hg38_chr17/4_analysis/distances-distributions/data
