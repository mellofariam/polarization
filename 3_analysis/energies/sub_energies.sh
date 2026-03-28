#!/bin/bash

#SBATCH -J energies
#SBATCH -o /scratch/mm146/Polarization/5_IMR-90_hg38_chr17/4_analysis/output/job.energies.%j.out
#SBATCH -n 4
#SBATCH -t 24:00:00
#SBATCH --account=ctbp-onuchic
#SBATCH --partition=ctbp-onuchic,ctbp-common,commons
#SBATCH --mem=8GB
#SBATCH --array=1-32

module load Mamba/23.1.0-4
source /opt/apps/software/Mamba/23.1.0-4/bin/activate
conda activate "$HOME/work"

replica=${SLURM_ARRAY_TASK_ID}

python compute-energies.py $condition $replica /scratch/mm146/Polarization/5_IMR-90_hg38_chr17/3_sampling /scratch/mm146/Polarization/5_IMR-90_hg38_chr17/4_analysis/energies/data
