#!/bin/bash

#SBATCH -J energies
#SBATCH -o /scratch/mm146/Polarization/3_IMR-90_hg38_fields-nb/4_analysis/output/job.energies.%j.out
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

python compute-energies.py $condition $replica
