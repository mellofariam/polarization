#!/bin/bash

#SBATCH -J lamina.contact
#SBATCH -o /scratch/mm146/Polarization/output/job.lamina-contact.%j.out
#SBATCH -n 4
#SBATCH -t 24:00:00
#SBATCH --account=ctbp-onuchic
#SBATCH --partition=ctbp-onuchic,ctbp-common,commons
#SBATCH --mem=8GB
#SBATCH --array=1-32

source /home/mm146/.conda/envs/work/bin/activate
conda activate work

which python

# remember to define `condition` and `strength` with --export=condition=$condition,strength=$strength when calling sbatch:
# sbatch --export=condition=$condition,strength=$strength sub_lamina-contact.sh
replica=${SLURM_ARRAY_TASK_ID}

nuclear_body_treatment="beads"  # "fields" or "beads"

chr_folder="/scratch/mm146/Polarization/4_IMR-90_hg38_chr1"

python compute-lamina-contact.py \
    $condition \
    $replica \
    $chr_folder/8_sampling/$strength \
    $chr_folder/9_analysis/$strength/lamina-contact/data \
    $nuclear_body_treatment
