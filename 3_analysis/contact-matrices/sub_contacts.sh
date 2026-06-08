#!/bin/bash

#SBATCH -J contacts
#SBATCH -o /scratch/mm146/Polarization/output/job.contacts.%j.out
#SBATCH -n 4
#SBATCH -t 24:00:00
#SBATCH --account=ctbp-onuchic
#SBATCH --partition=ctbp-onuchic,ctbp-common,commons
#SBATCH --mem=16GB
#SBATCH --array=1-96

source /home/mm146/.conda/envs/work/bin/activate
conda activate work

# remember to define condition with --export=condition=$condition when calling sbatch
replica=${SLURM_ARRAY_TASK_ID}

chromosome="1"

if [ "$chromosome" == "1" ]; then
    chr_folder="/scratch/mm146/Polarization/4_IMR-90_hg38_chr1"
elif [ "$chromosome" == "17" ]; then
    chr_folder="/scratch/mm146/Polarization/5_IMR-90_hg38_chr17"
else
    echo "Invalid chromosome: $chromosome"
    exit 1
fi

python compute-contact-matrices.py \
    $condition \
    $replica \
    $chr_folder/6_sampling \
    $chr_folder/7_analysis/contacts/data
