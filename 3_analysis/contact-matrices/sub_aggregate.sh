#!/bin/bash

#SBATCH -J contacts
#SBATCH -o /scratch/mm146/Polarization/output/job.agg-contacts.%j.out
#SBATCH -n 4
#SBATCH -t 24:00:00
#SBATCH --account=commons
#SBATCH --partition=commons
#SBATCH --mem=500GB

source /home/mm146/.conda/envs/work/bin/activate
conda activate work

python aggregate-contact-data.py 
