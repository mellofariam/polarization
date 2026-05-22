#!/bin/bash

#SBATCH -J dist_patch
#SBATCH -o /scratch/mm146/Polarization/output/job.dist_patch.%j.out
#SBATCH -n 4
#SBATCH -t 24:00:00
#SBATCH --account=ctbp-onuchic
#SBATCH --partition=ctbp-onuchic,ctbp-common,commons
#SBATCH --mem=16GB
#SBATCH --array=1-32

source /home/mm146/.conda/envs/work/bin/activate
conda activate work

which python

replica=${SLURM_ARRAY_TASK_ID}

CHROMOSOME="1"

if [ $CHROMOSOME == "1" ]; then
    patches="19 47 75 95 175 205 277 393 456" # A compartment
    tag="A"

    # patches="12 50 88 116 152 202 274 368 372 428 492" # B compartment
    # tag="B"

    traj_folder="/scratch/mm146/Polarization/4_IMR-90_hg38_chr1/6_sampling"
    output_folder="/scratch/mm146/Polarization/4_IMR-90_hg38_chr1/7_analysis/distances-distributions/data"
elif [ $CHROMOSOME == "17" ]; then
    # patches="" # A compartment
    # patches="" # B compartment

    traj_folder="/scratch/mm146/Polarization/5_IMR-90_hg38_chr17/6_sampling"
    output_folder="/scratch/mm146/Polarization/5_IMR-90_hg38_chr17/7_analysis/distances-distributions/data"
else
    echo "Invalid chromosome: $CHROMOSOME. Options are: 'chr1' and 'chr17'."
    exit 1
fi

python compute-patch-distance-distribution.py $condition $replica $traj_folder $output_folder $tag $patches
