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

CHROMOSOME="17"

if [ $CHROMOSOME == "1" ]; then
    traj_folder="/scratch/mm146/Polarization/4_IMR-90_hg38_chr1/6_sampling"
    output_folder="/scratch/mm146/Polarization/4_IMR-90_hg38_chr1/7_analysis/distances-distributions/data"
elif [ $CHROMOSOME == "17" ]; then
    traj_folder="/scratch/mm146/Polarization/5_IMR-90_hg38_chr17/6_sampling"
    output_folder="/scratch/mm146/Polarization/5_IMR-90_hg38_chr17/7_analysis/distances-distributions/data"
else
    echo "Invalid chromosome: $CHROMOSOME. Options are: 'chr1' and 'chr17'."
    exit 1
fi

python compute-distance-histogram.py $condition $replica $traj_folder $output_folder
    
