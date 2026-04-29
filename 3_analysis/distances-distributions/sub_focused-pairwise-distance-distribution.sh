#!/bin/bash

#SBATCH -J dist_ij
#SBATCH -o /scratch/mm146/Polarization/output/job.dist_ij.%j.out
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH --account=ctbp-onuchic
#SBATCH --partition=ctbp-onuchic,ctbp-common,commons
#SBATCH --mem=4GB
#SBATCH --array=1-32

source /home/mm146/.conda/envs/work/bin/activate
conda activate work

which python

replica=${SLURM_ARRAY_TASK_ID}

CHROMOSOME="17"

if [ $CHROMOSOME == "1" ]; then
    focused_pairs=(
        "159 3119"
        "469 1059"
        "859 1059"
        "859 1859"
    )
    traj_folder="/scratch/mm146/Polarization/4_IMR-90_hg38_chr1/6_sampling"
    output_folder="/scratch/mm146/Polarization/4_IMR-90_hg38_chr1/7_analysis/distances-distributions/data"
elif [ $CHROMOSOME == "17" ]; then
    focused_pairs=(
        "79 349"
        "349 859"
        "859 1629"
        "79 1629"
        "679 219"
        "679 1079"
        "219 1079"
        "219 1419"
    )
    traj_folder="/scratch/mm146/Polarization/5_IMR-90_hg38_chr17/6_sampling"
    output_folder="/scratch/mm146/Polarization/5_IMR-90_hg38_chr17/7_analysis/distances-distributions/data"
else
    echo "Invalid chromosome: $CHROMOSOME. Options are: 'chr1' and 'chr17'."
    exit 1
fi

# Loop over pairs
for pair in "${focused_pairs[@]}"; do
    i=$(echo $pair | awk '{print $1}')
    j=$(echo $pair | awk '{print $2}')

    python compute-focused-pairwise-distance-distribution.py $condition $replica $traj_folder $output_folder $i $j
done
