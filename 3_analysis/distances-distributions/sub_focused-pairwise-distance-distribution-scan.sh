#!/bin/bash

#SBATCH -J dist_ij
#SBATCH -o /scratch/mm146/Polarization/output/job.dist_ij.%j.out
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH --account=commons
#SBATCH --partition=commons
#SBATCH --mem=4GB
#SBATCH --array=1-32

source /home/mm146/.conda/envs/work/bin/activate
conda activate work

which python

replica=${SLURM_ARRAY_TASK_ID}

CHROMOSOME="1"

if [ $CHROMOSOME == "1" ]; then
    focused_pairs=(
        "159 3119"
        "469 1059"
        "859 1059"
        "859 1859"
    )
    traj_folder="/scratch/mm146/Polarization/4_IMR-90_hg38_chr1/8_sampling/$strength"
    output_folder="/scratch/mm146/Polarization/4_IMR-90_hg38_chr1/9_analysis/$strength/distances-distributions/data"
else
    echo "Invalid chromosome: $CHROMOSOME. Options are: 'chr1'."
    exit 1
fi

# Loop over pairs
for pair in "${focused_pairs[@]}"; do
    i=$(echo $pair | awk '{print $1}')
    j=$(echo $pair | awk '{print $2}')

    python compute-focused-pairwise-distance-distribution.py $condition $replica $traj_folder $output_folder $i $j
done
