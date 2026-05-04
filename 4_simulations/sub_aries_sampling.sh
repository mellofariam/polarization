#!/bin/bash

#SBATCH -J smpl                     # Job name
#SBATCH -o output/job.%j.out        # Output folder
#SBATCH -N 1                        # Number of tasks
#SBATCH -n 8                        # Number of tasks
#SBATCH -t 24:00:00                 # Run time
#SBATCH -a 1,9,17,25                # Array: 1,9,17,25

source /home/mm146/.conda/envs/work/bin/activate
which python

echo "Running sampling of the nucleus with chromosome 1."
echo "Launching 8 jobs on different GPUs..."

condition=$1
i=${SLURM_ARRAY_TASK_ID}

CHROMOSOME="1"

if [ "$CHROMOSOME" == "1" ]; then
    chr_folder="/work/cms16/mm146/Polarization/4_IMR-90_hg38_chr1"
elif [ "$CHROMOSOME" == "17" ]; then
    chr_folder="/work/cms16/mm146/Polarization/5_IMR-90_hg38_chr17"
else
    echo "Invalid chromosome: $CHROMOSOME"
    exit 1
fi

SCRIPT_DIR="/home/mm146/Polarization/polarization/4_simulations/sim_sampling.py"

INPUT_BASE="${chr_folder}/5_assembly"
OUTPUT_BASE="${chr_folder}/6_sampling"

mkdir -p $OUTPUT_BASE/$condition/output

export HIP_VISIBLE_DEVICES=0; srun -n 1 -o $OUTPUT_BASE/$condition/output/job.%J.out python $SCRIPT_DIR ${OUTPUT_BASE} $condition $(($i + 0)) ${INPUT_BASE} HIP &
export HIP_VISIBLE_DEVICES=1; srun -n 1 -o $OUTPUT_BASE/$condition/output/job.%J.out python $SCRIPT_DIR ${OUTPUT_BASE} $condition $(($i + 1)) ${INPUT_BASE} HIP &
export HIP_VISIBLE_DEVICES=2; srun -n 1 -o $OUTPUT_BASE/$condition/output/job.%J.out python $SCRIPT_DIR ${OUTPUT_BASE} $condition $(($i + 2)) ${INPUT_BASE} HIP &
export HIP_VISIBLE_DEVICES=3; srun -n 1 -o $OUTPUT_BASE/$condition/output/job.%J.out python $SCRIPT_DIR ${OUTPUT_BASE} $condition $(($i + 3)) ${INPUT_BASE} HIP &
export HIP_VISIBLE_DEVICES=4; srun -n 1 -o $OUTPUT_BASE/$condition/output/job.%J.out python $SCRIPT_DIR ${OUTPUT_BASE} $condition $(($i + 4)) ${INPUT_BASE} HIP &
export HIP_VISIBLE_DEVICES=5; srun -n 1 -o $OUTPUT_BASE/$condition/output/job.%J.out python $SCRIPT_DIR ${OUTPUT_BASE} $condition $(($i + 5)) ${INPUT_BASE} HIP &
export HIP_VISIBLE_DEVICES=6; srun -n 1 -o $OUTPUT_BASE/$condition/output/job.%J.out python $SCRIPT_DIR ${OUTPUT_BASE} $condition $(($i + 6)) ${INPUT_BASE} HIP &
export HIP_VISIBLE_DEVICES=7; srun -n 1 -o $OUTPUT_BASE/$condition/output/job.%J.out python $SCRIPT_DIR ${OUTPUT_BASE} $condition $(($i + 7)) ${INPUT_BASE} HIP &

echo "Job steps submitted..."
sleep 1
squeue -u `id -un` -s

wait

echo "All Steps completed."
