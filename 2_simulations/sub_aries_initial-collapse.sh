#!/bin/bash

#SBATCH -J clpse                    # Job name
#SBATCH -o output/job.%j.out        # Output folder
#SBATCH -N 1                        # Number of tasks
#SBATCH -n 8                        # Number of tasks
#SBATCH -t 24:00:00                 # Run time
#SBATCH -a 1,9,17,25                # Array: 1,9,17,25

source /home/mm146/.conda/envs/work/bin/activate
which python

echo "Running initial collapse for chromosome 17."

echo "Launching 8 jobs on different GPUs..."

i=${SLURM_ARRAY_TASK_ID}

SCRIPT_DIR="/home/mm146/Polarization/polarization/2_simulations/sim_initial-collapse.py"
OUTPUT_DIR="/work/cms16/mm146/Polarization/4_IMR-90_hg38_chr17/1_chromosome-collapse"


mkdir -p $OUTPUT_DIR/output

export HIP_VISIBLE_DEVICES=0; srun -n 1 -o $OUTPUT_DIR/output/job.%J.out python $SCRIPT_DIR ${OUTPUT_DIR}/$(($i + 0)) 17 HIP &
export HIP_VISIBLE_DEVICES=1; srun -n 1 -o $OUTPUT_DIR/output/job.%J.out python $SCRIPT_DIR ${OUTPUT_DIR}/$(($i + 1)) 17 HIP &
export HIP_VISIBLE_DEVICES=2; srun -n 1 -o $OUTPUT_DIR/output/job.%J.out python $SCRIPT_DIR ${OUTPUT_DIR}/$(($i + 2)) 17 HIP &
export HIP_VISIBLE_DEVICES=3; srun -n 1 -o $OUTPUT_DIR/output/job.%J.out python $SCRIPT_DIR ${OUTPUT_DIR}/$(($i + 3)) 17 HIP &
export HIP_VISIBLE_DEVICES=4; srun -n 1 -o $OUTPUT_DIR/output/job.%J.out python $SCRIPT_DIR ${OUTPUT_DIR}/$(($i + 4)) 17 HIP &
export HIP_VISIBLE_DEVICES=5; srun -n 1 -o $OUTPUT_DIR/output/job.%J.out python $SCRIPT_DIR ${OUTPUT_DIR}/$(($i + 5)) 17 HIP &
export HIP_VISIBLE_DEVICES=6; srun -n 1 -o $OUTPUT_DIR/output/job.%J.out python $SCRIPT_DIR ${OUTPUT_DIR}/$(($i + 6)) 17 HIP &
export HIP_VISIBLE_DEVICES=7; srun -n 1 -o $OUTPUT_DIR/output/job.%J.out python $SCRIPT_DIR ${OUTPUT_DIR}/$(($i + 7)) 17 HIP &


echo "Job steps submitted..."
sleep 1
squeue -u `id -un` -s

wait

echo "All Steps completed."

