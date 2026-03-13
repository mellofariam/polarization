import os
import sys

import freud
import h5py
import numpy as np
from chroma import structure
from OpenMiChroM.CndbTools import cndbTools

num_chr = sys.argv[1]
condition = sys.argv[2]
replica = int(sys.argv[3])
chromosome = int(sys.argv[4])

if num_chr == "one":
    FOLDER_PATH =  "/scratch/mm146/Polarization/3_IMR-90_hg38_fields-nb/3_single-chr-simulation"
elif num_chr == "two":
    FOLDER_PATH =  "/scratch/mm146/Polarization/3_IMR-90_hg38_fields-nb/3_single-chr-simulation"
else:
    raise ValueError("Invalid number of chromosomes!")

if condition not in [
    "control",
    "complete",
    "lamina",
    "nucleolus",
    "isolation",
    "nuclear-bodies",
]:
    raise ValueError(
        "Invalid condition. Options are: 'control',"
        "'complete', 'lamina', 'nucleolus', "
        "'isolation', and 'nuclear-bodies'."
    )

print(
    f"Computing the MSD for chromosome: {chromosome}",
    flush=True,
)
print(f"\tnumber of chromosomes: {num_chr}", flush=True)
print(f"\tcondititon: {condition}", flush=True)
print(f"\treplica: {replica}", flush=True)
print("", flush=True)

print(f"chromosome {chromosome}", flush=True)
filepath = os.path.join(
    FOLDER_PATH,
    condition,
    str(replica),
    f"nucleus_{chromosome-10}.cndb",
)

traj = cndbTools()
traj.load(fileName=filepath)

positions = traj.xyz(frames=range(0, 10_000, 1))

msd = freud.msd.MSD()
msd.compute(positions=positions)

os.makedirs(f"data/{num_chr}/{condition}", exist_ok=True)
with h5py.File(
    f"data/{num_chr}/{condition}/msd-all-beads-chr{chromosome}-replica-{replica}.h5",
    "w",
) as f:
    f.create_dataset("msd-per-bead", data=msd.particle_msd)

print("Done!", flush=True)
