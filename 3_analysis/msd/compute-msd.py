"""
Calculate the MSD for chromosome 1 for all the conditions and replicas. The MSD is calculated
for each bead separately and saved in a h5 file. The MSD is calculated using the freud library.
The center of mass is removed before calculating the MSD.
"""

import os
import sys

import freud
import h5py
import numpy as np
from chroma import structure
from OpenMiChroM.CndbTools import cndbTools

condition = sys.argv[1]
replica = int(sys.argv[2])

TRAJ_FOLDER = (
    "/scratch/mm146/Polarization/3_IMR-90_hg38_fields-nb/3_sampling"
)
OUTPUT_FOLDER = "/scratch/mm146/Polarization/3_IMR-90_hg38_fields-nb/4_analysis/msd/data"

NUCLEUS_RADIUS = 32.5
NUCLEOLI_RADIUS = (NUCLEUS_RADIUS) / 5 ** (1 / 3)

if condition not in [
    "control",
    "complete",
    "lamina",
    "nucleolus",
]:
    raise ValueError(
        "Invalid condition. Options are: 'control', 'complete', 'lamina', "
        "'nucleolus', 'isolation', and 'nuclear-bodies'."
    )

print(
    "Computing the MSD for chromosome 1",
    flush=True,
)
print(f"\tcondititon: {condition}", flush=True)
print(f"\treplica: {replica}", flush=True)
print("", flush=True)

filepath = os.path.join(
    TRAJ_FOLDER,
    condition,
    str(replica),
    "nucleus_0.cndb",
)

traj = cndbTools()
traj.load(fileName=filepath)

positions = traj.xyz(frames=range(0, 10_000, 1))
positions -= positions.mean(axis=1, keepdims=True)

msd = freud.msd.MSD()
msd.compute(positions=positions)

os.makedirs(os.path.join(OUTPUT_FOLDER, condition), exist_ok=True)
with h5py.File(
    os.path.join(
        OUTPUT_FOLDER, condition, f"msd-all-beads-replica{replica}.h5"
    ),
    "w",
) as f:
    f.create_dataset("msd-per-bead", data=msd.particle_msd)

print("Done!", flush=True)
