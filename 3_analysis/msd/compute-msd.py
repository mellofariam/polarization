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
traj_folder = sys.argv[3]
output_folder = sys.argv[4]

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
    traj_folder,
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

os.makedirs(os.path.join(output_folder, condition), exist_ok=True)
with h5py.File(
    os.path.join(
        output_folder, condition, f"msd-all-beads-replica{replica}.h5"
    ),
    "w",
) as f:
    f.create_dataset("msd-per-bead", data=msd.particle_msd)

print("Done!", flush=True)
