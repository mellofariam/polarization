import os
import sys
import time
from datetime import timedelta

import h5py
import mdtraj
import numpy as np
import pandas
from chroma import energy, structure
from OpenMiChroM.CndbTools import cndbTools

condition = sys.argv[1]
replica = int(sys.argv[2])
traj_folder = sys.argv[3]
output_folder = sys.argv[4]

NUCLEUS_RADIUS = 32.5
NUCLEOLI_RADIUS = (NUCLEUS_RADIUS) / 5 ** (1 / 3)

if condition not in [
    "complete",
    "nucleolus",
    "lamina",
    "control",
]:
    raise ValueError(
        f"Invalid condition: {condition}. Options are: 'complete', 'nucleolus', 'lamina', and 'control'."
    )

print(
    "Computing the SASA",
    flush=True,
)
print(f"\tcondition: {condition}", flush=True)
print(f"\treplica: {replica}", flush=True)
print("", flush=True)

print("Opening trajectory...", flush=True)

filepath = os.path.join(
    traj_folder, condition, str(replica), "nucleus_0.cndb"
)

cndb = cndbTools()
cndb.load(fileName=filepath)

positions = cndb.xyz(frames=range(0, 10_000, 1))

num_frames = positions.shape[0]
num_beads = positions.shape[1]

chr_sequence = np.array([x.decode("utf-8") for x in cndb.ChromSeq])
chr_sequence_compartment = np.array([x[0] for x in chr_sequence])

is_A = chr_sequence_compartment == "A"
is_B = chr_sequence_compartment == "B"

topfile = os.path.join(
    traj_folder, condition, str(replica), "after-anneal_0.pdb"
)

traj = mdtraj.load(topfile)
traj.xyz = positions
traj.time = np.arange(0, 0.002 * len(positions), 0.002)

start = time.perf_counter()
print("Calculating SASA...", flush=True)
sasa = mdtraj.shrake_rupley(
    traj,
    probe_radius=0.5,
    change_radii={
        "C": 0.5,
        "CA": 0.5,
    },
    mode="atom",
)
print(f"Done in {timedelta(seconds=time.perf_counter() - start)}.")

os.makedirs(os.path.join(output_folder, condition), exist_ok=True)

print("Saving files...", flush=True)
with h5py.File(
    os.path.join(
        output_folder, condition, f"sasa-replica{replica}.h5"
    ),
    "w",
) as saving_file:
    saving_file.create_dataset(
        name="all",
        data=sasa,
    )
    saving_file.create_dataset(
        name="A", data=sasa[:, is_A].sum(axis=1)
    )
    saving_file.create_dataset(
        name="B", data=sasa[:, is_B].sum(axis=1)
    )

print("Done!", flush=True)
