"""
Script to compute the pairwise distance distributions for chromosome 1,
for a select pair loci.
"""

import os
import sys

import h5py
import numpy as np
from chroma import structure
from OpenMiChroM.CndbTools import cndbTools

sys.path.append(
    "/scratch/mm146/Polarization/polarization/2_simulations"
)

import nuclear_bodies as nb

condition = sys.argv[1]
replica = int(sys.argv[2])

bead_i = int(sys.argv[3])  # this should be 0-indexed
bead_j = int(sys.argv[4])  # this should be 0-indexed

TRAJ_FOLDER = (
    "/scratch/mm146/Polarization/3_IMR-90_hg38_fields-nb/3_sampling"
)
OUTPUT_FOLDER = "/scratch/mm146/Polarization/3_IMR-90_hg38_fields-nb/4_analysis/distances-distributions/data"

NUCLEUS_RADIUS = 32.5
NUCLEOLI_RADIUS = (NUCLEUS_RADIUS) / 5 ** (1 / 3)
CHROMATIN_DENSITY = 0.30

if condition not in [
    "control",
    "complete",
    "lamina",
    "nucleolus",
]:
    raise ValueError(
        "Invalid condition. Options are: 'control', 'complete', 'lamina', "
        "'nucleolus'."
    )

print(
    "Computing the distances distributions for chromosome 1",
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

chr_sequence = np.array([x.decode("utf-8")[0] for x in traj.ChromSeq])
num_beads = len(chr_sequence)

print(
    f"Extracting positions for beads {bead_i} ({chr_sequence[bead_i]}) and {bead_j} ({chr_sequence[bead_j]})...",
    flush=True,
)

positions = traj.xyz(
    frames=range(0, 10_000, 1), beadSelection=[bead_i, bead_j]
)

distances = np.linalg.norm(
    np.diff(positions, axis=1), axis=2
).reshape(-1)

half_angle_conic_confinement = nb.calc_half_angle(
    density=CHROMATIN_DENSITY,
    nucleus_radius=NUCLEUS_RADIUS,
    nucleolus_radius=NUCLEOLI_RADIUS,
    num_chromatin_beads=num_beads,
)
max_distance = (
    2 * NUCLEUS_RADIUS * np.sin(half_angle_conic_confinement)
)
print(
    f"Estimated maximum distance between beads: {max_distance:.2f}",
    flush=True,
)

distance_bins = np.linspace(0, max_distance, 101)

## distance distribution

print("Computing distance distribution...", flush=True)

pairwise_distance_histogram, bin_edges = np.histogram(
    distances,
    bins=distance_bins,
    density=False,
)

os.makedirs(os.path.join(OUTPUT_FOLDER, condition), exist_ok=True)
with h5py.File(
    os.path.join(
        OUTPUT_FOLDER,
        condition,
        f"pairwise-distance-histogram-replica{replica}-beads-{bead_i}-{bead_j}.h5",
    ),
    "w",
) as saving_file:
    saving_file.create_dataset(
        name="bins",
        data=bin_edges,
    )
    saving_file.create_dataset(
        name="counts",
        data=pairwise_distance_histogram,
    )

print("Done!", flush=True)
