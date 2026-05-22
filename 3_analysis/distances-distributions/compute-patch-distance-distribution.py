"""
Script to compute the minimum distance distributions for chromosome 1,
for selected pairs of compartment patches.
"""

import itertools
import os
import sys

import h5py
import numpy as np
from OpenMiChroM.CndbTools import cndbTools
from chroma import structure

condition = sys.argv[1]
replica = int(sys.argv[2])
traj_folder = sys.argv[3]
output_folder = sys.argv[4]
tag = sys.argv[5]

selected_patches = sorted(set(map(int, sys.argv[6:])))


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

if len(selected_patches) < 2:
    raise ValueError("Please provide at least two patch ids.")

print(
    "Computing the distances distributions for chromosome 1",
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

chr_sequence = np.array([x.decode("utf-8")[0] for x in traj.ChromSeq])
num_beads = len(chr_sequence)

patches = []
start = 0
current = chr_sequence[0]

for i in range(1, len(chr_sequence)):
    if chr_sequence[i] != current:
        patches.append(
            {
                "patch": len(patches),
                "compartment": current,
                "start": start,
                "end": i - 1,
                "length": i - start,
            }
        )
        start = i
        current = chr_sequence[i]

patches.append(
    {
        "patch": len(patches),
        "compartment": current,
        "start": start,
        "end": len(chr_sequence) - 1,
        "length": len(chr_sequence) - start,
    }
)

num_patches = len(patches)
invalid_patches = [
    patch for patch in selected_patches if patch < 0 or patch >= num_patches
]
if invalid_patches:
    raise ValueError(f"Invalid patch ids: {invalid_patches}")

selected_pairs = list(itertools.combinations(selected_patches, 2))

print(f"\tnumber of beads: {num_beads}", flush=True)
print(f"\tnumber of patches: {num_patches}", flush=True)
print(f"\tselected patches: {selected_patches}", flush=True)
print("", flush=True)

print("Extracting positions...", flush=True)
positions = traj.xyz(frames=range(0, 10_000, 1))

max_distance = 10
# max_distance = (2 * NUCLEUS_RADIUS)

distance_bins = np.linspace(0, max_distance, 101)

os.makedirs(os.path.join(output_folder, condition), exist_ok=True)
with h5py.File(
    os.path.join(
        output_folder,
        condition,
        f"patch-distance-histogram-replica{replica}-{tag}.h5",
    ),
    "w",
) as saving_file:
    saving_file.create_dataset(
        name="bins",
        data=distance_bins,
    )

    for patch_i, patch_j in selected_pairs:
        patch1 = patches[patch_i]
        patch2 = patches[patch_j]

        print(
            (
                f"Processing patches {patch_i}-{patch_j}: "
                f"{patch1['compartment']} [{patch1['start']}, {patch1['end']}] "
                f"and "
                f"{patch2['compartment']} [{patch2['start']}, {patch2['end']}]"
            ),
            flush=True,
        )

        patch1_positions = positions[
            :, patch1["start"] : patch1["end"] + 1, :
        ]
        patch2_positions = positions[
            :, patch2["start"] : patch2["end"] + 1, :
        ]

        distances = structure.compute_distances_trajectory(
            patch1_positions, patch2_positions
        )
        min_distances = np.min(distances, axis=(1, 2))

        pairwise_distance_histogram, _ = np.histogram(
            min_distances,
            bins=distance_bins,
            density=False,
        )

        group = saving_file.create_group(f"{patch_i}-{patch_j}")
        group.create_dataset(
            name="counts",
            data=pairwise_distance_histogram,
        )
        group.create_dataset(
            name="distances",
            data=min_distances,
        )

print("Done!", flush=True)
