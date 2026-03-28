"""
Script to compute the pairwise distance distributions for chromosome 1,
for each pair of subcompartments. We also compute the mean and standard
deviation of the distance matrix across frames.
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
traj_folder = sys.argv[3]
output_folder = sys.argv[4]

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
    traj_folder,
    condition,
    str(replica),
    "nucleus_0.cndb",
)

traj = cndbTools()
traj.load(fileName=filepath)

positions = traj.xyz(frames=range(0, 10_000, 10))

num_frames = positions.shape[0]
num_beads = positions.shape[1]

chr_sequence = np.array([x.decode("utf-8")[0] for x in traj.ChromSeq])
i_idx, j_idx = np.triu_indices(num_beads, k=3)

is_A = chr_sequence == "A"
is_B = chr_sequence == "B"
is_N = chr_sequence == "N"


def get_pairs(m1, m2):
    "Logic: (m1[i] AND m2[j]) OR (m1[j] AND m2[i])"

    mask = (m1[i_idx] & m2[j_idx]) | (m1[j_idx] & m2[i_idx])
    return np.column_stack((i_idx[mask], j_idx[mask]))


indices = {
    "AA": get_pairs(is_A, is_A),
    "AB": get_pairs(is_A, is_B),
    "BB": get_pairs(is_B, is_B),
    "AN": get_pairs(is_A, is_N),
    "BN": get_pairs(is_B, is_N),
    "NN": get_pairs(is_N, is_N),
}

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
genomic_distance_bins = np.logspace(0, 5, 101)

## distance distribution

print("Computing distance distribution...", flush=True)

distance_histograms = {tag: np.zeros((100, 100)) for tag in indices}
distance_histograms["all"] = np.zeros((100, 100))

sum_distance_matrix = np.zeros((num_beads, num_beads))
sum_sq_distance_matrix = np.zeros((num_beads, num_beads))
total_frames = 0

CHUNK_SIZE = 25
for start_frame in range(0, num_frames, CHUNK_SIZE):
    end_frame = min(start_frame + CHUNK_SIZE, num_frames)

    # Compute distances for just this chunk (Shape: chunk, N, N)
    chunk_positions = positions[start_frame:end_frame]
    chunk_distances = structure.compute_distances_trajectory(
        chunk_positions, chunk_positions
    )

    # 1. Update "all" histogram
    distance_values_all = chunk_distances[:, i_idx, j_idx]
    genomic_distances_all = np.abs(i_idx - j_idx)
    h, xedges, yedges = np.histogram2d(
        distance_values_all.flatten(),
        np.broadcast_to(
            genomic_distances_all, distance_values_all.shape
        ).flatten(),
        bins=[distance_bins, genomic_distance_bins],
        density=False,
    )
    distance_histograms["all"] += h

    # 2. Update subcompartment histograms
    for tag, pairs in indices.items():
        if len(pairs) == 0:
            continue

        distance_values = chunk_distances[:, pairs[:, 0], pairs[:, 1]]
        genomic_distances = np.abs(pairs[:, 0] - pairs[:, 1])
        h, _, _ = np.histogram2d(
            distance_values.flatten(),
            np.broadcast_to(
                genomic_distances, distance_values.shape
            ).flatten(),
            bins=[distance_bins, genomic_distance_bins],
            density=False,
        )
        distance_histograms[tag] += h

    # 3. Update distance matrix sums for mean/std calculation
    sum_distance_matrix += np.sum(chunk_distances, axis=0)
    sum_sq_distance_matrix += np.sum(chunk_distances**2, axis=0)
    total_frames += chunk_distances.shape[0]

    print(
        f"  Processed frames {start_frame} to {end_frame}", flush=True
    )

# The result in distance_histograms is now the sum of all frames.

mean_distance_matrix = sum_distance_matrix / total_frames
mean_sq_distance_matrix = sum_sq_distance_matrix / total_frames

os.makedirs(os.path.join(output_folder, condition), exist_ok=True)
with h5py.File(
    os.path.join(
        output_folder,
        condition,
        f"distance-histograms-replica{replica}.h5",
    ),
    "w",
) as saving_file:
    saving_file.create_dataset(
        name="xbins",
        data=xedges,
    )
    saving_file.create_dataset(
        name="ybins",
        data=yedges,
    )
    for tag, distances in distance_histograms.items():
        saving_file.create_dataset(
            name=tag,
            data=distances,
        )
    saving_file.attrs.create("nframes", data=len(positions))

with h5py.File(
    os.path.join(
        output_folder,
        condition,
        f"distance-matrix-replica{replica}.h5",
    ),
    "w",
) as saving_file:
    saving_file.create_dataset(
        name="mean_distance_matrix",
        data=mean_distance_matrix,
    )
    saving_file.create_dataset(
        name="mean_sq_distance_matrix",
        data=mean_sq_distance_matrix,
    )
    saving_file.attrs.create("nframes", data=len(positions))

print("Done!", flush=True)
