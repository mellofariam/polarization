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

NUCLEUS_RADIUS = 28
NUCLEOLI_RADIUS = 16.4
SPECKLES_RADIUS = 1.5
SPECKLES_NUMBER = 40

if num_chr == "one":
    FOLDER_PATH = "/scratch/mm146/Polarization/3_IMR-90_hg38_fields-nb/3_single-chr-simulation"
elif num_chr == "two":
    FOLDER_PATH = "/scratch/mm146/Polarization/3_IMR-90_hg38_fields-nb/3_single-chr-simulation"
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
        "Invalid condition. Options are: 'control', 'complete', 'lamina', "
        "'nucleolus', 'isolation', and 'nuclear-bodies'."
    )

print(
    f"Computing the distances distributions for chromosome: {chromosome}",
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

positions = traj.xyz(frames=range(0, 10_000, 10))

chr_sequence = traj.ChromSeq
chr_sequence = np.array(
    list(map(lambda x: x.decode("utf-8"), chr_sequence))
)

indices = {
    "AA": [],
    "AB": [],
    "BB": [],
    "AN": [],
    "BN": [],
    "NN": [],
}

for i, annot1 in enumerate(chr_sequence):
    for j, annot2 in enumerate(chr_sequence):
        if abs(i - j) > 2:
            tag = f"{annot1[0]}{annot2[0]}"

            if tag in indices:
                indices[tag].append((i, j))

for tag in indices:
    indices[tag] = np.array(indices[tag])

os.makedirs(f"data/{num_chr}/{condition}", exist_ok=True)

## distance distribution

print("Computing distance distribution...", flush=True)
distance_matrix = structure.compute_distances_trajectory(
    positions, positions
)

histogram_bins = np.linspace(0, 32, 100)
distance_histograms = {}

i, j = np.triu_indices(positions.shape[1], k=3)

distance_histograms["all"], _ = np.histogram(
    distance_matrix[:, i, j],
    bins=histogram_bins,
)
for tag, pairs in indices.items():
    distance_histograms[tag], _ = np.histogram(
        distance_matrix[:, pairs[:, 0], pairs[:, 1]],
        bins=histogram_bins,
        density=False,
    )

with h5py.File(
    f"data/{num_chr}/{condition}/distance-histograms-chr{chromosome}-replica{replica}.h5",
    "w",
) as saving_file:
    saving_file.create_dataset(
        name="bins",
        data=histogram_bins,
    )
    for tag in distance_histograms:
        saving_file.create_dataset(
            name=f"{tag}",
            data=distance_histograms[tag],
        )
    saving_file.attrs.create("nframes", data=1000)

with h5py.File(
    f"data/{num_chr}/{condition}/distance-matrix-with-std-chr{chromosome}-replica{replica}.h5",
    "w",
) as saving_file:
    saving_file.create_dataset(
        name="distance_matrix_mean",
        data=np.mean(distance_matrix, axis=0),
    )
    saving_file.create_dataset(
        name="distance_matrix_std",
        data=np.std(distance_matrix, axis=0, ddof=1),
    )
    saving_file.attrs.create("nframes", data=1000)

print("Done!", flush=True)