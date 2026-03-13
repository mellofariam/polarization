import os
import sys

import h5py
import freud
import numpy as np
from OpenMiChroM.CndbTools import cndbTools

sys.path.insert(0, "/home/mm146")

from chroma import structure

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
        "Invalid condition. Options are: 'control', 'complete', 'lamina', 'nucleolus', 'isolation', and 'nuclear-bodies'."
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

positions = traj.xyz(frames=range(0, 10_000, 1))

os.makedirs(f"data/{num_chr}/{condition}", exist_ok=True)

## distance to the center of mass

print("Computing distances to CoM...", end=" ", flush=True)
distance_CoM = np.linalg.norm(
    positions - np.mean(positions, axis=1, keepdims=True), axis=2
)
with h5py.File(
    f"data/{num_chr}/{condition}/distance-CoM-chr{chromosome}-replica{replica}.h5",
    "w",
) as saving_file:
    saving_file.create_dataset(
        name="average-distance-CoM",
        data=np.mean(distance_CoM, axis=0),
    )
    dataset = saving_file.create_dataset(
        name="distance-CoM", data=distance_CoM[::10, :]
    )
    dataset.attrs.create("nframes", data=1000)

del distance_CoM
print("Done!", flush=True)


## average distance matrix

print("Computing average distance matrix...", flush=True)
distance_matrix = np.zeros((positions.shape[1], positions.shape[1]))
for idx, frame in enumerate(positions):
    if idx % 100 == 0:
        print(f"\tAnalyzing frame {idx}", flush=True)
    distance_matrix += structure.compute_distances(frame, frame)

distance_matrix /= positions.shape[0]

np.savetxt(
    f"data/{num_chr}/{condition}/average-distance-matrix-chr{chromosome}-replica{replica}.dat",
    distance_matrix,
)
print("Done!", flush=True)

del distance_matrix

if condition != "isolation":
    ## distance to the lamina

    print("Computing distances to the lamina...", end=" ", flush=True)

    distance_lamina = NUCLEUS_RADIUS - np.linalg.norm(
        positions, axis=2
    )
    average_distance_lamina = np.mean(
        distance_lamina, axis=0
    )  # on average how much each bead binds to the lamina

    with h5py.File(
        f"data/{num_chr}/{condition}/distance-lamina-chr{chromosome}-replica{replica}.h5",
        "w",
    ) as saving_file:
        saving_file.create_dataset(
            name="average-distance-lamina",
            data=average_distance_lamina,
        )
        dataset = saving_file.create_dataset(
            name="distance-lamina", data=distance_lamina[::10, :]
        )
        dataset.attrs.create("nframes", data=1000)

    del distance_lamina
    del average_distance_lamina

    print("Done!", flush=True)

    ## distance to speckles

    print("Computing distances to speckles...", flush=True)

    speckles_path = os.path.join(
        FOLDER_PATH,
        condition,
        str(replica),
        f"nucleus_{1 if num_chr == 'one' else 2}.cndb",
    )

    speckles = cndbTools()
    speckles.load(fileName=speckles_path)

    speckles_positions = speckles.xyz(frames=[0])[
        0
    ]  # all have the same positions
    print(
        f"Speckles positions shape: {speckles_positions.shape}",
        flush=True,
    )

    min_distance_speckles = []
    for idx, frame in enumerate(positions):
        if idx % 100 == 0:
            print(f"\tAnalyzing frame {idx}", flush=True)
        distance_speckles = structure.compute_distances(
            frame, speckles_positions
        ) - SPECKLES_RADIUS
        min_distance_speckles.append(
            np.min(distance_speckles, axis=1)
        )

    min_distance_speckles = np.asarray(min_distance_speckles)

    with h5py.File(
        f"data/{num_chr}/{condition}/distance-speckles-chr{chromosome}-replica{replica}.h5",
        "w",
    ) as saving_file:
        saving_file.create_dataset(
            name="average-distance-speckles",
            data=np.mean(min_distance_speckles, axis=0),
        )
        dataset = saving_file.create_dataset(
            name="distance-speckles",
            data=min_distance_speckles[::10, :],
        )
        dataset.attrs.create("nframes", data=1000)

    del speckles
    del distance_speckles
    del min_distance_speckles

    print("Done!", flush=True)

    ## distance to nucleolus

    print("Computing distances to nucleolus...", flush=True)

    nucleolus_path = os.path.join(
        FOLDER_PATH,
        condition,
        str(replica),
        f"nucleus_{2 if num_chr == 'one' else 3}.cndb",
    )

    nucleolus = cndbTools()
    nucleolus.load(fileName=nucleolus_path)

    nucleolus_positions = nucleolus.xyz(frames=[0])[
        0
    ]  # all have the same positions
    print(
        f"Nucleolus positions shape: {nucleolus_positions.shape}",
        flush=True,
    )

    min_distance_nucleolus = []
    for idx, frame in enumerate(positions):
        if idx % 100 == 0:
            print(f"\tAnalyzing frame {idx}", flush=True)
        distance_nucleolus = structure.compute_distances(
            frame, nucleolus_positions
        ) - NUCLEOLI_RADIUS
        min_distance_nucleolus.append(
            np.min(distance_nucleolus, axis=1)
        )

    min_distance_nucleolus = np.asarray(min_distance_nucleolus)

    with h5py.File(
        f"data/{num_chr}/{condition}/distance-nucleolus-chr{chromosome}-replica{replica}.h5",
        "w",
    ) as saving_file:
        saving_file.create_dataset(
            name="average-distance-nucleolus",
            data=np.mean(min_distance_nucleolus, axis=0),
        )
        dataset = saving_file.create_dataset(
            name="distance-nucleolus",
            data=min_distance_nucleolus[::10, :],
        )
        dataset.attrs.create("nframes", data=1000)

    del nucleolus
    del distance_nucleolus
    del min_distance_nucleolus

    print("Done!", flush=True)
