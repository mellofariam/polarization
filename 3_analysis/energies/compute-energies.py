import os
import sys

import h5py
import pandas
import numpy as np
from OpenMiChroM.CndbTools import cndbTools

from chroma import structure, energy

condition = sys.argv[1]
replica = int(sys.argv[2])

TRAJ_FOLDER = "/scratch/mm146/Polarization/3_IMR-90_hg38_fields-nb/3_sampling"
OUTPUT_FOLDER = "/scratch/mm146/Polarization/3_IMR-90_hg38_fields-nb/4_analysis/energies/data"

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
    f"Computing the distances distributions for chromosome 1",
    flush=True,
)
print(f"\tcondititon: {condition}", flush=True)
print(f"\treplica: {replica}", flush=True)
print("", flush=True)

filepath = os.path.join(
    TRAJ_FOLDER,
    condition,
    str(replica),
    "nucleus_0.cndb"
)

traj = cndbTools()
traj.load(fileName=filepath)

positions = traj.xyz(frames=range(0, 10_000, 1))

n_frames = positions.shape[0]

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

# load types parameters
alpha_matrix = pandas.read_csv("../../1_inputs/ff_compartments-and-nb.csv")
alpha_matrix.index = alpha_matrix.columns

alpha_ij = np.zeros((positions.shape[1], positions.shape[1]))
for i in range(alpha_ij.shape[0]):
    for j in range(alpha_ij.shape[1]):
        alpha_ij[i, j] = alpha_matrix.loc[
            chr_sequence[i], chr_sequence[j]
        ]

# open nucleolus when applicable
if condition in ["complete", "nucleolus"]:
    compute_nucleolus = True

    nucleolus_path = os.path.join(
        TRAJ_FOLDER,
        condition,
        str(replica),
        f"nucleus_1.cndb",
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
    alpha_i_nucleolus = np.array(
        [alpha_matrix.loc[annot, "NC"] for annot in chr_sequence]
    )
else:
    compute_nucleolus = False
    nucleolus_positions = None
    alpha_i_nucleolus = None

if condition in ["lamina", "complete"]:
    compute_lamina = True
else:
    compute_lamina = False


## compute energies

energies = {tag: [] for tag in indices.keys()}

if compute_nucleolus:
    energies["nucleolus"] = []

print("Computing energies...", flush=True)
for idx, frame in enumerate(positions):
    if idx % 100 == 0:
        print(f"\tAnalyzing frame {idx}", flush=True)

    distance_matrix = structure.compute_distances(frame, frame)

    energy_matrix = (
        alpha_ij
        * energy._contact_switch(distance_matrix)
        * np.heaviside(3 - distance_matrix, 0)
    )

    for tag, pairs in indices.items():
        energies[tag].append(
            np.sum(energy_matrix[pairs[:, 0], pairs[:, 1]])
        )

    if compute_nucleolus:
        distance_nucleolus = (
            structure.compute_distances(frame, nucleolus_positions)
            - NUCLEOLI_RADIUS
        )
        energy_nucleolus = (
            alpha_i_nucleolus[:, None]
            * energy._contact_switch(distance_nucleolus)
            * np.heaviside(3 - distance_nucleolus, 0)
        )
        energies["nucleolus"].append(np.sum(energy_nucleolus))

if compute_lamina:
    alpha_lamina = {
        "B1": -1.034611,
        "B2": -0.928730,
        "B3": -0.741693,
        "NA": -0.667915,
    }
    alpha_i_lamina = np.array(
        [alpha_lamina.get(annot, 0) for annot in chr_sequence]
    )
    distance_lamina = NUCLEUS_RADIUS - np.linalg.norm(
        positions, axis=2
    )
    energy_lamina = (
        alpha_i_lamina[None, :]
        * energy._contact_switch(distance_lamina)
        * np.heaviside(3 - distance_lamina, 0)
    )
    energies["lamina"] = np.sum(energy_lamina, axis=1)

if "nucleolus" not in energies:
    energies["nucleolus"] = np.zeros(n_frames)
if "lamina" not in energies:
    energies["lamina"] = np.zeros(n_frames)

print("Saving energies...", flush=True)

os.makedirs(os.path.join(OUTPUT_FOLDER, condition), exist_ok=True)

with h5py.File(
    os.path.join(OUTPUT_FOLDER, condition, f"contact-energies-replica{replica}.h5"),
    "w",
) as saving_file:
    for name, energy_values in energies.items():

        print(f"\t{name}: {len(energy_values)} frames", flush=True)

        saving_file.create_dataset(
            name=name,
            data=np.array(energy_values),
        )

print("Done!", flush=True)
