import os
import sys

import h5py
import pandas
import numpy as np
from OpenMiChroM.CndbTools import cndbTools

from chroma import structure, energy

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
    "Computing the distances distributions",
    flush=True,
)
print(f"\tcondition: {condition}", flush=True)
print(f"\treplica: {replica}", flush=True)
print("", flush=True)

filepath = os.path.join(
    traj_folder,
    condition,
    str(replica),
    "nucleus_0.cndb"
)

traj = cndbTools()
traj.load(fileName=filepath)

positions = traj.xyz(frames=range(0, 10_000, 1))

num_frames = positions.shape[0]
num_beads = positions.shape[1]

chr_sequence = np.array([x.decode("utf-8") for x in traj.ChromSeq])
chr_sequence_compartment = np.array([x[0] for x in chr_sequence])

i_idx, j_idx = np.triu_indices(num_beads, k=3)

is_A = chr_sequence_compartment == "A"
is_B = chr_sequence_compartment == "B"
is_N = chr_sequence_compartment == "N"


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

# load types parameters
alpha_matrix = pandas.read_csv("../../1_inputs/ff_nucleus-compartment-level.csv")
alpha_matrix.index = alpha_matrix.columns

alpha_ij = np.zeros((positions.shape[1], positions.shape[1]))
for i in range(alpha_ij.shape[0]):
    for j in range(alpha_ij.shape[1]):
        alpha_ij[i, j] = alpha_matrix.loc[
            chr_sequence[i], chr_sequence[j]
        ]

mol_index = 1
if condition in ["lamina", "complete"]:
    compute_lamina = True

    traj_lamina = cndbTools()
    traj_lamina.load(
        fileName=os.path.join(
            traj_folder,
            condition,
            str(replica),
            f"nucleus_{mol_index}.cndb",
        )
    )
    positions_lamina = traj_lamina.xyz(frames=[0])[0]
    mol_index += 1

    print(
        f"Lamina positions shape: {positions_lamina.shape}",
        flush=True,
    )
    alpha_i_lamina = np.array(
        [alpha_matrix.loc[annot, "LM"] for annot in chr_sequence]
    )
else:
    compute_lamina = False
    positions_lamina = None
    alpha_i_lamina = None


# open nucleolus when applicable
if condition in ["complete", "nucleolus"]:
    compute_nucleolus = True

    traj_nucleolus = cndbTools()
    traj_nucleolus.load(
        fileName=os.path.join(
            traj_folder,
            condition,
            str(replica),
            f"nucleus_{mol_index}.cndb",
        )
    )
    positions_nucleolus = traj_nucleolus.xyz(frames=[0])[0]

    print(
        f"Nucleolus positions shape: {positions_nucleolus.shape}",
        flush=True,
    )
    alpha_i_nucleolus = np.array(
        [alpha_matrix.loc[annot, "NC"] for annot in chr_sequence]
    )
else:
    compute_nucleolus = False
    positions_nucleolus = None
    alpha_i_nucleolus = None


## compute energies

energies = {tag: [] for tag in indices.keys()}

if compute_nucleolus:
    energies["nucleolus"] = []
else:
    energies["nucleolus"] = np.zeros(num_frames)

if compute_lamina:
    energies["lamina"] = []
else:
    energies["lamina"] = np.zeros(num_frames)

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
        distance_nucleolus = structure.compute_distances(frame, positions_nucleolus)
            
        energy_nucleolus = (
            alpha_i_nucleolus[:, None]
            * energy._contact_switch(distance_nucleolus)
            * np.heaviside(3 - distance_nucleolus, 0)
        )
        energies["nucleolus"].append(np.sum(energy_nucleolus))

    if compute_lamina:
        distance_lamina = structure.compute_distances(frame, positions_lamina)
            
        energy_lamina = (
            alpha_i_lamina[:, None]
            * energy._contact_switch(distance_lamina)
            * np.heaviside(3 - distance_lamina, 0)
        )
        energies["lamina"].append(np.sum(energy_lamina))

print("Saving energies...", flush=True)

os.makedirs(os.path.join(output_folder, condition), exist_ok=True)

with h5py.File(
    os.path.join(output_folder, condition, f"contact-energies-replica{replica}.h5"),
    "w",
) as saving_file:
    for name, energy_values in energies.items():

        print(f"\t{name}: {len(energy_values)} frames", flush=True)

        saving_file.create_dataset(
            name=name,
            data=np.array(energy_values),
        )

print("Done!", flush=True)
