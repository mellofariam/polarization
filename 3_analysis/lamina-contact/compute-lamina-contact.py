"""
Calculate the contact of each bead with the lamina.
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
nuclear_body_treatment = sys.argv[5]

NUCLEUS_RADIUS = 32.5
NUCLEOLI_RADIUS = (NUCLEUS_RADIUS) / 5 ** (1 / 3)
CONTACT_THRESHOLD = 2.0

if condition not in [
    "complete",
    "lamina",
]:
    raise ValueError("Invalid condition. Options are: 'complete' and 'lamina'.")

if nuclear_body_treatment not in ["fields", "beads"]:
    raise ValueError(
        "Invalid nuclear body treatment. Options are: 'fields' and 'beads'."
    )

print(
    f"Computing the lamina contacts for {traj_folder.split('/')[-1]}",
    flush=True,
)
print(f"\tcondititon: {condition}", flush=True)
print(f"\treplica: {replica}", flush=True)
print(f"\tnuclear_body_treatment: {nuclear_body_treatment}", flush=True)
print("", flush=True)

filepath = os.path.join(
    traj_folder,
    condition,
    str(replica),
    "nucleus_0.cndb",
)

traj = cndbTools()
traj.load(fileName=filepath)

chr_sequence = [x.decode("utf-8")[0] for x in traj.ChromSeq]
B_indices = np.where(np.array(chr_sequence) == "B")[0]

positions = traj.xyz(frames=range(0, 10_000, 1), beadSelection=B_indices.tolist())

if nuclear_body_treatment == "fields":
    # calculate the distance to the lamina as the distance to the center of the nucleus
    distances_to_lamina = NUCLEUS_RADIUS - np.linalg.norm(positions, axis=2)
elif nuclear_body_treatment == "beads":
    lamina_traj = cndbTools()
    lamina_traj.load(
        os.path.join(
            traj_folder,
            condition,
            str(replica),
            "nucleus_1.cndb",
        )
    )
    lamina_positions = lamina_traj.xyz(frames=[0])[0]

    distances_to_lamina = []
    for frame in positions:
        frame_distances_to_lamina = structure.compute_distances(frame, lamina_positions)
        distances_to_lamina.append(frame_distances_to_lamina.min(axis=1))
    distances_to_lamina = np.array(distances_to_lamina)
else:
    raise ValueError(
        "Invalid nuclear body treatment. Options are: 'fields' and 'beads'."
    )

os.makedirs(os.path.join(output_folder, condition), exist_ok=True)
with h5py.File(
    os.path.join(output_folder, condition, f"lamina-contact-replica{replica}.h5"),
    "w",
) as f:
    f.create_dataset("lamina-contact", data=(distances_to_lamina <= CONTACT_THRESHOLD).astype(int))
    f.create_dataset("B-indices", data=B_indices)
