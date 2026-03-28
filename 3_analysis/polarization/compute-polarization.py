"""
Script to calculate the compartment polarization of the simulated chromosome.
"""

import os
import sys

import h5py
import numpy as np
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
    "Computing the polarization",
    flush=True,
)
print(f"\tcondititon: {condition}", flush=True)
print(f"\treplica: {replica}", flush=True)
print("", flush=True)

os.makedirs(os.path.join(output_folder, condition), exist_ok=True)
with h5py.File(
    os.path.join(
        output_folder, condition, f"polarization-replica{replica}.h5"
    ),
    "w",
) as saving_file:

    filepath = os.path.join(
        traj_folder,
        condition,
        str(replica),
        "nucleus_0.cndb",
    )

    traj = cndbTools()
    traj.load(fileName=filepath)

    subcompartments = list(
        map(lambda x: str(x, "utf-8"), traj.ChromSeq)
    )

    positions = traj.xyz(frames=range(0, 10_000, 1))
    center_of_mass = positions.mean(axis=1)

    positions -= positions.mean(axis=1, keepdims=True)

    q = {"A": 1, "B": -1}

    q_mask = np.array(
        [q.get((subcmpt[0]), 0) for subcmpt in subcompartments]
    )

    dipole_moment = np.sum(
        positions * q_mask.reshape((1, -1, 1)), axis=1
    )  # sum along beads -> shape result Nframes x 3
    dipole_magnitude = np.linalg.norm(
        dipole_moment, axis=1
    )  # norm along xyz axis (1)

    # Normalize each vector in both groups
    dipole_norm = (
        dipole_moment
        / np.linalg.norm(dipole_moment, axis=1)[:, np.newaxis]
    )
    com_norm = (
        center_of_mass
        / np.linalg.norm(center_of_mass, axis=1)[:, np.newaxis]
    )

    # Compute cosine of the angle between each pair
    cos_theta = np.einsum("ij,ij->i", dipole_norm, com_norm)

    # Numerical stability: clip values to [-1, 1]
    cos_theta = np.clip(cos_theta, -1.0, 1.0)

    # Convert to angles
    dipole_alignment = np.arccos(cos_theta)

    print("Done! Saving in file...", flush=True)

    saving_file.create_dataset(
        name="dipole_moment",
        data=dipole_moment,
    )
    saving_file.create_dataset(
        name="center_of_mass",
        data=center_of_mass,
    )
    saving_file.create_dataset(
        name="dipole_magnitude",
        data=dipole_magnitude,
    )
    saving_file.create_dataset(
        name="dipole_alignment",
        data=dipole_alignment,
    )
    print("Done!", flush=True)
