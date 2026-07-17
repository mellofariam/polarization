"""Load CNDB trajectory and compute polar histogram."""

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
    "Computing localization of compartments in the nucleus",
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


def spherical_coordinates(xyz):
    if xyz.ndim == 2:
        xyz = xyz[None, :, :]

    r = np.sqrt(np.sum(xyz**2, axis=2))
    theta = np.arccos(xyz[:, :, 2] / r)
    phi = np.arctan2(xyz[:, :, 1], xyz[:, :, 0])
    return r, theta, phi


def align_vector_to_z_axis(xyz, vectors):
    """Rotate each frame so the provided vector points along +z."""
    if xyz.ndim == 2:
        xyz = xyz[None, :, :]
    if vectors.ndim == 1:
        vectors = vectors[None, :]

    z_axis = np.array([0.0, 0.0, 1.0])
    rotated = np.empty_like(xyz)

    for frame_idx, (frame_xyz, vector) in enumerate(zip(xyz, vectors)):
        norm = np.linalg.norm(vector)
        if norm == 0:
            rotated[frame_idx] = frame_xyz
            continue

        unit_vector = vector / norm
        if np.allclose(unit_vector, z_axis):
            rotated[frame_idx] = frame_xyz
            continue

        if np.allclose(unit_vector, -z_axis):
            rotation_matrix = np.array(
                [
                    [-1.0, 0.0, 0.0],
                    [0.0, -1.0, 0.0],
                    [0.0, 0.0, 1.0],
                ]
            )
        else:
            rotation_axis = np.cross(unit_vector, z_axis)
            rotation_axis /= np.linalg.norm(rotation_axis)
            angle = np.arccos(np.clip(np.dot(unit_vector, z_axis), -1.0, 1.0))

            skew_axis = np.array(
                [
                    [0.0, -rotation_axis[2], rotation_axis[1]],
                    [rotation_axis[2], 0.0, -rotation_axis[0]],
                    [-rotation_axis[1], rotation_axis[0], 0.0],
                ]
            )
            rotation_matrix = (
                np.eye(3)
                + np.sin(angle) * skew_axis
                + (1.0 - np.cos(angle)) * (skew_axis @ skew_axis)
            )

        rotated[frame_idx] = frame_xyz @ rotation_matrix.T

    return rotated


com = np.mean(positions, axis=1)
positions = align_vector_to_z_axis(positions, com)

r, theta, phi = spherical_coordinates(positions)

theta_bins = np.linspace(0, np.pi / 4, 21)
r_bins = np.linspace(NUCLEOLI_RADIUS + 0.5, NUCLEUS_RADIUS - 0.5, 11)
# r_bins = np.arange(20, 32.5, 1.0)

counts = {}
for annotation in ["A1", "A2", "B1", "B2", "B3", "NA"]:
    is_annotation = chr_sequence == annotation
    r_annotation = r[:, is_annotation]
    theta_annotation = theta[:, is_annotation]

    counts[annotation], _, _ = np.histogram2d(
        r_annotation.flatten(),
        theta_annotation.flatten(),
        bins=[r_bins, theta_bins],
    )

os.makedirs(os.path.join(output_folder, condition), exist_ok=True)
with h5py.File(
    os.path.join(
        output_folder,
        condition,
        f"localization-replica{replica}.h5",
    ),
    "w",
) as f:
    for annotation, hist in counts.items():
        f.create_dataset(f"{annotation}", data=hist)
    f.create_dataset("r_bins", data=r_bins)
    f.create_dataset("theta_bins", data=theta_bins)
