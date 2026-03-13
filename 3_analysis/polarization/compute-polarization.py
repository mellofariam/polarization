import os
import sys

import h5py
import numpy as np
from OpenMiChroM.CndbTools import cndbTools

num_chr = sys.argv[1]
condition = sys.argv[2]
replica = int(sys.argv[3])

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
        "'nucleolus', 'isoaltion', and 'nuclear-bodies'."
    )

if num_chr == "one":
    chromosomes = [10]
elif num_chr == "two":
    chromosomes = [10, 11]
else:
    raise ValueError("Invalid number of chromosomes!")

print(
    f"Computing the polarization for chromosomes: {chromosomes}",
    flush=True,
)
print(f"\tnumber of chromosomes: {num_chr}", flush=True)
print(f"\tcondititon: {condition}", flush=True)
print(f"\treplica: {replica}", flush=True)
print("", flush=True)

os.makedirs(f"data/{num_chr}/{condition}", exist_ok=True)
with h5py.File(
    f"data/{num_chr}/{condition}/polarization-replica{replica}.h5",
    "w",
) as saving_file:

    two_chr_positions = []
    two_chr_sequence = []
    for chromosome in chromosomes:
        print(f"chromosome {chromosome}", flush=True)
        filepath = os.path.join(
            FOLDER_PATH,
            FOLDER_PATH,
            condition,
            str(replica),
            f"nucleus_{chromosome-10}.cndb",
        )

        traj = cndbTools()
        traj.load(fileName=filepath)

        subcompartments = list(
            map(lambda x: str(x, "utf-8"), traj.ChromSeq)
        )

        positions = traj.xyz(frames=range(0, 10_000, 1))
        center_of_mass = positions.mean(axis=1)
        if num_chr == "two":
            two_chr_positions.append(positions)
            two_chr_sequence.extend(subcompartments)
        positions -= positions.mean(axis=1, keepdims=True)

        q = {"A": 1, "B": -1}

        q_mask = np.array(
            [q.get((subcmpt[0]), 0) for subcmpt in subcompartments]
        )

        dipole_moment = np.sum(
            positions * q_mask.reshape(1, -1, 1), axis=1
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
        saving_file.create_group(f"chr{chromosome}")
        saving_file.create_dataset(
            name=f"chr{chromosome}/dipole_moment", data=dipole_moment
        )
        saving_file.create_dataset(
            name=f"chr{chromosome}/center_of_mass", data=center_of_mass
        )
        saving_file.create_dataset(
            name=f"chr{chromosome}/dipole_magnitude", data=dipole_magnitude
        )
        saving_file.create_dataset(
            name=f"chr{chromosome}/dipole_alignment", data=dipole_alignment
        )
        print("Done!", flush=True)

    if two_chr_positions:
        print("both chromosomes", flush=True)

        positions = np.concatenate(two_chr_positions, axis=1)
        center_of_mass = positions.mean(axis=1)
        positions -= positions.mean(axis=1, keepdims=True)

        subcompartments = two_chr_sequence

        q = {"A": 1, "B": -1}

        q_mask = np.array(
            [q.get((subcmpt[0]), 0) for subcmpt in subcompartments]
        )

        dipole_moment = np.sum(
            positions * q_mask.reshape(1, -1, 1), axis=1
        )  # sum along beads -> shape result Nframes x 3
        dipole_magnitude = np.linalg.norm(
            dipole_moment, axis=1
        )  # norm along xyz axis (1)

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
        saving_file.create_group("both")
        saving_file.create_dataset(
            name="both/dipole_moment", data=dipole_moment
        )
        saving_file.create_dataset(
            name="both/center_of_mass", data=center_of_mass
        )
        saving_file.create_dataset(
            name="both/dipole_magnitude", data=dipole_magnitude
        )
        saving_file.create_dataset(
            name="both/dipole_alignment", data=dipole_alignment
        )
        print("Done!", flush=True)
