"""
Concatenate the MSD for all beads in each compartment across all replicas,
and save the average MSD per compartment. One single .h5 is generated with
the average MSD per compartment for each condition.
"""

import os
from functools import reduce

import h5py
import numpy as np


def read_sequence(filename):
    """Extract subcompartment information from file."""
    sequence = np.loadtxt(filename, dtype=str)

    return sequence[:, 1].astype(str)


CHROMOSOME = 1
FOLDER_NUM = 4 if CHROMOSOME == 1 else 5
OUTPUT_FOLDER = f"/scratch/mm146/Polarization/{FOLDER_NUM}_IMR-90_hg38_chr{CHROMOSOME}/7_analysis/msd/data"
NUM_REPLICAS = 32

conditions = [
    "complete",
    "nucleolus",
    # "lamina",
    # "control",
]

chr_sequence = read_sequence(
    F"../../1_inputs/chr{CHROMOSOME}_subcompartments.txt"
)

compartment_indices = {}
for subcompartment in ["A1", "A2", "B1", "B2", "B3", "NA"]:
    compartment_indices[subcompartment] = np.where(
        chr_sequence == subcompartment
    )[0]

compartment_indices["A"] = reduce(
    np.union1d,
    (
        compartment_indices["A1"],
        compartment_indices["A2"],
    ),
)
compartment_indices["B"] = reduce(
    np.union1d,
    (
        compartment_indices["B1"],
        compartment_indices["B2"],
        compartment_indices["B3"],
    ),
)

with h5py.File(
    os.path.join(OUTPUT_FOLDER, "msd-average-per-compartment.h5"), "w"
) as f_out:
    for condition in conditions:
        condition_group = f_out.create_group(condition)

        print(
            f"Concatenating the MSD for condition: {condition}",
            flush=True,
        )

        msd_all_replicas = []
        for replica in range(1, NUM_REPLICAS + 1):
            print("\tProcessing replica:", replica, flush=True)
            with h5py.File(
                os.path.join(
                    OUTPUT_FOLDER,
                    condition,
                    f"msd-all-beads-replica{replica}.h5",
                ),
                "r",
            ) as f:
                msd_all_replicas.append(f["msd-per-bead"][()])

        msd_all_replicas = np.array(msd_all_replicas)

        for label, indices in compartment_indices.items():
            msd_compartment = (
                msd_all_replicas[:, :, indices]
                .transpose(1, 0, 2)
                .reshape(10000, -1)
                .mean(axis=1)
            )

            condition_group.create_dataset(
                label,
                data=msd_compartment,
            )

print("Done!", flush=True)
