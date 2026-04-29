"""
Concatenate the MSD for B beads across all replicas according to whether they
interact with the lamina, and save the average MSD per compartment. One single
.h5 is generated with the average MSD per compartment for the complete
condition, split into in-contact and not-in-contact groups.
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
LAMINA_CONTACT_FOLDER = f"/scratch/mm146/Polarization/{FOLDER_NUM}_IMR-90_hg38_chr{CHROMOSOME}/7_analysis/lamina-contact/data"
NUM_REPLICAS = 32
CONDITION = "complete"
LAMINA_CONTACT_THRESHOLD = 0.20


def read_lamina_contact(filename):
    """Extract the lamina contact and B indices from file."""
    with h5py.File(filename, "r") as saved_file:
        lamina_contact = saved_file["lamina-contact"][()]
        B_indices = saved_file["B-indices"][()]

    return lamina_contact, B_indices


chr_sequence = read_sequence(
    F"../../1_inputs/chr{CHROMOSOME}_subcompartments.txt"
)

compartment_indices = {}
for subcompartment in ["B1", "B2", "B3"]:
    compartment_indices[subcompartment] = np.where(
        chr_sequence == subcompartment
    )[0]

compartment_indices["B"] = reduce(
    np.union1d,
    (
        compartment_indices["B1"],
        compartment_indices["B2"],
        compartment_indices["B3"],
    ),
)

subcompartment_masks = {
    subcompartment: chr_sequence == subcompartment
    for subcompartment in ["B1", "B2", "B3"]
}
subcompartment_masks["B"] = np.isin(chr_sequence, ["B1", "B2", "B3"])

with h5py.File(
    os.path.join(OUTPUT_FOLDER, "msd-average-per-compartment-lamina.h5"), "w"
) as f_out:
    print(
        f"Concatenating the MSD for condition: {CONDITION}",
        flush=True,
    )

    msd_groups = {
        "in-contact": {label: [] for label in compartment_indices},
        "not-in-contact": {label: [] for label in compartment_indices},
    }

    for replica in range(1, NUM_REPLICAS + 1):
        print("\tProcessing replica:", replica, flush=True)

        with h5py.File(
            os.path.join(
                OUTPUT_FOLDER,
                CONDITION,
                f"msd-all-beads-replica{replica}.h5",
            ),
            "r",
        ) as saved_file:
            msd_per_bead = saved_file["msd-per-bead"][()]
        num_lags = msd_per_bead.shape[0]

        lamina_contact, B_indices = read_lamina_contact(
            os.path.join(
                LAMINA_CONTACT_FOLDER,
                CONDITION,
                f"lamina-contact-replica{replica}.h5",
            )
        )

        average_lamina_contact = np.mean(lamina_contact, axis=0)
        selected_indices = {
            "in-contact": B_indices[
                average_lamina_contact > LAMINA_CONTACT_THRESHOLD
            ],
            "not-in-contact": B_indices[
                average_lamina_contact == 0
            ],
        }

        for group_label, group_indices in selected_indices.items():
            for label in compartment_indices:
                indices_to_average = group_indices[
                    subcompartment_masks[label][group_indices]
                ]
                if indices_to_average.size == 0:
                    continue

                msd_groups[group_label][label].append(
                    msd_per_bead[:, indices_to_average]
                )

    for group_label, grouped_msd in msd_groups.items():
        group = f_out.create_group(group_label)

        for label, msd_values in grouped_msd.items():
            if len(msd_values) == 0:
                group.create_dataset(
                    label,
                    data=np.full(num_lags, np.nan),
                )
                continue

            msd_compartment = np.concatenate(msd_values, axis=1).mean(axis=1)
            group.create_dataset(label, data=msd_compartment)

print("Done!", flush=True)
