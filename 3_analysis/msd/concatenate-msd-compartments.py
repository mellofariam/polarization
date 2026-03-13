import os
import sys

import h5py
import freud
import numpy as np
from OpenMiChroM.CndbTools import cndbTools
from functools import reduce

from chroma import structure


def read_sequence(filename):
    """Extract subcompartment information from file."""
    sequence = np.loadtxt(filename, dtype=str)

    return sequence[:, 1].astype(str)


conditions = [
    # "control",
    "complete",
    # "lamina",
    # "nucleolus",
    # "isolation",
    # "nuclear-bodies",
]

chr_sequence = {
    10: read_sequence(
        "/scratch/mm146/Polarization/3_IMR-90_hg38_fields-nb/1_inputs/chr10_IMR90_hg38.subcmpt"
    ),
    11: read_sequence(
        "/scratch/mm146/Polarization/3_IMR-90_hg38_fields-nb/1_inputs/chr11_IMR90_hg38.subcmpt"
    ),
}

compartment_indices = {}
for chromosome in [10, 11]:
    compartment_indices[chromosome] = {}
    for subcompartment in ["A1", "A2", "B1", "B2", "B3", "NA"]:
        compartment_indices[chromosome][subcompartment] = np.where(
            chr_sequence[chromosome] == subcompartment
        )[0]

    compartment_indices[chromosome]["A"] = reduce(
        np.union1d,
        (
            compartment_indices[chromosome]["A1"],
            compartment_indices[chromosome]["A2"],
        ),
    )
    compartment_indices[chromosome]["B"] = reduce(
        np.union1d,
        (
            compartment_indices[chromosome]["B1"],
            compartment_indices[chromosome]["B2"],
            compartment_indices[chromosome]["B3"],
        ),
    )

with h5py.File("data/msd-average-per-compartment.h5", "w") as f_out:
    for num_chr in [
        "one",
        # "two",
    ]:
        f_out.create_group(num_chr)
        if num_chr == "one":
            FOLDER_PATH = "/scratch/mm146/Polarization/3_IMR-90_hg38_fields-nb/3_single-chr-simulation"
            chromosomes = [10]
        elif num_chr == "two":
            FOLDER_PATH = "/scratch/mm146/Polarization/3_IMR-90_hg38_fields-nb/3_single-chr-simulation"
            chromosomes = [10, 11]
        else:
            raise ValueError("Invalid number of chromosomes!")

        for condition in conditions:
            f_out[num_chr].create_group(condition)
            for chromosome in chromosomes:
                f_out[num_chr][condition].create_group(
                    f"chr{chromosome}"
                )
                print(
                    f"Concatenating the MSD for chromosome: {chromosome}",
                    flush=True,
                )
                print(
                    f"\tnumber of chromosomes: {num_chr}", flush=True
                )
                print(f"\tcondititon: {condition}", flush=True)
                print("", flush=True)

                msd_all_replicas = []
                for replica in range(1, 31):
                    print(
                        "\tProcessing replica:", replica, flush=True
                    )
                    with h5py.File(
                        f"data/{num_chr}/{condition}/msd-all-beads-chr{chromosome}-replica-{replica}.h5",
                        "r",
                    ) as f:
                        msd_all_replicas.append(f["msd-per-bead"][()])

                msd_all_replicas = np.array(msd_all_replicas)

                for label in compartment_indices[chromosome].keys():
                    indices = compartment_indices[chromosome][label]
                    msd_compartment = (
                        msd_all_replicas[:, :, indices]
                        .transpose(1, 0, 2)
                        .reshape(10000, -1)
                        .mean(axis=1)
                    )

                    f_out[num_chr][condition][
                        f"chr{chromosome}"
                    ].create_dataset(
                        label,
                        data=msd_compartment,
                    )


print("Done!", flush=True)
