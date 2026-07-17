import sys

import h5py
import numpy as np
import pandas

condition = sys.argv[1]
replica = sys.argv[2]
threshold = sys.argv[3]

NUM_REPLICAS = 96
NUM_BEADS = 4980
NUCLEUS_RADIUS = 32.5
NUCLEOLI_RADIUS = (NUCLEUS_RADIUS) / 5 ** (1 / 3)
CHROMATIN_DENSITY = 0.30

OUTPUT_FOLDER = f"/scratch/mm146/Polarization/4_IMR-90_hg38_chr1/7_analysis/contacts/data-{threshold}"

annotations = np.loadtxt(
    "../../1_inputs/chr1_compartments.txt", dtype=object
)[:, 1]


vmax = {
    "0.50": 20,
    "0.90": 15,
}
bins = np.arange(0, vmax[threshold] + 1, 1)
counts = np.zeros((3, len(bins) - 1), dtype=np.float64)

with h5py.File(
    f"{OUTPUT_FOLDER}/{condition}/contact-data-replica{replica}.h5",
    "a",
) as f:

    num_contacts_per_bead = f["num_contacts_per_bead"][
                    :, :, :
                ].sum(axis=-1)
    
    for iannot, annot in enumerate(["A", "B", "NA"]):
        mask = annotations == annot
        hist, _ = np.histogram(
            num_contacts_per_bead[:, mask], bins=bins
        )
        counts[iannot] += hist

    if "hist_total_contacts" in f:
        del f["hist_total_contacts"]

    f.create_group("hist_total_contacts")

    f["hist_total_contacts"].create_dataset("bins", data=bins) 

    for iannot, annot in enumerate(["A", "B", "NA"]):
        f["hist_total_contacts"].create_dataset(
            f"{annot}", data=counts[iannot]
        )
    
