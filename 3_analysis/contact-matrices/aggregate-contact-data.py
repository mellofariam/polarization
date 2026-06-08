import os
import sys
import shutil

import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import yaml

from chroma import structure, plotting
from OpenMiChroM.CndbTools import cndbTools
import seaborn as sns
import pandas

from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

plt.style.use("chroma.tex")

conditions = ["complete", "nucleolus"]

NUM_REPLICAS = 96
NUM_BEADS = 4980
NUCLEUS_RADIUS = 32.5
NUCLEOLI_RADIUS = (NUCLEUS_RADIUS) / 5 ** (1 / 3)
CHROMATIN_DENSITY = 0.30

OUTPUT_FOLDER = "/scratch/mm146/Polarization/4_IMR-90_hg38_chr1/7_analysis/contacts/data"

contacts = {}

for condition in conditions:
    contacts[condition] = {}
    for replica in range(1, NUM_REPLICAS + 1):
        print(f"Processing {condition} replica {replica}")
        with h5py.File(
            f"{OUTPUT_FOLDER}/{condition}/contact-data-replica{replica}.h5",
            "r",
        ) as f:
            for key in f.keys():
                if key not in contacts[condition]:
                    contacts[condition][key] = []
                contacts[condition][key].append(f[key][()])

    contacts[condition]["contact_probability"] = np.mean(
        contacts[condition]["contact_probability"], axis=0
    )
    contacts[condition]["num_contacts_per_bead"] = np.array(
        contacts[condition]["num_contacts_per_bead"]
    )
    contacts[condition]["num_contacts_per_patch"] = np.array(
        contacts[condition]["num_contacts_per_patch"]
    )

with h5py.File(
    os.path.join(OUTPUT_FOLDER, "contact-data-aggregated.h5"), "w"
) as f:
    for condition in contacts.keys():
        grp = f.create_group(condition)
        for key in contacts[condition].keys():
            grp.create_dataset(key, data=contacts[condition][key])
