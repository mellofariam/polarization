import os
import shutil

import chroma
from chroma import convert

CONDITIONS = [
    "complete",
    # "lamina",
    "nucleolus",
    # "control",
]
INPUT_PATH = "/scratch/mm146/Polarization/4_IMR-90_hg38_chr1/6_sampling"
OUTPUT_PATH = (
    "/scratch/mm146/Polarization/4_IMR-90_hg38_chr1/7_analysis/visualization"
)

REPLICA = 1

os.makedirs(OUTPUT_PATH, exist_ok=True)

for condition in CONDITIONS:
    print(f"Processing condition: {condition}")
    output_folder = os.path.join(OUTPUT_PATH, condition)
    os.makedirs(output_folder, exist_ok=True)

    for molecule in range(3):
        try:
            print(f"  Processing molecule: {molecule}")
            topfile = os.path.join(
                INPUT_PATH,
                condition,
                f"{REPLICA}/after-anneal_{molecule}.pdb",
            )
            shutil.copy(
                topfile,
                os.path.join(output_folder, f"nucleus_{molecule}.pdb"),
            )
        except:
            pass
        if molecule == 0:
            traj_file = os.path.join(
                INPUT_PATH,
                condition,
                f"{REPLICA}/nucleus_{molecule}.cndb",
            )

            convert.cndb2xtc(
                cndb_file=traj_file,
                topfile=topfile,
                filename=os.path.join(output_folder, f"nucleus_{molecule}"),
                frames=range(0, 10_000, 20),
            )
