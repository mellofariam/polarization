import os
import shutil

import chroma
from chroma import convert

CONDITIONS = ["complete", "lamina", "nuclear-bodies", "control"]
INPUT_PATH = "/scratch/mm146/Polarization/3_IMR-90_hg38_fields-nb/3_single-chr-simulation"
OUTPUT_PATH = "/scratch/mm146/Polarization/3_IMR-90_hg38_fields-nb/4_analysis/visualization"

REPLICA = 65

for condition in CONDITIONS:
    print(f"Processing condition: {condition}")
    output_folder = os.path.join(OUTPUT_PATH, condition)
    os.makedirs(output_folder, exist_ok=True)

    for molecule in range(3):
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

        traj_file = os.path.join(
            INPUT_PATH,
            condition,
            f"{REPLICA}/nucleus_{molecule}.cndb",
        )

        convert.cndb2xtc(
            cndb_file=traj_file,
            topfile=topfile,
            filename=os.path.join(
                output_folder, f"nucleus_{molecule}"
            ),
            frames=range(0, 10_000, 20),
        )
