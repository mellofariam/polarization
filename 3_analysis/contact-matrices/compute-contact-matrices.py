"""Load trajectory for chr1 and compute contact probability."""
import os
import time
import sys

from chroma import structure
from chroma import energy
import h5py
import numpy as np
from OpenMiChroM.CndbTools import cndbTools

traj_folder = sys.argv[1]
condition = sys.argv[2]
replica = int(sys.argv[3])
output_folder = sys.argv[4]

start = time.time()

traj = cndbTools()
traj.load(
    os.path.join(
        traj_folder, condition, str(replica), "nucleus_0.cndb"
    )
)
positions = traj.xyz(frames=range(10_000))
num_frames = positions.shape[0]

contact_probability = np.zeros(
    (positions.shape[1], positions.shape[1]), dtype=float
)

for frame in range(positions.shape[0]):
    contact_probability += energy._contact_switch(
        structure.compute_distances(
            positions[frame], positions[frame]
        )
    )
contact_probability /= num_frames

os.makedirs(os.path.join(output_folder, condition), exist_ok=True)

out_path = os.path.join(
    output_folder,
    condition,
    f"contact-probability-replica{replica}.h5",
)

with h5py.File(out_path, "w") as f:
    f.create_dataset(
        "contact_probability", data=contact_probability
    )

end = time.time()
print(
    f"Contact probability done in {end - start:.2f} seconds"
) 