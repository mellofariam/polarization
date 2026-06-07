"""Load trajectory for chr1 and compute contact probability."""

import datetime
import os
import sys
import time

import h5py
import numba
import numpy as np
import pandas
from chroma import energy, structure
from OpenMiChroM.CndbTools import cndbTools

condition = sys.argv[1]
replica = int(sys.argv[2])
traj_folder = sys.argv[3]
output_folder = sys.argv[4]

CONTACT_THRESHOLD = 0.90

start = time.perf_counter()

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
    "Computing contacts data for",
    flush=True,
)
print(f"\tcondition: {condition}", flush=True)
print(f"\treplica: {replica}", flush=True)
print("", flush=True)


traj = cndbTools()
traj.load(
    os.path.join(
        traj_folder, condition, str(replica), "nucleus_0.cndb"
    )
)
positions = traj.xyz(frames=range(0, 10_000, 1))

num_frames = positions.shape[0]
num_beads = positions.shape[1]

chr_sequence = np.array([x.decode("utf-8") for x in traj.ChromSeq])
chr_sequence_compartment = np.array([x[0] for x in chr_sequence])

is_A = chr_sequence_compartment == "A"
is_B = chr_sequence_compartment == "B"
is_N = chr_sequence_compartment == "N"


def remove_n_neighbors(matrix, diagonals_to_delete=[-2, -1, 0, 1, 2]):
    for k in diagonals_to_delete:
        i = np.arange(matrix.shape[0] - abs(k))
        if k >= 0:
            matrix[i, i + k] = 0
        else:
            matrix[i - k, i] = 0

    return matrix


@numba.njit(fastmath=True, parallel=True)
def get_patch_contacts(matrix, starts, ends):
    """
    Refer to https://github.com/numba/numba-scipy/issues/38
    """

    m = len(starts)
    C = np.zeros((m, m), matrix.dtype)

    for i in range(m):
        for j in range(i + 1, m):
            mx = -np.inf

            for r in range(starts[i], ends[i] + 1):
                for c in range(starts[j], ends[j] + 1):
                    if matrix[r, c] > mx:
                        mx = matrix[r, c]

            C[i, j] = mx
            C[j, i] = mx

    return C


df_patches = pandas.read_pickle(
    "/scratch/mm146/Polarization/polarization/1_inputs/chr1_patches.pkl"
)
num_patches = len(df_patches)

contact_probability = np.zeros(
    (positions.shape[1], positions.shape[1]), dtype=float
)
num_contacts = np.zeros(
    (positions.shape[0], positions.shape[1], 3), dtype=int
)
num_patch_contacts = np.zeros(
    (positions.shape[0], num_patches, 3),
    dtype=int,
)

for frame in range(positions.shape[0]):
    if frame % 100 == 0:
        print(f"Computing frame {frame}...", flush=True)

    frame_probabilities = energy._contact_switch(
        structure.compute_distances(
            positions[frame], positions[frame]
        )
    )

    # contact probability matrix
    contact_probability += frame_probabilities

    frame_contacts = (frame_probabilities >= CONTACT_THRESHOLD).astype(int)
    frame_contacts = remove_n_neighbors(frame_contacts)

    # number of contacts
    num_contacts[frame, :, 0] = frame_contacts[:, is_A].sum(axis=1)
    num_contacts[frame, :, 1] = frame_contacts[:, is_B].sum(axis=1)
    num_contacts[frame, :, 2] = frame_contacts[:, is_N].sum(axis=1)

    # number of patch contacts
    frame_patch_contacts = get_patch_contacts(
        frame_contacts,
        starts=df_patches["start"].values,
        ends=df_patches["end"].values,
    )

    num_patch_contacts[frame, :, 0] = frame_patch_contacts[
        :, df_patches["compartment"] == "A"
    ].sum(axis=1)
    num_patch_contacts[frame, :, 1] = frame_patch_contacts[
        :, df_patches["compartment"] == "B"
    ].sum(axis=1)
    num_patch_contacts[frame, :, 2] = frame_patch_contacts[
        :, df_patches["compartment"] == "NA"
    ].sum(axis=1)

contact_probability /= num_frames

print("Saving files...", flush=True)

os.makedirs(os.path.join(output_folder, condition), exist_ok=True)

out_path = os.path.join(
    output_folder,
    condition,
    f"contact-data-replica{replica}.h5",
)
with h5py.File(out_path, "w") as f:
    f.create_dataset("contact_probability", data=contact_probability)
    f.create_dataset("num_contacts_per_bead", data=num_contacts)
    f.create_dataset(
        "num_contacts_per_patch", data=num_patch_contacts
    )

end = time.perf_counter()
print(
    f"Contact probability done in {datetime.timedelta(seconds=end - start)}."
)
