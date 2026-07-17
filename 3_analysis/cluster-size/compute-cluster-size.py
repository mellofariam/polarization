"""Load CNDB trajectory and compute cluster sizes."""

import datetime
import os
import sys
import time

import h5py
import numpy as np
from OpenMiChroM.CndbTools import cndbTools
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

condition = sys.argv[1]
replica = int(sys.argv[2])
traj_folder = sys.argv[3]
output_folder = sys.argv[4]

CONTACT_THRESHOLD = float(sys.argv[5])
MU = 3.22
RC = 1.78

rc = RC - 1 / MU * np.arctanh(2 * CONTACT_THRESHOLD - 1)

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
    "Computing cluster sizes for",
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

chr_sequence = np.array([x.decode("utf-8") for x in traj.ChromSeq])
chr_sequence_compartment = np.array([x[0] for x in chr_sequence])

is_A = chr_sequence_compartment == "A"
is_B = chr_sequence_compartment == "B"
is_N = chr_sequence_compartment == "N"

positions = positions[:, is_B, :]
num_beads = positions.shape[1]

backbone_skip = 2
fmax = np.zeros(num_frames)
for frame in range(num_frames):
    if frame % 100 == 0:
        print(f"Computing frame {frame}...", flush=True)

    tree = cKDTree(positions[frame])

    pairs = tree.query_pairs(rc)

    pairs = [(i, j) for i, j in pairs if abs(i - j) > backbone_skip]

    graph = coo_matrix(
        (
            np.ones(len(pairs)),
            ([p[0] for p in pairs], [p[1] for p in pairs]),
        ),
        shape=(num_beads, num_beads),
    )

    graph = graph + graph.T

    n_clusters, labels = connected_components(graph, directed=False)

    cluster_sizes = np.bincount(labels)

    fmax[frame] = cluster_sizes.max() / num_beads

print("Saving files...", flush=True)

os.makedirs(os.path.join(output_folder, condition), exist_ok=True)

out_path = os.path.join(
    output_folder,
    condition,
    f"cluster-size-replica{replica}.h5",
)
with h5py.File(out_path, "w") as f:
    f.create_dataset("cluster_size", data=fmax)

end = time.perf_counter()
print(
    f"Cluster size done in {datetime.timedelta(seconds=end - start)}."
)
