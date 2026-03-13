import collections
import itertools
import os
import pickle
import sys
import time
from itertools import combinations_with_replacement

import freud
import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
import pandas
import seaborn as sns
from matplotlib import gridspec
from OpenMiChroM.CndbTools import cndbTools as ctools
from mpl_toolkits.axes_grid1 import make_axes_locatable

import chroma

plt.style.use("chroma.paper")

import numpy as np

condition = sys.argv[1]


def analyze_points(pos, sphere_center):
    """
    pos: array of shape (t, N, 3)
    sphere_center: array-like shape (3,)
    Returns:
        com           : (t, 3)
        radial_proj   : (t, N)     projections on center→COM direction
        lateral_proj  : (t, N)     magnitudes of perpendicular projection
    """
    pos = np.asarray(pos)
    C = np.asarray(sphere_center)

    # 1. Center of mass for each frame
    com = pos.mean(axis=1)  # (t, 3)

    # 2. Direction vector center → COM per frame
    dir_vec = com - C  # (t, 3)
    dir_norm = np.linalg.norm(dir_vec, axis=1, keepdims=True)
    unit_dir = dir_vec / dir_norm  # (t, 3)

    # 3. Recenter coordinates relative to COM
    rel = pos - com[:, None, :]  # (t, N, 3)

    # 4. Projection of each point on radial direction
    radial_proj = np.sum(rel * unit_dir[:, None, :], axis=2)  # (t, N)

    # 5. Projection perpendicular to radial direction
    # component along direction
    parallel_component = (
        radial_proj[..., None] * unit_dir[:, None, :]
    )  # (t, N, 3)

    # subtract to get perpendicular component and take magnitude
    perp = rel - parallel_component
    lateral_proj = np.linalg.norm(perp, axis=2)  # (t, N)

    return com, radial_proj, lateral_proj


radial_projection = []
lateral_projection = []
angle_with_zaxis = []
angle_with_com = []

for replica in range(31, 61):
    print(f"Processing replica {replica}...", flush=True)
    traj = ctools()
    traj.load(
        f"/scratch/mm146/Polarization/3_IMR-90_hg38_fields-nb/3_single-chr-simulation/{condition}/{replica}/nucleus_0.cndb"
    )
    positions = traj.xyz(frames=range(0, 10_000, 100))

    CoM, radial_proj, lateral_proj = analyze_points(
        positions, sphere_center=[0, 0, 0]
    )

    angle_with_zaxis.extend(
        np.max(
            np.arccos(
                positions[:, :, 2] / np.linalg.norm(positions, axis=2)
            ),
            axis=1,
        )
    )
    angle_with_com.extend(
        np.max(
            np.arccos(
                np.clip(
                    np.sum(
                        np.mean(positions, axis=1, keepdims=True)
                        * positions,
                        axis=2,
                    )
                    / (
                        np.linalg.norm(
                            np.mean(positions, axis=1, keepdims=True),
                            axis=2,
                        )
                        * np.linalg.norm(positions, axis=2)
                    ),
                    -1.0,
                    1.0,
                )
            ),
            axis=1,
        )
    )

    radial_projection.extend(np.max(radial_proj, axis=1))
    radial_projection.extend(np.abs(np.min(radial_proj, axis=1)))
    lateral_projection.extend(np.max(lateral_proj, axis=1))

radial_projection = np.array(radial_projection)
lateral_projection = np.array(lateral_projection)
angle_with_zaxis = np.array(angle_with_zaxis)
angle_with_com = np.array(angle_with_com)

np.savez(
    f"projections-{condition}-center.npz",
    radial_projection=radial_projection,
    lateral_projection=lateral_projection,
    angle_with_zaxis=angle_with_zaxis,
    angle_with_com=angle_with_com,
)
