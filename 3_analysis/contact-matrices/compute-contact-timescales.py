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

df_patches = pandas.read_pickle("../../1_inputs/chr1_patches.pkl")


def contact_durations(x):
    """
    Parameters
    ----------
    x : array-like, shape (Ntime,) or (Nreplicas, Ntime)
        Binary trajectories: 0 = no contact, 1 = contact.

    Returns
    -------
    t_on : np.ndarray
        All consecutive 1-run lengths across all replicas.
    t_off : np.ndarray
        All consecutive 0-run lengths across all replicas.
    """
    x = np.asarray(x, dtype=np.int8)

    if x.ndim == 1:
        x = x[None, :]

    t_on_all = []
    t_off_all = []

    for row in x:
        if row.size == 0:
            continue

        changes = np.flatnonzero(np.diff(row)) + 1
        boundaries = np.r_[0, changes, row.size]

        run_lengths = np.diff(boundaries)
        run_values = row[boundaries[:-1]]

        t_on_all.append(run_lengths[run_values == 1])
        t_off_all.append(run_lengths[run_values == 0])

    t_on = (
        np.concatenate(t_on_all)
        if t_on_all
        else np.array([], dtype=int)
    )
    t_off = (
        np.concatenate(t_off_all)
        if t_off_all
        else np.array([], dtype=int)
    )

    return t_on, t_off


def data_moments(data):
    """
    Compute the first four moments of a 1D dataset.

    Parameters
    ----------
    data : array_like
        One-dimensional array containing the data.

    Returns
    -------
    tuple
        (mean, variance, skewness, kurtosis)

        - mean: arithmetic mean
        - variance: population variance
        - skewness: standardized third central moment
        - kurtosis: Pearson kurtosis (not excess kurtosis)

    Notes
    -----
    The variance is computed with denominator N (ddof=0), matching the
    histogram_moment() function. Kurtosis is the Pearson kurtosis, so a
    normal distribution has kurtosis equal to 3.
    """

    data = np.asarray(data, dtype=np.float64)

    if data.ndim != 1:
        raise ValueError("data must be a 1D array")
    if data.size == 0:
        return (np.nan, np.nan, np.nan, np.nan)

    mean = np.mean(data)

    centered = data - mean
    variance = np.mean(centered**2)

    if variance == 0:
        return (mean, variance, 0.0, 0.0)

    skewness = np.mean(centered**3) / variance**1.5
    kurtosis = np.mean(centered**4) / variance**2

    return mean, variance, skewness, kurtosis


with h5py.File(
    f"{OUTPUT_FOLDER}/{condition}/contact-data-replica{replica}.h5",
    "a",
) as f:
    in_contact_per_bead = (f["num_contacts_per_bead"][()] > 0).astype(
        int
    )

    if "t_on" in f:
        del f["t_on"]
    if "t_off" in f:
        del f["t_off"]

    f.create_group("t_on")
    f.create_group("t_off")

    f["t_on"].create_dataset("bins_A", data=np.linspace(0, 5000, 51)) # 1000
    f["t_off"].create_dataset("bins_A", data=np.linspace(0, 200, 51)) # 300
    f["t_on"].create_dataset("bins_B", data=np.linspace(0, 10000, 51)) # 3000
    f["t_off"].create_dataset("bins_B", data=np.linspace(0, 600, 51)) # 600

    for idx, cmpt in enumerate(["A", "B"]):
        cmpt_contacts = in_contact_per_bead[
            :, annotations == cmpt, idx
        ]
        cmpt_contacts = cmpt_contacts.transpose()
        t_on, t_off = contact_durations(cmpt_contacts)

        f["t_on"].create_dataset(cmpt, data=t_on)
        f["t_on"].create_dataset(
            f"hist_{cmpt}",
            data=np.histogram(t_on, bins=f["t_on"][f"bins_{cmpt}"][()])[0],
        )
        f["t_on"].create_dataset(
            f"moments_{cmpt}", data=data_moments(t_on)
        )

        f["t_off"].create_dataset(cmpt, data=t_off)
        f["t_off"].create_dataset(
            f"hist_{cmpt}",
            data=np.histogram(t_off, bins=f["t_off"][f"bins_{cmpt}"][()])[0],
        )
        f["t_off"].create_dataset(
            f"moments_{cmpt}", data=data_moments(t_off)
        )

    del in_contact_per_bead, cmpt_contacts, t_on, t_off

    in_contact_per_patch = (
        f["num_contacts_per_patch"][()] > 0
    ).astype(int)

    f["t_on"].create_dataset(
        "bins_patch", data=np.linspace(0, 10000, 101)
    )
    f["t_off"].create_dataset(
        "bins_patch", data=np.linspace(0, 10000, 51)
    )
    for idx, cmpt in enumerate(["A", "B"]):
        cmpt_contacts = in_contact_per_patch[
            :, df_patches["compartment"] == cmpt, idx
        ]
        cmpt_contacts = cmpt_contacts.transpose()
        t_on, t_off = contact_durations(cmpt_contacts)

        f["t_on"].create_dataset(f"patch_{cmpt}", data=t_on)
        f["t_on"].create_dataset(
            f"hist_patch_{cmpt}",
            data=np.histogram(t_on, bins=f["t_on"]["bins_patch"][()])[
                0
            ],
        )

        f["t_off"].create_dataset(f"patch_{cmpt}", data=t_off)
        f["t_off"].create_dataset(
            f"hist_patch_{cmpt}",
            data=np.histogram(
                t_off, bins=f["t_off"]["bins_patch"][()]
            )[0],
        )
