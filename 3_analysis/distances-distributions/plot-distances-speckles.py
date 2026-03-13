# %%
import sys

import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas
from mpl_toolkits.axes_grid1 import make_axes_locatable

plt.style.use("~/chroma/chroma/paper.mplstyle")

import chroma

plt.style.use("chroma.paper")

# %%
conditions = [
    # "control",
    "complete",
    # "lamina",
    # "nucleolus",
    # "isolation",
    # "nuclear-bodies",
]

# %%
subcompartments = {}

for chr in [10, 11]:
    subcompartments[chr] = pandas.read_csv(
        f"/scratch/mm146/Polarization/2_IMR-90_hg38/1_inputs/chr{chr}_IMR90_hg38.subcmpt",
        sep=" ",
        header=None,
        dtype=str,
        keep_default_na=False,
        names=["index", "annotation"],
    ).annotation.values

# %%
set_subcompartments = ["A1", "A2", "B1", "B2", "B3", "NA"]


def plot_subcompartment_annotations(ax, subcompartments):
    colors = {
        "A1": "red",
        "A2": "orange",
        "B1": "blue",
        "B2": "cyan",
        "B3": "green",
        "B4": "lightgreen",
        "NA": "grey",
    }

    color_sequence = list(map(colors.get, subcompartments))

    for i in range(len(subcompartments)):
        rect1 = mpl.patches.Rectangle(
            (i, 0),
            1,
            1,
            facecolor=color_sequence[i],
            edgecolor="None",
        )
        ax.add_patch(rect1)
    ax.set_xlim(0, len(subcompartments))
    ax.tick_params(
        bottom=False,
        left=False,
        top=False,
        labelbottom=False,
        labeltop=False,
        labelleft=False,
        labelright=False,
    )


# %% [markdown]
# # Distance to the CoM


# %%
def plot_average_distance_speckles(
    average_distance_speckles,
    subcompartments,
    num_chr,
    condition,
    chr,
):
    fig, ax = plt.subplots(1, 1, figsize=(3, 2))

    chr_size = len(subcompartments)

    ax.plot(
        np.arange(chr_size) * 50e3 / 1e6,
        average_distance_speckles,
        color="black",
    )

    divider = make_axes_locatable(ax)

    ax_top = divider.append_axes("top", size=0.05, pad=0.05)

    plot_subcompartment_annotations(ax_top, subcompartments)

    if num_chr == "one":
        title_head = "single chr. simulation"
    elif num_chr == "two":
        title_head = "two chr. simulation"

    title = f"{title_head} | {condition} | chr{chr}"

    ax_top.set_title(title)

    ax.set_ylabel(r"Distance to speckles ($\sigma$)")
    ax.set_xlabel("Genomic position (Mb)")

    ax.set_xlim((0, chr_size * 50e3 / 1e6))

    fig.tight_layout()
    fig.savefig(
        f"figures/distance-speckles/fig-average-distance-speckles-{num_chr}-{condition}-chr{chr}.svg",
    )

    plt.close(fig=fig)


def plot_distribution_distance_speckles(
    distance_speckles, num_chr, condition, chr
):

    fig, ax = plt.subplots(1, 1, figsize=(3, 2), layout="constrained")

    ax.hist(
        distance_speckles["A1"] + distance_speckles["A2"],
        bins=50,
        color="red",
        density=True,
        zorder=4,
        label="A",
        histtype="step",
    )
    ax.hist(
        distance_speckles["B1"]
        + distance_speckles["B2"]
        + distance_speckles["B3"],
        bins=50,
        color="blue",
        density=True,
        zorder=3,
        label="B",
        histtype="step",
    )
    ax.hist(
        distance_speckles["NA"],
        bins=50,
        color="grey",
        density=True,
        zorder=1,
        label="NA",
        histtype="step",
    )

    ax.legend()

    if num_chr == "one":
        title_head = "single chr. simulation"
    elif num_chr == "two":
        title_head = "two chr. simulation"

    title = f"{title_head} | {condition} | chr{chr}"

    ax.set_title(title)

    ax.set_xlabel(r"Distance to speckles ($\sigma$)")
    ax.set_ylabel("Probability Density")

    fig.savefig(
        f"figures/distance-speckles/fig-distribution-distance-speckles-{num_chr}-{condition}-chr{chr}.svg",
    )

    plt.close(fig=fig)


# %%
for num_chr in ["one"]:
    for condition in conditions:
        if num_chr == "one":
            chromosomes = [10]
        elif num_chr == "two":
            chromosomes = [10, 11]
        for chr in chromosomes:
            print
            average_distance_speckles = []
            distance_speckles = {
                annot: [] for annot in set_subcompartments
            }
            for replica in range(1, 31):
                with h5py.File(
                    f"data/{num_chr}/{condition}/distance-speckles-chr{chr}-replica{replica}.h5"
                ) as saved_file:
                    average_distance_speckles.append(
                        saved_file["average-distance-speckles"][()]
                    )
                    distance_speckles_all_beads = saved_file[
                        "distance-speckles"
                    ][()]

                    for annot in set_subcompartments:
                        distance_speckles[annot].extend(
                            distance_speckles_all_beads[
                                :, subcompartments[chr] == annot
                            ].reshape(-1)
                        )

            average_distance_speckles = np.mean(
                average_distance_speckles, axis=0
            )

            plot_average_distance_speckles(
                average_distance_speckles,
                subcompartments[chr],
                num_chr,
                condition,
                chr,
            )

            plot_distribution_distance_speckles(
                distance_speckles,
                num_chr,
                condition,
                chr,
            )
