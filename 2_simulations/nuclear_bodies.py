"""
Functions to position nuclear bodies (speckles and nucleoli) inside the nucleus
and away from chromatin.

These implement the nuclear bodies as 'fields', meaning that each body is
represented by a single point in space, and the interaction is calculated
as a function of the distance from that point.

Authors: Matheus Mello
Date: Oct 2025
"""

import numpy as np
from scipy.spatial import distance


def _calc_rotation_matrix(vec1, vec2):
    """
    Find the rotation matrix that aligns vec1 to vec2.
    See: https://stackoverflow.com/questions/45142959/calculate-rotation-matrix-to-align-two-vectors-in-3d-space

    Args:
        vec1 (np.array): A 3d "source" vector
        vec2 (np.array): A 3d "destination" vector
    
    Returns:
        np.array: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    """

    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (
        vec2 / np.linalg.norm(vec2)
    ).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    kmat = np.array(
        [[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]]
    )
    rotation_matrix = (
        np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s**2))
    )

    return rotation_matrix


def fix_chromatin(chromatin, nucleus_radius, thresh=2.0):
    """
    Bring chromatin close to the lamina by aligning the furthest bead
    to the z axis and translating it close to the lamina.

    Args:
        chromatin (np.array): Nx3 array of chromatin bead positions.
        nucleus_radius (float): Radius of the nucleus.
        thresh (float): Minimum distance between chromatin bead and lamina.
    Returns:
        chromatin (np.array): Nx3 array of transformed chromatin bead positions.
    """

    print("Bringing chromatin close to the lamina...", flush=True)

    distance_from_center = [
        np.linalg.norm(bead) for bead in chromatin
    ]

    max_distance = np.max(distance_from_center)
    max_bead = np.argmax(distance_from_center)

    # align in the z direction
    rotation_matrix = _calc_rotation_matrix(
        vec1=chromatin[max_bead], vec2=[0, 0, 1]
    )

    rotated_chromatin = np.array(
        [rotation_matrix.dot(bead) for bead in chromatin]
    )
    rotated_chromatin += np.array(
        [0, 0, nucleus_radius - max_distance - thresh]
    )

    print("Chromatin positions fixed!", flush=True)

    return rotated_chromatin


def position_nucleoli(
    chromatin, nucleoli_radius, nucleus_radius, thresh=2.0
):
    """
    Position nucleoli inside the nucleus away from chromatin.
    To be called after `fix_chromatin`, such that chromatin is
    aligned to the z axis.

    Args:
        chromatin (np.array): Nx3 array of chromatin bead positions.
        nucleoli_radius (float): Radius of the nucleolus.
        nucleus_radius (float): Radius of the nucleus.
        thresh (float): Minimum distance between nucleolus surface and chromatin.
    Returns:
        nucleoli (np.array): 1x3 array of nucleolus center positions.
    """

    print("Positioning nucleolus in the nucleus...", flush=True)
 
    zmin = np.min(chromatin, axis=0)[2]  # in z direction

    if nucleus_radius + zmin < thresh + 2 * nucleoli_radius + 1.0:
        raise ValueError(
            "Nucleoli cannot be positioned inside the nucleus!"
        )

    nucleoli = np.array([[0, 0, zmin - nucleoli_radius - thresh]])

    print("Nucleolus positioned!", flush=True)

    return nucleoli


def _random_point_in_nucleus(nucleus_radius, thresh=2.0):
    """
    Return a random 3D point uniformly distributed inside the nucleus.
    
    Args:
        nucleus_radius (float): Radius of the nucleus.
        thresh (float): Minimum distance from the lamina.
    Returns:
        np.array: 1x3 array of the random point coordinates.
    """

    # Generate a random direction
    vec = np.random.normal(0, 1, 3)
    vec /= np.linalg.norm(vec)

    # Generate a random radius with cubic scaling for uniform volume distribution
    r = (nucleus_radius - thresh) * np.cbrt(np.random.rand())

    return r * vec


def position_speckles(
    num,
    speckles_radius,
    nucleus_radius,
    chromatin,
    nucleoli_radius,
    nucleoli_CoM,
):
    """
    Distribute speckles inside the nucleus away from chromatin and nucleoli.
    To be called after `fix_chromatin`, such that chromatin is
    aligned to the z axis.

    Args:
        num (int): Number of speckles to distribute.
        nucleus_radius (float): Radius of the nucleus.
        chromatin (np.array): Nx3 array of chromatin bead positions.
        nucleoli_radius (float): Radius of the nucleolus.
        nucleoli_CoM (np.array): 1x3 array of nucleolus center position.
    Returns:
        positions (np.array): numx3 array of speckle center positions.
    """
    print("Distributing speckles...", flush=True)

    if len(nucleoli_CoM.shape) == 1:
        nucleoli_CoM = nucleoli_CoM.reshape(1, 3)

    positions = []

    i = 0
    while i < num:
        new_speckle = _random_point_in_nucleus(
            nucleus_radius,
            thresh=speckles_radius
            + 1.0,  # 0.5 of speckle + 0.5 of lamina
        ).reshape(1, 3)

        if (
            np.min(distance.cdist(new_speckle, chromatin), axis=None)
            > speckles_radius
            + 1.0  # 0.5 of speckle + 0.5 of chr bead
        ) and (
            np.min(
                distance.cdist(new_speckle, nucleoli_CoM), axis=None
            )
            > speckles_radius
            + nucleoli_radius
            + 1.0  # 0.5 of speckle + 0.5 of chr bead
        ):
            if i >= 1:
                if np.min(
                    distance.cdist(
                        new_speckle, np.concatenate(positions, axis=0)
                    )
                    > 2 * speckles_radius + 1.0, # 0.5 per speckle
                    axis=None,
                ):
                    positions.append(new_speckle)
                    i += 1
            else:
                positions.append(new_speckle)
                i += 1

    positions = np.concatenate(positions, axis=0)

    print("Speckles positioned away from chromatin!", flush=True)

    return positions


def add_nuclear_bodies(
    simulation, chromatin, speckles, nucleoli, mass=0
):
    """
    This function is to be used after MiChroM.initStructure. It adds
    nuclear bodies (speckles and nucleoli) to the simulation object,
    updating the chains and genome annotations accordingly.

    Args:
        simulation (MiChroM): MiChroM simulation object.
        chromatin (np.array): Nx3 array of chromatin bead positions.
        speckles (np.array): Mx3 array of speckle center positions.
        nucleoli (np.array): Kx3 array of nucleolus center positions.
        mass (float, optional): Mass of the nuclear body beads. Defaults to 0.
    Returns:
        positions (np.array): (N+M+K)x3 array of all bead positions.
        masses (list): List of masses for all beads.
    """

    # concatenating positions
    positions = np.concatenate(
        (chromatin, speckles, nucleoli), axis=0
    )

    # fixing chains
    simulation.chains.append(
        (
            simulation.chains[-1][1] + 1,
            simulation.chains[-1][1] + 1 + len(speckles) - 1,
            0,
        )
    )
    simulation.chains.append(
        (
            simulation.chains[-1][1] + 1,
            simulation.chains[-1][1] + 1 + len(nucleoli) - 1,
            0,
        )
    )

    # fixing genome annotations
    simulation.type_list_letter.extend(
        ["SP" for _ in range(len(speckles))]
    )
    simulation.type_list_letter.extend(
        ["NC" for _ in range(len(nucleoli))]
    )

    simulation.diff_types.update(["SP", "NC"])

    # setting masses
    masses = []
    for _ in range(len(chromatin)):
        masses.append(1)
    for _ in range(len(speckles)):
        masses.append(mass)
    for _ in range(len(nucleoli)):
        masses.append(mass)

    return positions, masses