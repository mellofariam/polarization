"""
Functions to position nuclear bodies (speckles and nucleoli) inside the nucleus
and away from chromatin.

These implement the nuclear bodies as beads, meaning that each body is
represented by a set of beads positioned at a spherical surface, and the
interaction is calculated in a similar way as the Type-to-Type interactions
for the chromatin beads.

Authors: Matheus Mello
Date: Apr 2026
"""

import numpy as np
import random
from scipy.spatial import distance, cKDTree
from OpenMiChroM.ChromDynamics import MiChroM


def mean_nearest_neighbor_distance(points):
    tree = cKDTree(points)
    distances, _ = tree.query(points, k=2)
    return np.mean(distances[:, 1])


def _estimate_shell_point_count(radius, target_spacing):
    """
    Estimate how many surface points are needed to achieve a target spacing.

    Uses the area per site of a hexagonal lattice as a good first-order
    approximation for evenly distributed points on a sphere.
    """
    sphere_area = 4.0 * np.pi * radius**2
    hex_cell_area = 0.5 * np.sqrt(3.0) * target_spacing**2
    return max(12, int(np.round(sphere_area / hex_cell_area)))


def fibonacci_sphere(n_points, radius=1.0, epsilon=0.36):
    """
    Evenly distribute points on a sphere. The poles of the spiral are
    placed along the X-axis (in the XY plane).

    epsilon = 0.50 -> Canonical lattice
    epsilon = 0.36 -> Off-lattice (Optimized for average nearest-neighbor)
    """
    i = np.arange(n_points, dtype=float)
    phi = (1.0 + np.sqrt(5.0)) / 2.0
    golden_angle = 2.0 * np.pi * (1.0 - 1.0 / phi)

    # Offset lattice equation from the Extreme Learning article
    # This dictates how close the first/last points get to the poles
    t = (i + epsilon) / (n_points - 1 + 2.0 * epsilon)
    spiral_axis = 1.0 - 2.0 * t

    r_eq = np.sqrt(1.0 - spiral_axis * spiral_axis)
    theta = golden_angle * i

    # By assigning the spiral axis to X, the poles are at (+1, 0, 0) and (-1, 0, 0)
    # The positive Z region now lies along the highly uniform "equator" of the spiral
    x = spiral_axis
    y = r_eq * np.cos(theta)
    z = r_eq * np.sin(theta)

    return radius * np.column_stack((x, y, z))


def _build_shell_with_spacing(radius, n_points):
    shell = fibonacci_sphere(n_points, radius=radius)
    mean_spacing = mean_nearest_neighbor_distance(shell)
    return shell, mean_spacing


def create_nuclear_body(radius: float, target_spacing: float = 1.5):
    """
    Build a Fibonacci shell whose mean nearest-neighbor spacing is as close as
    possible to the requested target.

    The point count is estimated analytically, then refined with a bracketed
    binary search to find the best integer point count.

    Args:
        radius (float, required):
            Radius of the nuclear body (distance from center to surface).
        target_spacing (float, required):
            Desired mean nearest-neighbor spacing between points on the shell.

    Returns:
        shell (np.array): Nx3 array of point coordinates on the shell.
        mean_spacing (float): Mean nearest-neighbor spacing of the returned shell.
    """
    if radius <= 0:
        raise ValueError("`radius` must be positive.")
    if target_spacing <= 0:
        raise ValueError("`target_spacing` must be positive.")

    estimate = _estimate_shell_point_count(radius, target_spacing)

    best_n_points = estimate
    best_shell, best_mean_spacing = _build_shell_with_spacing(radius, estimate)
    best_error = abs(best_mean_spacing - target_spacing)

    if best_mean_spacing > target_spacing:
        low = estimate
        high = max(estimate * 2, estimate + 1)
        high_shell, high_spacing = _build_shell_with_spacing(radius, high)

        while high_spacing > target_spacing:
            low = high
            high *= 2
            high_shell, high_spacing = _build_shell_with_spacing(radius, high)
    else:
        high = estimate
        high_shell, high_spacing = best_shell, best_mean_spacing
        low = max(12, estimate // 2)
        low_shell, low_spacing = _build_shell_with_spacing(radius, low)

        while low > 12 and low_spacing < target_spacing:
            high = low
            high_shell, high_spacing = low_shell, low_spacing
            low = max(12, low // 2)
            low_shell, low_spacing = _build_shell_with_spacing(radius, low)

        if low == 12 and low_spacing < target_spacing:
            return low_shell, low_spacing

    candidate_cache = {
        low: (low_shell, low_spacing),
        high: (high_shell, high_spacing),
        best_n_points: (best_shell, best_mean_spacing),
    }

    for n_points, (shell, mean_spacing) in candidate_cache.items():
        error = abs(mean_spacing - target_spacing)
        if error < best_error:
            best_n_points = n_points
            best_shell = shell
            best_mean_spacing = mean_spacing
            best_error = error

    while high - low > 1:
        mid = (low + high) // 2
        mid_shell, mid_spacing = _build_shell_with_spacing(radius, mid)
        candidate_cache[mid] = (mid_shell, mid_spacing)

        error = abs(mid_spacing - target_spacing)
        if error < best_error:
            best_n_points = mid
            best_shell = mid_shell
            best_mean_spacing = mid_spacing
            best_error = error

        if mid_spacing > target_spacing:
            low = mid
        else:
            high = mid

    return best_shell, best_mean_spacing


def _calc_rotation_matrix(vec1, vec2):
    """Find the rotation matrix that aligns vec1 to vec2
    :param vec1: A 3d "source" vector
    :param vec2: A 3d "destination" vector
    :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    See: https://stackoverflow.com/questions/45142959/calculate-rotation-matrix-to-align-two-vectors-in-3d-space
    """

    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (
        vec2 / np.linalg.norm(vec2)
    ).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s**2))

    return rotation_matrix


def fix_chromatin(
    chromatin,
    nucleus_radius,
    translate=True,
    thresh=2.0,
):
    """
    Bring chromatin close to the lamina by aligning the furthest bead
    to the z axis and translating it close to the lamina.

    Args:
        chromatin (np.array): Nx3 array of chromatin bead positions.
        nucleus_radius (float): Radius of the nucleus.
        translate (bool): Whether to translate chromatin close to the lamina.
        thresh (float): Minimum distance between chromatin bead and lamina,
            to be used when `translate` is `True`.
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

    if translate:
        rotated_chromatin += np.array(
            [0, 0, nucleus_radius - max_distance - thresh]
        )

    print("Chromatin positions fixed!", flush=True)

    return rotated_chromatin


def fix_nucleoli(
    chromatin, nucleoli, nucleoli_radius, nucleus_radius, thresh=2.0
):
    """
    Position a nucleolus shell inside the nucleus away from chromatin.
    To be called after `fix_chromatin`, such that chromatin is
    aligned to the z axis.

    Unlike `position_nucleoli`, which returns a single nucleolus center,
    this function translates an existing set of beads distributed on the
    nucleolus surface.

    Args:
        chromatin (np.array): Nx3 array of chromatin bead positions.
        nucleoli (np.array): Mx3 array of nucleolus bead positions.
        nucleoli_radius (float): Radius of the nucleolus.
        nucleus_radius (float): Radius of the nucleus.
        thresh (float): Minimum distance between nucleolus surface and chromatin.
    Returns:
        nucleoli (np.array): Mx3 array of translated nucleolus bead positions.
        nucleoli_center (np.array): 3D coordinates of the nucleolus center after translation.
    """

    print("Positioning nucleolus in the nucleus...", flush=True)

    zmin = np.min(chromatin, axis=0)[2]

    if nucleus_radius + zmin < thresh + 2 * nucleoli_radius + 1.0:
        raise ValueError(
            "Nucleoli cannot be positioned inside the nucleus!"
        )

    target_center = np.array([0.0, 0.0, zmin - nucleoli_radius - thresh])
    nucleoli_center = np.mean(nucleoli, axis=0)
    translated_nucleoli = nucleoli + (target_center - nucleoli_center)

    nucleolus_distance_from_center = np.linalg.norm(
        translated_nucleoli, axis=1
    )
    if np.max(nucleolus_distance_from_center) >= nucleus_radius - 1:
        raise ValueError("Nucleoli could not be positioned inside the nucleus!")

    print("Nucleolus positioned!", flush=True)

    return translated_nucleoli, target_center


def calc_half_angle(
    density: float,
    nucleus_radius: float,
    nucleolus_radius: float,
    num_chromatin_beads: int,
    bead_radius: float = 0.5,
):
    """
    Calculates the half-opening angle of the conic collapse in radians.

    Args:
        density (float): Density of the chromatin in the conic region.
        nucleus_radius (float): Radius of the nucleus.
        nucleolus_radius (float): Radius of the nucleolus.
        num_chromatin_beads (int): Number of chromatin beads.
        bead_radius (float): Radius of the chromatin beads.
    Returns: float
        Half-opening angle of the conic collapse in radians.
    """
    # 1. Total volume of all beads
    v_total_beads = num_chromatin_beads * (4 / 3) * np.pi * (bead_radius**3)

    # 2. Volume needed to satisfy the density requirement
    v_region_required = v_total_beads / density

    # 3. Geometric factor for the spherical shell: (2/3) * pi * (R^3 - r^3)
    shell_factor = (
        (2 / 3)
        * np.pi
        * (np.power(nucleus_radius, 3) - np.power(nucleolus_radius, 3))
    )

    # 4. Solve for cos(alpha)
    # V = shell_factor * (1 - cos_alpha)
    one_minus_cos_alpha = v_region_required / shell_factor
    cos_alpha = 1 - one_minus_cos_alpha

    # Clip cos_alpha to [-1, 1] to avoid NaNs from floating point errors
    cos_alpha = np.clip(cos_alpha, -1.0, 1.0)

    # 5. Return alpha in radians
    alpha_rad = np.arccos(cos_alpha)

    return alpha_rad


def filter_points_in_cone(
    points,
    half_angle,
    border_angle=0.0,
):
    """
    Keep only the points inside a cone aligned with the positive z axis.

    The cone apex is at the origin and its effective half-opening angle is
    `half_angle + border_angle`, allowing the user to keep a small angular
    border outside the chromatin confinement region.

    Args:
        points (np.array): Nx3 array of point coordinates.
        half_angle (float): Half-opening angle of the chromatin cone in radians.
        border_angle (float): Extra angle in radians added to `half_angle`.
    Returns:
        np.array: Mx3 array containing only the points inside the enlarged cone.
    """
    if half_angle < 0:
        raise ValueError("`half_angle` must be non-negative.")
    if border_angle < 0:
        raise ValueError("`border_angle` must be non-negative.")

    print("Filtering nuclear body points inside the cone...", flush=True)

    effective_half_angle = min(half_angle + border_angle, np.pi)

    radial_distance = np.linalg.norm(points, axis=1)
    cos_theta = np.ones(len(points), dtype=float)
    nonzero = radial_distance > 0
    cos_theta[nonzero] = points[nonzero, 2] / radial_distance[nonzero]
    cos_theta = np.clip(cos_theta, -1.0, 1.0)

    inside_cone = np.arccos(cos_theta) <= effective_half_angle
    filtered_points = points[inside_cone]

    print("Nuclear body points filtered!", flush=True)

    return filtered_points



def _random_point_in_nucleus(nucleus_radius, thresh=2.5):
    r = random.random() * (nucleus_radius - thresh)
    theta = random.random() * 2 * np.pi
    phi = random.random() * np.pi

    return [
        r * np.sin(phi) * np.cos(theta),
        r * np.sin(phi) * np.sin(theta),
        r * np.cos(phi),
    ]


def add_nuclear_bodies(
    simulation: MiChroM,
    chromatin: np.ndarray,
    lamina: None | np.ndarray = None,
    speckles: None | np.ndarray = None,
    nucleoli: None | np.ndarray = None,
    mass=0,
):
    """
    This function is to be used after MiChroM.initStructure. It adds
    nuclear bodies (lamina, speckles, and nucleoli) to the simulation object,
    updating the chains and genome annotations accordingly.

    Args:
        simulation (MiChroM): MiChroM simulation object.
        chromatin (np.array): Nx3 array of chromatin bead positions.
        lamina (np.array): Lx3 array of lamina bead positions.
        speckles (np.array): Mx3 array of speckle center positions.
        nucleoli (np.array): Kx3 array of nucleolus center positions.
        mass (float, optional): Mass of the nuclear body beads. Defaults to 0.
    Returns:
        positions (np.array): (N+L+M+K)x3 array of all bead positions.
        masses (list): List of masses for all beads.
    """

    if lamina is None and speckles is None and nucleoli is None:
        raise ValueError(
            "At least one of the nuclear bodies must be provided!"
        )

    positions = [chromatin]
    masses = [1 for _ in range(len(chromatin))]

    print("Adding nuclear bodies to the system...", flush=True)

    if lamina is not None:
        print(
            f"\t- lamina: chain {len(simulation.chains)}",
            flush=True,
        )
        positions.append(lamina)
        simulation.chains.append(
            (
                simulation.chains[-1][1] + 1,
                simulation.chains[-1][1] + 1 + len(lamina) - 1,
                0,
            )
        )
        simulation.type_list_letter.extend(
            ["LM" for _ in range(len(lamina))]
        )
        masses.extend([mass for _ in range(len(lamina))])
        simulation.diff_types.update(["LM"])

    if speckles is not None:
        print(
            f"\t- speckles: chain {len(simulation.chains)}",
            flush=True,
        )
        positions.append(speckles)
        simulation.chains.append(
            (
                simulation.chains[-1][1] + 1,
                simulation.chains[-1][1] + 1 + len(speckles) - 1,
                0,
            )
        )
        simulation.type_list_letter.extend(
            ["SP" for _ in range(len(speckles))]
        )
        masses.extend([mass for _ in range(len(speckles))])
        simulation.diff_types.update(["SP"])

    if nucleoli is not None:
        print(
            f"\t- nucleoli: chain {len(simulation.chains)}",
            flush=True,
        )
        positions.append(nucleoli)
        simulation.chains.append(
            (
                simulation.chains[-1][1] + 1,
                simulation.chains[-1][1] + 1 + len(nucleoli) - 1,
                0,
            )
        )
        simulation.type_list_letter.extend(
            ["NC" for _ in range(len(nucleoli))]
        )
        masses.extend([mass for _ in range(len(nucleoli))])
        simulation.diff_types.update(["NC"])

    positions = np.concatenate(positions, axis=0)

    return positions, masses
