"""
This script samples the chromosome in the nucleus with the updated
explicit-bead representation of the lamina and nuclear bodies.

The spherical lamina confinement and the nucleolus excluded volume remain
field-based, while the attractive interactions are mediated by explicit
lamina and nucleolus beads when enabled by the selected condition.
"""

import os
import sys

import numpy as np
from OpenMiChroM.ChromDynamics import MiChroM
from scipy.spatial import distance

import nuclear_bodies as nb

output_base = sys.argv[1]
condition = sys.argv[2]
replicaID = int(sys.argv[3])
input_base = sys.argv[4]
lamina_divide_by = sys.argv[5]
platform = sys.argv[6] if len(sys.argv) > 6 else "CUDA"

OUTPUT_FOLDER = os.path.join(
    output_base, lamina_divide_by, condition, str(replicaID)
)
os.makedirs(OUTPUT_FOLDER, exist_ok=True)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
TYPES_TABLE = os.path.join(
    SCRIPT_DIR,
    f"../1_inputs/ff_nucleus-compartment-level_divided-by-{lamina_divide_by}.csv",
)

## Defining parameters
NUCLEUS_RADIUS = 32.5
NUCLEOLI_RADIUS = NUCLEUS_RADIUS / 5 ** (1 / 3)
NUCLEAR_BODY_SPACING = 1.0
CHROMATIN_DENSITY = 0.30

###

nucleus = MiChroM(name="nucleus", temperature=1.0, timeStep=0.01)
nucleus.setup(platform=platform)

# Output folder:
nucleus.saveFolder(OUTPUT_FOLDER)

# Loading chromatin chain saved by the assembly step:
chromosomes = nucleus.initStructure(
    mode="pdb",
    CoordFiles=[
        os.path.join(
            input_base,
            condition,
            str(replicaID),
            "initial-state_0.pdb",
        )
    ],
)
num_beads = len(chromosomes)

chromosomes = nb.fix_chromatin(
    chromosomes, nucleus_radius=NUCLEUS_RADIUS, translate=False
)
print("zmax:", np.max(chromosomes[:, 2]))
print("zmin:", np.min(chromosomes[:, 2]))

half_angle_conic_confinement = nb.calc_half_angle(
    density=CHROMATIN_DENSITY,
    nucleus_radius=NUCLEUS_RADIUS,
    nucleolus_radius=NUCLEOLI_RADIUS,
    num_chromatin_beads=num_beads,
)

include_lamina = condition in ["complete", "lamina"]
include_nucleoli = condition in [
    "complete",
    "nuclear-bodies",
    "nucleolus",
]

lamina = None
if include_lamina:
    lamina, lamina_spacing = nb.create_nuclear_body(
        radius=NUCLEUS_RADIUS,
        target_spacing=NUCLEAR_BODY_SPACING,
    )

    lamina = nb.filter_points_in_cone(
        lamina,
        half_angle_conic_confinement,
        0.25 * half_angle_conic_confinement,
    )

    print(
        "Lamina created with",
        len(lamina),
        "points and mean spacing",
        lamina_spacing,
        flush=True,
    )

nucleoli_center = np.zeros(3, dtype=float)
nucleoli = None
if include_nucleoli:
    nucleoli, nucleoli_spacing = nb.create_nuclear_body(
        radius=NUCLEOLI_RADIUS,
        target_spacing=NUCLEAR_BODY_SPACING,
    )
    nucleoli -= np.mean(nucleoli, axis=0)

    nucleoli = nb.filter_points_in_cone(
        nucleoli,
        half_angle_conic_confinement,
        0.25 * half_angle_conic_confinement,
    )

    print(
        "Nucleolus created with",
        len(nucleoli),
        "points and mean spacing",
        nucleoli_spacing,
        flush=True,
    )

chromatin_chains = list(range(len(nucleus.chains)))
positions, masses = nb.add_nuclear_bodies(
    simulation=nucleus,
    chromatin=chromosomes,
    lamina=lamina,
    nucleoli=nucleoli,
)

if nucleoli is not None:
    print("Min distance between nucleoli and chromatin:", flush=True)
    print(
        np.min(distance.cdist(nucleoli, chromosomes)),
        flush=True,
    )

# Loading all particles in the simulation context
nucleus.loadStructure(positions, center=False, massList=masses)

## Including the forces in the system

# Homopolymer Potentials
nucleus.addFENEBonds(chainIndices=chromatin_chains)
nucleus.addAngles(chainIndices=chromatin_chains)
nucleus.addSelfAvoidance(chainIndices=chromatin_chains)

# Chromosome and explicit nuclear-body bead potentials
nucleus.addCustomTypes(
    mu=3.22,
    rc=1.78,
    TypesTable=TYPES_TABLE,
)

nucleus.addMultiChainIC(
    mu=3.22, rc=1.78, dinit=3, dend=500, chainIndex=0
)

# Add confinement
nucleus.addSphericalConfinementTanh(radius=NUCLEUS_RADIUS, Ecut=20)


nucleus.addConicConfinement(
    k=5e-3, halfOpeningAngle=half_angle_conic_confinement
)

nucleus.addNuclearBodiesExcludedVolumeFromPoint(
    nuclearBodyRadius=NUCLEOLI_RADIUS,
    nuclearBodyCenter=nucleoli_center,
    chromatinChainIndices=chromatin_chains,
    nuclearBodyName="nucleoli",
    Ecut=20,
)

# Set up the simulation
nucleus.createSimulation()
nucleus.saveStructure(fileName="setup", mode="pdb")

nucleus.run(
    nsteps=500 * 10**3,
    report=True,
    interval=10**4,
    totalSteps=13_000 * 10**3,
)

## Running the annealing to and back from T = 2.0

for T in np.linspace(1.1, 2.0, num=100):
    nucleus.integrator.setTemperature(T / 0.008314)
    nucleus.run(
        nsteps=10 * 10**3,
        report=True,
        interval=10**4,
        totalSteps=13_000 * 10**3,
    )

## Running the annealing to return to T = 1.0

for T in np.linspace(1.99, 1.00, num=100):
    nucleus.integrator.setTemperature(T / 0.008314)
    nucleus.run(
        nsteps=10 * 10**3,
        report=True,
        interval=10**4,
        totalSteps=13_000 * 10**3,
    )

nucleus.run(
    nsteps=500 * 10**3,
    report=True,
    interval=10**4,
    totalSteps=13_000 * 10**3,
)

nucleus.saveStructure(fileName="after-anneal", mode="pdb")

### Starting sampling part of the simulation

# Initiating .cndb file
print("Initiating sampling...", flush=True)
nucleus.createReporters(
    statistics=True,
    energyComponents=True,
    traj=True,
    interval=10**3,
)

nucleus.run(
    nsteps=10_000 * 10**3,
    report=True,
    interval=10**4,
    totalSteps=13_000 * 10**3,
)

nucleus.saveStructure(fileName="last-frame", mode="pdb")

print("Forces in the end:", flush=True)
nucleus.printForces()

print("Files closed! All set!", flush=True)
