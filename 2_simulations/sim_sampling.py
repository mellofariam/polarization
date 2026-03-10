"""
This script simulated the chromosome in the nucleus, in the presence of
nuclear bodies and the lamina.

Here, we include only interactions with the nucleolus.
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
platform = sys.argv[5] if len(sys.argv) > 5 else "CUDA"

OUTPUT_FOLDER = os.path.join(output_base, condition, str(replicaID))
os.makedirs(OUTPUT_FOLDER, exist_ok=True)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
TYPES_TABLE = os.path.join(
    SCRIPT_DIR, "../1_inputs/ff_compartments-and-nb.csv"
)

## Defining parameters
NUCLEUS_RADIUS = 32.5
NUCLEOLI_RADIUS = (NUCLEUS_RADIUS) / 5 ** (1 / 3)
SPECKLES_RADIUS = 1.5
SPECKLES_NUMBER = 40
CHROMATIN_DENSITY = 0.30

###

nucleus = MiChroM(name="nucleus", temperature=1.0, timeStep=0.01)
nucleus.setup(platform=platform)

# Output folder:
nucleus.saveFolder(OUTPUT_FOLDER)

# Loading individual chromosomes:
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

# nucleolus positioned int the center of the nucleus
nucleolus_xyz = np.array([[0.0, 0.0, 0.0]])

chromatin_chains = list(range(len(nucleus.chains)))
positions, masses = nb.add_nuclear_bodies(
    simulation=nucleus,
    chromatin=chromosomes,
    nucleoli=nucleolus_xyz,
)
nucleoli_chains = [chromatin_chains[-1] + 1]

print("Chromatin chains:", chromatin_chains)
print("Nucleoli chains:", nucleoli_chains)

print("Min distance between nucleoli and chromatin:", flush=True)
print(
    np.min(distance.cdist(nucleolus_xyz, chromosomes)),
    flush=True,
)

# Loading chromosomes in the simulation context
nucleus.loadStructure(positions, center=False, massList=masses)

## Including the forces in the system

# Homopolymer Potentials
nucleus.addFENEBonds(chainIndices=chromatin_chains)
nucleus.addAngles(chainIndices=chromatin_chains)

nucleus.addSelfAvoidance(chainIndices=chromatin_chains)

# Chromosome Potentials
nucleus.addCustomTypes(
    mu=3.22,
    rc=1.78,
    TypesTable=TYPES_TABLE,
    chainIndices=chromatin_chains,
)

nucleus.addMultiChainIC(
    mu=3.22, rc=1.78, dinit=3, dend=500, chainIndex=0
)

# Nuclear Bodies Potentials

## Nucleolus
if condition in ["complete", "nuclear-bodies", "nucleolus"]:
    nucleus.addNuclearBodiesInteraction(
        nuclearBodyRadius=NUCLEOLI_RADIUS,
        chromatinChainIndices=chromatin_chains,
        nuclearBodyChainIndices=nucleoli_chains,
        forceName="Nucleoli",
        forceNumber=1,
        mu=3.22,
        rc=1.78,
        TypesTable=TYPES_TABLE,
    )
nucleus.addNuclearBodiesExcludedVolume(
    nuclearBodyRadius=NUCLEOLI_RADIUS,
    chromatinChainIndices=chromatin_chains,
    nuclearBodyChainIndices=nucleoli_chains,
    forceName="Nucleoli_EV",
    forceNumber=1,
)

# Add confinement
nucleus.addSphericalConfinementTanh(radius=NUCLEUS_RADIUS)

half_angle_conic_confinement = nb.calc_half_angle(
    density=CHROMATIN_DENSITY,
    nucleus_radius=NUCLEUS_RADIUS,
    nucleolus_radius=NUCLEOLI_RADIUS,
    num_chromatin_beads=num_beads,
)

nucleus.addConicConfinement(
    k=5e-3, halfOpeningAngle=half_angle_conic_confinement
)

# Add lamina Interaction
if condition in ["complete", "lamina"]:
    nucleus.addLaminaInteraction(
        radius=NUCLEUS_RADIUS,
        eLam={
            "B1": -1.034611,
            "B2": -0.928730,
            "B3": -0.741693,
            "NA": -0.667915,
        },
        subcompartments=["B1", "B2", "B3", "NA"],
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

print("Simulation ended. Closing files...", flush=True)

print("Files closed! All set!", flush=True)
