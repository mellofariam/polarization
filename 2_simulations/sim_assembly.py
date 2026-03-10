"""
This script positions the already collapsed chromosome in the nucleus.

It also includes the nuclear bodies in the simulation.

Here we will include only the nucleolus and the chromatin.

The nucleolus is positioned in the center of the nucleus.
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
            str(replicaID),
            f"chr1_init_{condition}.pdb",
        )
    ],
)

chromosomes = nb.fix_chromatin(
    chromosomes, nucleus_radius=NUCLEUS_RADIUS, thresh=1.15
)
print("zmax:", np.max(chromosomes[:, 2]))
print("zmin:", np.min(chromosomes[:, 2]))

nucleolus_xyz = nb.position_nucleoli(
    chromatin=chromosomes,
    nucleoli_radius=NUCLEOLI_RADIUS,
    nucleus_radius=NUCLEUS_RADIUS,
    thresh=1.15,
)

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
)  # chromosome chain

# Nuclear Bodies Potentials

## Nucleolus
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

# Add spherical confinement
nucleus.addSphericalConfinementTanh(radius=NUCLEUS_RADIUS)

# Add lamina Interaction
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

print("Simulation created!", flush=True)

nucleus.saveStructure(fileName="setup", mode="pdb")

print("Starting simulation...", flush=True)
print("\tInitial equilibration...", flush=True)
nucleus.run(
    nsteps=500 * 10**3,
    report=True,
    interval=10**4,
    totalSteps=3_500 * 10**3,
)
print("\tDone!", flush=True)

current_nucleoli_position = nucleolus_xyz[0, 2]
nucleoli_chain_info = nucleus.chains[nucleoli_chains[0]]

Z_TARGET = 0
dz = (Z_TARGET - current_nucleoli_position) / 99

print("\tMoving nucleolus...", flush=True)

for i in range(100):
    nucleus.run(
        nsteps=5 * 10**3,
        report=True,
        interval=10**4,
        totalSteps=3_500 * 10**3,
    )
    current_positions = nucleus.getPositions()
    print(current_positions[-1], flush=True)
    current_positions[
        nucleoli_chain_info[0] : nucleoli_chain_info[1] + 1, 2
    ] += dz

    if i < 99:
        nucleus.context.setPositions(current_positions)

nucleus.saveStructure(fileName="nucleolus-moved", mode="pdb")

print("\tNucleolus moved! Forces after movement:", flush=True)
nucleus.printForces()

## Running the annealing to and back from T = 2.0

print("\tAnnealing...", flush=True)

for T in np.linspace(1.1, 2.0, num=100):
    nucleus.integrator.setTemperature(T / 0.008314)
    nucleus.run(
        nsteps=10 * 10**3,
        report=True,
        interval=10**4,
        totalSteps=3_500 * 10**3,
    )

## Running the annealing to return to T = 1.0

for T in np.linspace(1.99, 1.00, num=100):
    nucleus.integrator.setTemperature(T / 0.008314)
    nucleus.run(
        nsteps=10 * 10**3,
        report=True,
        interval=10**4,
        totalSteps=3_500 * 10**3,
    )

print("\tDone!", flush=True)
print("\tFinal equilibration...", flush=True)

nucleus.run(
    nsteps=500 * 10**3,
    report=True,
    interval=10**4,
    totalSteps=3_500 * 10**3,
)
print("\tDone!", flush=True)

print("Saving files...", flush=True, end=" ")
nucleus.saveStructure(fileName="initial-state", mode="pdb")
print("Done!", flush=True)

print("Forces in the end:", flush=True)
nucleus.printForces()

print("Simulation ended. Closing files...", flush=True)

print("Files closed! All set!", flush=True)
