import os
import sys

import numpy as np
from OpenMiChroM.ChromDynamics import MiChroM

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

output_folder = sys.argv[1]
platform = sys.argv[2] if len(sys.argv) > 2 else "CUDA"

nucleus = MiChroM(name="chr1", temperature=1.0, timeStep=0.01)
nucleus.setup(platform=platform)


# Output folder:
nucleus.saveFolder(output_folder)

# Loading individual chromosomes:
chromosomes = nucleus.initStructure(
    mode="spring",
    ChromSeq=os.path.join(SCRIPT_DIR, "../1_inputs/chr1_subcompartments.txt"),
)

# Loading chromosomes in the simulation context
nucleus.loadStructure(chromosomes, center=True)

## Including the forces in the system

# Homopolymer Potentials
nucleus.addFENEBonds()
nucleus.addAngles()
nucleus.addSelfAvoidance()

# Chromosome Potentials
nucleus.addCustomTypes(
    mu=3.22,
    rc=1.78,
    TypesTable=os.path.join(SCRIPT_DIR, "../1_inputs/ff_compartments-and-nb.csv"),
)
nucleus.addIdealChromosome(mu=3.22, rc=1.78, dinit=3, dend=500)

# Collapse Potential
nucleus.addFlatBottomHarmonic(nRad=10)

nucleus.createSimulation()

nucleus.saveStructure(
    fileName="initial-structure", mode="pdb"
)

nucleus.createReporters(
    statistics=True, energyComponents=True, traj=False, interval=10**3
)

conditions = [
    "control",
    "lamina",
    "nucleolus",
    "complete",
    "isolation",
    "nuclear-bodies",
]
TOTAL_STEPS = (
    500_000
    + len(conditions) * 3_000 * 10**3
    + (len(conditions) - 1) * 1_000 * 10**3
)

print("Running initial collapse...", flush=True)
nucleus.run(
    nsteps=500 * 10**3,
    report=True,
    interval=10**4,
    totalSteps=TOTAL_STEPS,
)
print("Done!", flush=True)

for i, condition in enumerate(conditions):
    print(
        f"Generating initial structure for {condition} condition...",
        flush=True,
    )

    # Running the annealing to T = 2.0
    for T in np.linspace(1.1, 2.0, num=100):
        nucleus.integrator.setTemperature(T / 0.008314)
        nucleus.run(
            nsteps=10 * 10**3,
            report=True,
            interval=10**4,
            totalSteps=TOTAL_STEPS,
        )

    # Running the annealing to go to T = 0.1
    for T in np.linspace(1.99, 0.01, num=200):
        nucleus.integrator.setTemperature(T / 0.008314)
        nucleus.run(
            nsteps=10 * 10**3,
            report=True,
            interval=10**4,
            totalSteps=TOTAL_STEPS,
        )

    nucleus.saveStructure(
        fileName=f"chr1_init_{condition}", mode="pdb"
    )

    if i < len(conditions) - 1:
        # Returning temperature to T = 1.0
        for T in np.linspace(0.02, 1.0, num=100):
            nucleus.integrator.setTemperature(T / 0.008314)
            nucleus.run(
                nsteps=10 * 10**3,
                report=True,
                interval=10**4,
                totalSteps=TOTAL_STEPS,
            )
    print("Done!", flush=True)
