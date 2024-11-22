from md import MolecularDynamics
import matplotlib.pyplot as plt
import numpy as np

md = MolecularDynamics(
    numberOfParticles=100, 
    temperature=0.5, 
    dt=0.005, 
    boxSize=5.0, 
    logLevel=0,
    seed=12
)

# Running the simulation to equilibrate and then to gather data.
md.run(5000, equilibrate=True, outputPDB=False)
md.run(50000, equilibrate=False, outputPDB=False)