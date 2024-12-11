from _molsim import HardDisks, MonteCarlo, MolecularDynamics
from .mullerBrown import (
    mullerBrownPotential,
    mullerBrownGradient_dx,
    mullerBrownGradient_dy,
    mullerBrownHeatmap,
    mullerBrownPotentialAndGradient,
    plot_muller_brown_heatmap,
)
from .hardDisks import dynamicHardDisks, staticHardDisks, hardDiskTrajectory

import matplotlib as mpl

mpl.rcParams["font.size"] = 16
mpl.rcParams["figure.figsize"] = (8, 6)
