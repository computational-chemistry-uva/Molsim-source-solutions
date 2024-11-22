import matplotlib.pyplot as plt
import numpy as np
import sys

# load data
dist = np.genfromtxt("results.dat")
rms = np.genfromtxt("rms.dat")

fig, ax = plt.subplots(2)

# set all style
ax[0].set_xlabel("Position relative to starting point")
ax[0].set_ylabel("Distribution")
ax[0].grid(True)
ax[1].set_xlabel("Time")
ax[1].set_ylabel("Mean square displacement")
ax[1].grid(True)

# plot results and theoretical line
ax[0].plot(*dist[:, :2].T)
ax[1].plot(*rms.T)
fig.tight_layout()
plt.savefig("plot.png")
plt.show()