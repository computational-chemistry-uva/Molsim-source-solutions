import matplotlib.pyplot as plt
import numpy as np
import sys

# load data
data = np.genfromtxt("results.dat")

fig, ax = plt.subplots()

# set all style
ax.set_xlabel("Energy levels")
ax.set_ylabel("Distribution")
ax.grid(True)

# plot results and theoretical line
ax.plot(*data.T)
ax.legend()
plt.savefig("plot.png")
plt.show()