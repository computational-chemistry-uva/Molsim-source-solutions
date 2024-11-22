import matplotlib.pyplot as plt
import numpy as np
import sys

# load data
rms = np.genfromtxt("rms.dat")

fig, ax = plt.subplots(1)

# set all style
ax.set_xlabel("Time")
ax.set_ylabel("Mean square displacement")
ax.grid(True)

# fit line
p = np.polyfit(rms[:, 0], np.average(rms[:, 1:], axis=1), deg=1)

# plot results and theoretical line
ax.plot(rms[:, 0], rms[:, 1], label='x')
ax.plot(rms[:, 0], rms[:, 2], label='y')
ax.plot(rms[:, 0], p[0] * rms[:, 0] + p[1], label="fit", c="black", ls="--")
ax.legend()

fig.tight_layout()
plt.savefig("plot.png")
plt.show()