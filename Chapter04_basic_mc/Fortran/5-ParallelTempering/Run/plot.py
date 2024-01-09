import matplotlib.pyplot as plt
import numpy as np
import sys

with open("energy_distribution.dat") as f:
    energy_data = f.read().split("\n\n")[:-1]
energy_data = [np.fromstring(d, sep=' ').reshape(-1,2) for d in energy_data]

with open("position_distribution.dat") as f:
    position_data = f.read().split("\n\n")[:-1]
position_data = [np.fromstring(d, sep=' ').reshape(-1,2) for d in position_data]

fig, ax = plt.subplots(2)

# fill temperatures manually!
if len(sys.argv) != len(energy_data) + 1:
    print("Usage: python plot.py 0.1 0.2 0.3 ... , for all temperatures")
    sys.exit(0)
temperatures = sys.argv[1:]

labels = [f"T={t}" for t in temperatures]
for i, d in enumerate(position_data):
    ax[0].plot(*d.T, label=labels[i])
ax[0].set_xlabel("Fractional box position x")
ax[0].set_ylabel("Probability")
ax[0].set_xlim(0,1)
ax[0].set_ylim(0, 0.02)
ax[0].legend()
for i, d in enumerate(energy_data):
    ax[1].plot(*d.T, label=labels[i])
ax[1].set_xlabel("Energy")
ax[1].set_ylabel("Probability")
ax[1].legend()
ax[1].set_xlim(0.0, 1.0)
plt.show()