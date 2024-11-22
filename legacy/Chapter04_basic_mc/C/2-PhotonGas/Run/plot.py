import matplotlib.pyplot as plt
import numpy as np
import sys

# load data
if len(sys.argv) < 2:
    print("Usage: python plot.py filename")
    data = np.empty((1,2))
else:
    data = np.genfromtxt(sys.argv[1], comments="#", delimiter=" ")

    # catch empty file (except header)
    if not len(data):
        data = np.empty((1,2))

fig, ax = plt.subplots()

# set all style
ax.set_xlabel("Inverse temperature")
ax.set_ylabel("Average occupancy")
ax.grid(True)
ax.set_xlim(0.0, 5.0)
ax.set_ylim(0, 2)

# plot results and theoretical line
x = np.linspace(0.0, 5.0, 1000)
ax.plot(x, 1.0 / (np.exp(x)-1), label="theoretical, 1/(exp(x)-1)")
ax.scatter(*data.T, label="simulation", c="black", marker="x")
ax.legend()
plt.savefig("plot.png")
plt.show()