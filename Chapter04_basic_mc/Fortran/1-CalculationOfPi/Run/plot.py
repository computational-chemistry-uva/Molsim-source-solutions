import matplotlib.pyplot as plt
import numpy as np
import sys

# load data
if len(sys.argv) < 2:
    print("Usage: python plot.py filename")
    data = np.empty((1,2))
else:
    data = np.loadtxt(sys.argv[1])

fig, ax = plt.subplots()

# set all style
ax.set_xlabel("Number of trials (x1000)")
ax.set_ylabel("Relative error")
ax.grid(True)
ax.set_xlim(90, 110000)
ax.set_ylim(1e-5, 1e-2)
ax.set_xscale("log")
ax.set_yscale("log")

# plot results and theoretical line
x = np.linspace(90, 110000)
ax.plot(x, 0.013/np.sqrt(x), label="1/sqrt(n)")
ax.scatter(*data.T, label="error", c="black", marker="x")
ax.legend()
plt.show()