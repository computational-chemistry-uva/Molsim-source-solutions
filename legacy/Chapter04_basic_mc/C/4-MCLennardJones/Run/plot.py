import matplotlib.pyplot as plt
import numpy as np
import glob
import sys

fig, ax = plt.subplots()

# set all style
ax.set_xlabel("Density")
ax.set_ylabel("Pressure")
ax.grid(True)
#ax.set_xlim(, 5.0)
#ax.set_ylim(0, 4.0)
ax.set_xscale("log")
ax.set_yscale("log")

# plot results and theoretical line
for fname in glob.glob('scatter-*.dat'):
    data=np.loadtxt(fname, ndmin=2)
    x, y = data[:,0], data[:,1]
    plt.plot(x, y, c="blue", marker=".", markersize=1)
for fname in glob.glob('results-*.dat'):
    data=np.loadtxt(fname, ndmin=2)
    x, y = data[:,0], data[:,1]
    plt.plot(x, y, c="black", marker="x", markersize=10)
plt.savefig("plot.png")
plt.show()
