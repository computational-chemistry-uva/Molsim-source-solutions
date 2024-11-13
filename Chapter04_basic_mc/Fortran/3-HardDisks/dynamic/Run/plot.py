import matplotlib.pyplot as plt
import numpy as np
import sys

# load data
data = {}
if len(sys.argv) < 2:
    print("Usage: python plot.py filename1 filename2 filename3 ...")
    sys.exit(0)
else:
    for fn in sys.argv[1:]:
        data[fn] = np.genfromtxt(fn, comments="#", delimiter=" ")

fig, ax = plt.subplots()

# set all style
ax.set_xlabel("Radial distance r")
ax.set_ylabel("Radial distribution g(r)")
ax.grid(True)
ax.set_xlim(0.0, 5.0)
ax.set_ylim(0, 4.0)

# plot results and theoretical line
for fn, d in data.items():
    ax.plot(*d.T, label=fn)
ax.legend()
plt.savefig("plot.png")
plt.show()