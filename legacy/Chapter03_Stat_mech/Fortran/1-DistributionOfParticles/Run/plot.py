import matplotlib.pyplot as plt
import numpy as np
import sys

def read_multiple_runs(fn):
    # reads multiple results in one file separated by whitespace
    # returns array of shape (n, k, c) for n results with k lines and c columns
    with open(fn, "r") as f:
        # split on whiteline
        data = f.read().split("\n\n")[:-1]
    # split on newline then whitespace
    data = np.array([[d.split(" ") for d in entry.split("\n")] for entry in data])
    return data.astype(float)

# load data
data = read_multiple_runs("results.dat")
analytical = np.genfromtxt("analytical.dat")

fig, ax = plt.subplots()

# set all style
ax.set_xlabel("Number of particles")
ax.set_ylabel("Distribution")
ax.grid(True)

# plot results and theoretical line
for i, d in enumerate(data):
    ax.plot(*d.T, label=f"Box {i}")
ax.plot(*analytical.T, label="analytical", c="black", ls="--")
ax.legend()
plt.savefig("plot.png")
plt.show()