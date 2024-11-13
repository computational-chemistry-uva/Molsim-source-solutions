import matplotlib.pyplot as plt
import numpy as np

msd = np.genfromtxt("msd.dat")
rdf = np.genfromtxt("rdf.dat")
vacf = np.genfromtxt("vacf.dat")

fig, ax = plt.subplots(1, 3, figsize=(12,5))


# plot msd
ax[0].set_title("Mean square displacement")
ax[0].set_xlabel("time")
ax[0].set_ylabel("msd")
ax[0].grid(True)
ax[0].plot(msd[:, 0], msd[:, 1])

# plot diff coef fit
ax_inset_msd = fig.add_axes([0.15, 0.2, 0.15, 0.15])
ax_inset_msd.plot(msd[:, 0], msd[:, 2])
ax_inset_msd.set_title("Diffusion coefficient")
ax_inset_msd.set_xlabel("time")
ax_inset_msd.set_xticks([0, 37.5])

# plot rdf
ax[1].set_title("Radial distribution function")
ax[1].set_xlabel("r")
ax[1].set_ylabel("g(r)")
ax[1].grid(True)
ax[1].plot(rdf[:, 0], rdf[:, 1])

# plot velocity autocorrelation
ax[2].set_title("Velocity autocorrelation function")
ax[2].set_xlabel("time")
ax[2].set_ylabel("vacf")
ax[2].grid(True)
ax[2].plot(vacf[:, 0], vacf[:, 1])

# plot diff coef fit
ax_inset_vacf = fig.add_axes([0.78, 0.4, 0.2, 0.2])
ax_inset_vacf.plot(vacf[:, 0], vacf[:, 2])
ax_inset_vacf.set_title("Integral <v(t)v(0)>")
ax_inset_vacf.set_xticks([0, 37.5])

fig.tight_layout()
plt.savefig("plot.png")
plt.show()