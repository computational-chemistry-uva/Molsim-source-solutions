import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.animation import FuncAnimation
from IPython.display import HTML



# Function to compute periodic offsets
def get_periodic_offsets(boxSize):
    return [
        np.array([0, 0]),  # Original position
        np.array([boxSize, 0]),  # Right
        np.array([-boxSize, 0]),  # Left
        np.array([0, boxSize]),  # Up
        np.array([0, -boxSize]),  # Down
        np.array([boxSize, boxSize]),  # Top-right
        np.array([-boxSize, boxSize]),  # Top-left
        np.array([boxSize, -boxSize]),  # Bottom-right
        np.array([-boxSize, -boxSize]),  # Bottom-left
    ]


def hardDiskAnimation(x,boxSize, periodic_boundary_conditions, dynamic=True, pIdx=None):
    n_frames = np.shape(x)[0]
    n_particles = np.shape(x)[1]

    # Initialize the plot
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_xlim(-0.25 * boxSize, boxSize * 1.25)
    ax.set_ylim(-0.25 * boxSize, boxSize * 1.25)
    ax.axhline(0)
    ax.axhline(boxSize)
    ax.axvline(0)
    ax.axvline(boxSize)
    ax.set_aspect("equal")

    # Pre-compute offsets for periodic images
    if periodic_boundary_conditions:
        offsets = get_periodic_offsets(boxSize)
    else:
        offsets = [np.array([0, 0])]
    n_offsets = len(offsets)
    # Pre-create all Circle objects
    circle_patches = []
    for _ in range(n_particles * n_offsets):  # Total number of circles
        circle = Circle((0, 0), 0.5, edgecolor="blue", fill=False)
        circle_patches.append(circle)
        ax.add_patch(circle)

    def update(frame):
        positions = x[frame] % boxSize  # Get positions for the current frame
        modifiedIdx = [int(pIdx[frame, 0]) * len(offsets) + i for i in range(len(offsets))]
        index = 0
        for position in positions:
            for offset in offsets:
                circle_patches[index].center = position + offset
                if dynamic:
                    if index in modifiedIdx:
                        circle_patches[index].set_edgecolor("red")
                        if pIdx[frame, 1]:
                            circle_patches[index].set_fill(True)
                    else:
                        circle_patches[index].set_edgecolor("blue")
                        circle_patches[index].set_fill(False)

                index += 1


    # Create the animation
    ani = FuncAnimation(fig, update, frames=n_frames, repeat=True)
    # Prevent the static figure from displaying
    plt.close(fig)
    # Display the animation inline in Jupyter Notebook

    return HTML(ani.to_jshtml())