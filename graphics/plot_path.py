"""
Plot_path
"""
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D  # <-- Note the capitalization!

def plot_path(F, V, path):
    fig = plt.figure()
    ax = Axes3D(fig)

    if V.shape[1] == 3:

        ax.plot_trisurf(V[:, 0], V[:, 1], 0 * V[:, 1], triangles=F, color=(0, 0, 0, 0), edgecolor=[0.2, 0.2, 0.2] , linewidth=0.5)
        ax.plot(V[path, 0], V[path, 1], V[path, 2], 'r-', linewidth=3)

    if V.shape[1] == 2:
        ax.plot_trisurf(V[:,0], V[:,1], V[:,1]*0, triangles=F,  color=(0, 0, 0, 0), edgecolor=[0.2, 0.2, 0.2] , linewidth = 0.5)
        ax.plot(V[path, 0], V[path, 1], V[path, 1]*0, 'r-', linewidth=3)
        ax.set_aspect('equal')

    ax.set_axis_off()

    # Draw now
    plt.draw()

    plt.show()
