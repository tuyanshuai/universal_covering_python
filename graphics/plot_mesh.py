"""
Plot_mesh
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from mpl_toolkits.mplot3d import axes3d, Axes3D  # <-- Note the capitalization!
import subprocess

def plot_mesh(F, V):
    fig = plt.figure()
    ax = Axes3D(fig)
    if V.shape[1] == 2:
        ax.plot_trisurf(V[:, 0], V[:, 1], 0*V[:, 1], triangles=F, color=(0,0,0,0),edgecolor=(0,0,0), linewidth = 0.5)
        ax.view_init(90, -90)
        ax.set_axis_off()
        ax.set_aspect('equal')
    if V.shape[1] == 3:
        ax.plot_trisurf(V[:, 0], V[:, 1], V[:, 2], triangles=F, cmap=plt.cm.Spectral)

    # Draw now
    plt.show()



