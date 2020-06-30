"""
Plot_mesh
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from mpl_toolkits.mplot3d import axes3d, Axes3D  # <-- Note the capitalization!
import subprocess

def plot_mesh(F, V):
    fig = plt.figure(figsize=(8, 8))
    ax = Axes3D(fig)
    if V.shape[1] == 2:
        ax.plot_trisurf(V[:, 0], V[:, 1], 0*V[:, 1], triangles=F, color=(0,0,0,0),edgecolor=(0,0,1), linewidth = 0.5)
        ax.view_init(90, -90)
        ax.set_axis_off()
    if V.shape[1] == 3:
        ax.plot_trisurf(V[:, 0], V[:, 1], V[:, 2], triangles=F, cmap=plt.cm.Spectral)

    # Draw now
    plt.show()

# draw mesh with unit circle
def plot_mesh_Circle(F, V):
    fig = plt.figure(figsize=(8, 8))
    ax = Axes3D(fig)
    if V.shape[1] == 2:
        ax.plot_trisurf(V[:, 0], V[:, 1], 0*V[:, 1], triangles=F, color=(0,0,0,0),edgecolor=(0,0,1), linewidth = 0.5)
        ax.view_init(90, -90)
        ax.set_axis_off()
    if V.shape[1] == 3:
        ax.plot_trisurf(V[:, 0], V[:, 1], V[:, 2], triangles=F, cmap=plt.cm.Spectral)

    theta = np.array(range(0,360))/360.0*2*np.pi
    ax.plot(np.cos(theta), np.sin(theta), 'r-')


    # Draw now

    # Create cubic bounding box to simulate equal aspect ratio
    if V.shape[1] == 2:
        Z = V[:, 0]*0
    else:
        Z = V[:, 2]
    X = V[:, 0]
    Y = V[:, 1]

    max_range = np.array([X.max() - X.min(), Y.max() - Y.min(), Z.max() - Z.min()]).max()
    Xb = 0.5 * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][0].flatten() + 0.5 * (X.max() + X.min())
    Yb = 0.5 * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][1].flatten() + 0.5 * (Y.max() + Y.min())
    Zb = 0.5 * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][2].flatten() + 0.5 * (Z.max() + Z.min())
    # Comment or uncomment following both lines to test the fake bounding box:
    for xb, yb, zb in zip(Xb, Yb, Zb):
        ax.plot([xb], [yb], [zb], 'w')

    plt.grid()
    plt.show()
    ax.view_init(elev=90, azim=0)






