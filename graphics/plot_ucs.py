"""
Plot_mesh
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from mpl_toolkits.mplot3d import axes3d, Axes3D  # <-- Note the capitalization!

from algebra import *


def colorpool(i, j):
    r = (0.8, 0.2, 0.2)
    g = (0.2, 0.8, 0.2)
    b = (0.2, 0.2, 0.8)
    y = (0.8, 0.8, 0.2)
    c = (0.2, 0.8, 0.8)
    m = (0.8, 0.2, 0.8)
    s = (0.6, 0.6, 0.6)
    p = (1, 0.6, 0.8)

    cl = list([s, r, g, b, y, c, m])
    ci = (i + j) % 7
    return cl[ci]


def plot_ucs(face, ucs):
    fig = plt.figure(figsize=(8, 8))
    ax = Axes3D(fig)

    bd = compute_bd(face)

    for i in range(len(ucs)):
        ui = ucs[i]
        for j in range(len(ui)):
            if len(ui[j].shape) == 1:  # complex number
                uij = ui[j]
                nv = uij.shape[0]
                vertex = np.zeros((nv, 3))
                vertex[:, 0] = real(uij)
                vertex[:, 1] = imag(uij)
                vertex[:, 2] = real(uij) * 0
            else:
                vertex = ui[j]

            ax.plot_trisurf(vertex[:, 0], vertex[:, 1], vertex[:, 2], triangles=face,
                            edgecolor=colorpool(i, j), color=colorpool(i, j),
                            linewidth=0.1, antialiased=False)

    # Draw now
    ax.view_init(elev=90, azim=0)
    theta = np.array(range(0, 360)) / 360.0 * 2 * np.pi
    ax.plot(np.cos(theta), np.sin(theta), 'r-')
    ax.set_axis_off()
    plt.show()
