"""
Plot_ vertex ring
"""
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D  # <-- Note the capitalization!
import numpy as np


def plot_vertex_ring(F, V, vr):
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.view_init(90, -90)

    if V.shape[1] == 3:
        ax.plot_trisurf(V[:,0], V[:,1], V[:,2], triangles=F,  color=[0.2, 0.2, 0.2, 0] ,  edgecolor=[0,0,0], linewidth = 0.5)
        for vri in vr:
            ax.plot(V[vri, 0], V[vri, 1], V[vri, 2], 'r-', linewidth=3)


    if V.shape[1] == 2:
        ax.plot_trisurf(V[:, 0], V[:, 1], 0*V[:, 0], triangles=F, color=[1.0, 1.0, 1.0, 0.0], edgecolor=[0,0,0],linewidth=0.5, alpha=0.05)
        for i in range(len(vr)):
            cur_color =  np.random.rand(3,)
            vri = vr[i]
            ax.plot(V[list([i]), 0], V[list([i]), 1], 0*V[list([i]), 0], 'o', linewidth=2, color =cur_color)
            ax.plot(V[vri, 0], V[vri, 1], 0 * V[vri, 0], '-', linewidth=2, color=cur_color)




    ax.set_aspect('equal')


    plt.show()