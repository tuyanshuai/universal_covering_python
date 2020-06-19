"""

function vr = compute_vertex_ring(face,vertex,vc,ordered)

"""
import numpy as np
from algebra import *
from dbgtool.dbgtool import *

from scipy.spatial import ConvexHull


# find next point along the ccw
def find_next_point(faces, eg):
    for fi in faces:
        for j in range(3):
            if fi[j] == eg[0] and fi[(j + 1) % 3] == eg[1]:
                return fi[(j + 2) % 3]
    return None


def compute_vertex_ring(face, vertex, vc=np.array([]), ordered=None):
    nv = np.max(np.max(face)) + 1
    if vc.size == 0:
        vc = np.array(range(nv))
    if ordered == None:
        ordered = False

    vr = [[] for i in range(vc.size)]
    bd = compute_bd(face)
    isbd = np.zeros((nv,)).astype(int)
    isbd[bd] = 1

    if not ordered:
        # if not ordered:
        am, _ = compute_adjacency_matrix(face)
        IJ = np.argwhere(am[:, vc])
        I = IJ[:, 0]
        J = IJ[:, 1]
        ijsid = np.argsort(I + J * nv)
        I = I[ijsid]
        J = J[ijsid]
        for i in range(J.size):
            vr[J[i]].append(I[i])

    # order the neighbour
    if ordered:
        vvif, nvif, pvif = compute_connectivity(face)
        for i in range(vc.size):
            fs = vvif[vc[i], :]
            v1 = np.argwhere(fs)[0][1]

            if isbd[vc[i]]:
                while vvif[v1, vc[i]]:
                    f2 = vvif[v1, vc[i]] - 1
                    v1 = pvif[f2, v1] - 1

            vi = list([v1])
            v0 = v1

            while vvif[vc[i], v1]:
                f1 = vvif[vc[i], v1] - 1
                v1 = nvif[f1, v1] - 1
                vi.append(v1)

                if v0 == v1:
                    break
            vr[i] = vi

    return vr
