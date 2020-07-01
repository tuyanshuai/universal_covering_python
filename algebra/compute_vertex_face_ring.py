import numpy as np
from algebra.compute_vertex_ring import *
from algebra.compute_halfedge import *


# TODO Check
def compute_vertex_face_ring(face, vertex=np.array([]), vc=np.array([]), ordered=None):
    nv = np.max(np.max(face)) + 1
    if vc.size == 0:
        vc = np.array(range(nv))
    if ordered == None:
        ordered = False

    vr = compute_vertex_ring(face, vertex, vc, ordered)
    he, heif = compute_halfedge(face)

    eifs = sp.csr_matrix((heif+1, (he[:, 0], he[:, 1])))

    vfr = [list()] * nv

    for i in range(nv):
        vri = vr[i]
        for vj in vri:
            if eifs[i,vj] >0:
                vfr[i].append(eifs[i,vj]-1)
    return vfr
