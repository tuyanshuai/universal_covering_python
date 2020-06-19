# slice_mesh
# Slice mesh open along a collection of edges ee, which usually comes from 
# cut_graph(directly), or compute_greedy_homotopy_basis and
# compute_homology_basis (need to form edges from closed loops in basis).
# ee can form a single closed loops or multiple closed loops. 
# 
# Syntax
#   [face_new,vertex_new,father] = slice_mesh(face,vertex,ee)
#
# Description
#  face  : double array, nf x 3, connectivity of mesh
#  vertex: double array, nv x 3, vertex of mesh
#  ee    : double array, n x 2, a collection of edges, each row is an edge on 
#          mesh, may not be in consecutive order. 
# 
#  face_new  : double array, nf x 3, connectivity of new mesh after slice
#  vertex_new: double array, nv' x 3, vertex of new mesh, vertex number is
#              more than original mesh, since slice mesh will separate each
#              vertex on ee to two vertices or more.
#  father    : double array, nv' x 1, father indicates the vertex on original
#              mesh that new vertex comes from.
#
from algebra import *


def slice_mesh(face, vertex, ee):
    nv = vertex.shape[0]
    _, amd = compute_adjacency_matrix(face)

    if isinstance(ee[0], list):
        ees = list([list(), list()])
        for ei in ee:
            ees[0] = ees[0] + ei
            ees[1] = ees[1] + ei[1::] + [ei[0]]
        ee = ees
    ee = np.array(ee).transpose()

    G = sp.csr_matrix((np.ones((ee.shape[0], 1)).flatten(), (ee[:, 0], ee[:, 1])), shape=(nv, nv))
    G = G + G.transpose()

    ev = np.unique(ee.flatten())

    vre = compute_vertex_ring(face, vertex, ev, ordered=True)

    face_new = face.copy()
    vertex2 = np.zeros((ee.shape[0] * 2, 3))
    father2 = np.zeros((ee.shape[0] * 2, 1))
    k = 0
    for i in range(ev.shape[0]):
        evr = vre[i]
        for i0 in range(len(evr)):
            if G[evr[i0], ev[i]]:
                break

        if evr[0] == evr[-1]:  # interior point
            evr = evr[i0:-1] + evr[::i0 + 1]
        else:
            evr = evr[i0:-1] + evr[::i0]

        for j in range(1, len(evr)):
            fi = amd[evr[j], ev[i]] - 1
            if fi >= 0:
                fij = face_new[fi, :] == ev[i]
                face_new[fi, fij] = nv + k

            if G[ev[i], evr[j]]:
                vertex2[k, :] = vertex[ev[i], :]
                father2[k] = ev[i]
                k = k + 1

    vertex_new = np.concatenate((vertex, vertex2), axis=0)
    father = np.array(range(nv))
    father = np.concatenate((father, father2), axis=0)

    fu = np.unique(face_new)
    index = np.zeros((np.max(fu), 1))
    index[fu] = np.array(range(fu.shape[0]))
    face_new = index[face_new]
    vertex_new = vertex_new[fu, :]
    father = father[fu]
    return face_new, vertex_new, father
