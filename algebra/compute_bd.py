"""
Compute the boundary given the face

bd = compute_bd(face)
face: double array, nf x 3, connectivity of mesh
bd: double array, n x 1, consecutive boundary vertex list in ccw order
or dist n x 1, each cell is one closed boundary
"""
import numpy as np
from algebra.compute_adjacent_matrix import compute_adjacency_matrix


def compute_bd(face):
    # first compute the adjacent matrix
    am, amd = compute_adjacency_matrix(face)
    md = am - (amd > 0) * 2
    I = np.argwhere(md.transpose() == -1)[:,1] # boundary eij
    Ii = np.argsort(I)
    nbd = I.shape[0]
    bd = np.zeros(nbd, dtype=int)
    k = 0
    for i in range(0, nbd):
        bd[i] = I[k]
        k = Ii[k]
    return bd.astype(int)