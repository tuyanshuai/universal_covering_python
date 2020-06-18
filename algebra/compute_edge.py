"""
[edge, eif] = compute_edge(face)
"""
import numpy as np
from scipy import sparse
from algebra.compute_adjacent_matrix import compute_adjacency_matrix

def compute_edge(face):
    am, amd = compute_adjacency_matrix(face)
    IJ = np.argwhere(am)
    I = IJ[:, 1]
    J = IJ[:, 0]
    ind = I < J
    edge = np.array([I[ind], J[ind]]).transpose()

    eif = np.zeros(edge.shape)
    neg = edge.shape[0]
    for i in range(0,neg):
        eif[i, 0] = amd[edge[i, 0], edge[i, 1]]
        eif[i, 1] = amd[edge[i, 1], edge[i, 0]]

    eif = eif-1
    return edge.astype('int'), eif.astype('int')
