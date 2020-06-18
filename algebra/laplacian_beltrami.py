"""Laplace Beltrami operator on the mesh.
L = laplace_beltrami(face,vertex)
"""
import numpy as np
from algebra.compute_edge import *
from scipy.linalg import norm
from scipy import sparse
def cot2(pi, pj, pk):
    a = norm(pj - pk, axis=1)
    b = norm(pk - pi, axis=1)
    c = norm(pi - pj, axis=1)
    cs =np.divide((np.multiply(b,b) + np.multiply(c,c) - np.multiply(a,a)) , (2*np.multiply(b,c)) )
    ss2 = 1 - np.multiply(cs,cs)
    ss2[ss2 < 0] = 0
    ss2[ss2 > 1] = 1
    ss = np.sqrt(ss2)
    ct = np.divide(cs, ss)
    return ct


def laplace_beltrami(face,vertex):
    edge, eif = compute_edge(face)
    ne = edge.shape[0]
    ew = np.zeros((ne, 1))
    ind = eif[:,0] >= 0

    ev1 = np.sum(face[eif[ind,0],:], axis=1) - np.sum(edge[ind,:],axis = 1)

    ct1 = cot2(vertex[ev1,:], vertex[edge[ind, 0],:], vertex[edge[ind, 1],:])
    ew[ind] = (ew[ind].transpose() + ct1).transpose()
    ind = eif[:, 1] >= 0
    ev2 = np.sum(face[eif[ind, 1],:], axis=1) - np.sum(edge[ind,:], axis = 1)
    ct2 = cot2(vertex[ev2,:], vertex[edge[ind, 0],:], vertex[edge[ind, 1],:])
    ew[ind] = (ew[ind].transpose() + ct2).transpose()

    I = np.array([edge[:,0], edge[:,1] ])
    J = np.array([edge[:,1], edge[:,0] ])
    V = np.array([ew/2, ew/2])

    A = sparse.coo_matrix((V.flatten(), (I.flatten(), J.flatten())))
    A = A.tolil()

    sA = np.sum(A, axis=1)
    for i in range(A.shape[0]):
        A[i,i] -= sA[i,0]
    return A

