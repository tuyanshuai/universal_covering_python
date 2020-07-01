# # disk area - preserving  mapping 
# # Contribution
# Author: Yanshuai  Tu
# Created: 2018 / 05 / 07
# Revised:
# CIDSE, ASU, http: // gsl.lab.asu.ed

from algebra.face_area import *
from algebra.compute_bd import *
from algebra.compute_edge import *
import scipy.sparse as sp

def dot(A,B, axis = 0):
    return np.sum(A.conj() * B, axis=axis)

# # Laplace with area on
def simLap(face, vertex, f):
    edge, eif = compute_edge(face)
    ne = edge.shape[0]
    ew = np.zeros(ne)
    ind = eif[:, 0] > 0
    ev1 = np.sum(face[eif[ind, 0],:], axis=1) - np.sum(edge[ind,:], axis=1)

    fi = f[edge[ind, 0],:]
    fj = f[edge[ind, 1],:]
    fk = f[ev1,:]

    vi = vertex[edge[ind, 0],:]
    vj = vertex[edge[ind, 1],:]
    vk = vertex[ev1,:]

    ew[ind] = ew[ind] + dot(fi - fk, fj - fk, axis=1) / (2.0 * triangle_area(vi, vj, vk))
    ind = eif[:, 1] > 0
    ev2 = np.sum(face[eif[ind, 1],:], axis=1) - np.sum(edge[ind,:], axis=1)

    fi = f[edge[ind, 0],:]
    fj = f[edge[ind, 1],:]
    fk = f[ev2,:]

    vi = vertex[edge[ind, 0],:]
    vj = vertex[edge[ind, 1],:]
    vk = vertex[ev2,:]


    ew[ind] = ew[ind] + dot(fi - fk, fj - fk, axis=1) / (2.0 * triangle_area(vi, vj, vk))

    I = np.concatenate((edge[:, 0], edge[:, 1]))
    J = np.concatenate((edge[:, 1], edge[:, 0]))
    V=  np.concatenate((ew/2.0, ew/2.0))
    A = sp.csc_matrix((V,(I,J)))
    sA = np.sum(A, axis = 1)
    A = A - np.diag(sA)
    return A


def disk_area_mapping(face, vertex, uv):



    # sum of surface face area
    # sumfa = np.sum(face_area(face, vertex))
    B = compute_bd(face)

    # set landmark and desired location

    # B: Boundary inds
    # I: Interier
    # N: size(f, 1)
    N = np.array(range(uv.shape[0]))
    I = np.setdiff1d(N, B)

    uv_new = np.copy(uv)
    L = simLap(face, vertex, uv)

    for k in range(20):

        f1 = np.copy(uv_new)
        fbd = sp.linalg.spsolve(-L[B, B], L[B, I].dot(uv_new[I,:]))
        fbd = fbd - np.mean(fbd, axis=1)

        dl = np.sqrt(np.linalg.norm(fbd, axis=1))
        uv_new[B,0] =  fbd[:,0] / dl
        uv_new[B, 1] = fbd[:, 1] / dl

        fin = sp.linalg.spsolve(-L[I, I], L[I, B].dot(uv_new[B, :]))
        uv_new[I,:] = fin
        L = simLap(face, vertex, uv_new)

        maxdf = np.max(np.linalg.norm(uv_new - f1, axis=1))
        print('maxdf = %f\n' % maxdf)
        if maxdf < 1e-6:
            break

    return uv_new