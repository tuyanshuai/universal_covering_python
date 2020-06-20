# compute uniformization metric using Ricci flow method
from algebra import *
from math import *
import numpy as np
import numpy_groupies as npg
import scipy.sparse as sp

def dump(a):
    if sp.issparse(a):
        np.savetxt("var.csv", a.todense(), delimiter=",")
    else:
        np.savetxt("var.csv", a, delimiter=",")

def hyperbolic_ricci_flow(face,vertex):
    edge, eif = compute_edge(face)
    bd = compute_bd(face)
    nv = vertex.shape[0]
    ne = edge.shape[0]
    nb = bd.shape[0]

    # hyperbolic_cosine_law
    def cosine_law(li, lj, lk):
        cs = (np.cosh(lj) * np.cosh(lk) - np.cosh(li)) / (np.sinh(lj) * np.sinh(lk))
        cs = np.arccos(cs)
        return cs


    def edge_length(e, u):
        # hyperbolic_edge_length
        el = -np.log(np.tanh(-u[e] / 2))
        return el

    def edge_weight(u):
        # hyperbolic_edge_weight
        ew = np.zeros((ne, 1))
        el = np.sum(edge_length(edge, u), axis=1)
        r = edge_length(face, u)
        w = np.sqrt(np.sinh(r[:, 0]) * np.sinh(r[:, 1]) *np.sinh(r[:,2]) / np.sinh(r[:,0]+r[:,1]+r[:,2]))
        ind = eif[:,0] >= 0
        ew[ind] = ew[ind] + w[eif[ind, 0]] / np.sinh(el[ind])
        ind = eif[:,1] >= 0
        ew[ind] = ew[ind] + w[eif[ind, 1]] / np.sinh(el[ind])
        return ew

    def set_target_curvature():
        vtk = np.zeros(shape=(nv,1))
        if nb > 0 :
            vtk[bd] = 2*pi/nb
        return vtk


    def calculate_corner_angle(u):
        ca = np.zeros(face.shape)
        r = edge_length(face, u)
        eli = r[:, 1] + r[:, 2]
        elj = r[:, 2] + r[:, 0]
        elk = r[:, 0] + r[:, 1]
        ca[:, 0] = cosine_law(eli, elj, elk)[:,0]
        ca[:, 1] = cosine_law(elj, elk, eli)[:,0]
        ca[:, 2] = cosine_law(elk, eli, elj)[:,0]
        return ca


    def calculate_vertex_curvature(u):
        vk = np.ones((nv, 1)) * pi * 2
        vk[bd] = pi
        ca = calculate_corner_angle(u)
        vk = vk - np.transpose(np.array(npg.aggregate(face.flatten(), ca.flatten()),ndmin =2))
        return vk


    def calculate_metric():

        u = np.ones(shape=(nv,1))* np.log(np.tanh(0.5))
        vtk = set_target_curvature()

        while True:
            ew = edge_weight(u)
            vk = calculate_vertex_curvature(u)

            err = abs(vk-vtk)
            print('current error is %.10f' % np.max(err))

            if np.max(err) < 1e-10  or isnan(np.max(err)):
                break

            b = vtk - vk
            I = np.concatenate((edge[:, 0], edge[:, 1]))
            J = np.concatenate((edge[:, 1], edge[:, 0]))
            V = -np.concatenate((ew, ew)).flatten()
            A = sp.csr_matrix((V, (I, J)), shape=(nv, nv))

            el = np.sum(edge_length(edge, u), axis=1)
            V2 = -np.concatenate((ew * np.cosh(el), ew * np.cosh(el))).flatten()
            A2 = sp.csr_matrix((V2, (I, J)), shape=(nv, nv))
            A = A - sp.spdiags(np.sum(A2, axis=1).flatten(), 0, nv, nv)
            x = sp.linalg.spsolve(A,b)
            u = u + x.reshape(nv,1)
        return u

    return calculate_metric().flatten()