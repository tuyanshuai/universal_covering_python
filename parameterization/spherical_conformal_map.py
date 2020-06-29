"""
spherical conformal map
Spherical Conformal Map of a closed genus-0 surface. Methodology details
please refer to [1]. Some code is derived from Gary Choi.

K.C. Lam, P.T. Choi and L.M. Lui, FLASH: Fast landmark aligned
  spherical harmonic parameterization for genus-0 closed brain surfaces,
  UCLA CAM Report, ftp://ftp.math.ucla.edu/pub/camreport/cam13-79.pdf?

 Syntax
   uvw = spherical_conformal_map(face,vertex)

Description
  face  : double array, nf x 3, connectivity of mesh
  vertex: double array, nv x 3, vertex of mesh

  uvw: double array, nv x 3, spherical uvw coordinates of vertex on 3D unit sphere


Init the checking by yanshuai 06/28/2020
"""
import numpy as np
from algebra.laplacian_beltrami import *
from algebra.compute_bc import *
from algebra.linear_beltrami_solver import *
import scipy.sparse as sp
from dbgtool.dbgtool import *



def spherical_conformal_map(face, vertex):
    dv1 = vertex[face[:, 1], :] - vertex[face[:, 2], :]
    e1 = np.linalg.norm(dv1, axis=1)

    dv2 = vertex[face[:, 2], :] - vertex[face[:, 0], :]
    e2 = np.linalg.norm(dv2, axis=1)

    dv3 = vertex[face[:, 0], :] - vertex[face[:, 1], :]
    e3 = np.linalg.norm(dv3, axis=1)

    e123 = e1 + e2 + e3
    regularity = np.abs(e1 / e123 - 1 / 3.0) + np.abs(e2 / e123 - 1 / 3.0) + np.abs(e3 / e123 - 1 / 3.0)

    # choose vertex bi as the big triangle
    bi = np.argmin(regularity, axis=0)
    nv = vertex.shape[0]
    A = laplace_beltrami(face, vertex)
    fi = face[bi, :]
    Afi = A[fi, :]
    IJ = np.argwhere(Afi)
    I = IJ[:, 0]
    J = IJ[:, 1]
    V = Afi[I, J]

    # A = A - sparse(fi(I), J, V, nv, nv) + sparse(fi, fi, [1, 1, 1], nv, nv);
    A = A - sp.csr_matrix((V.toarray()[0], (fi[I], J)), shape=(nv, nv)) + sp.csr_matrix((np.array([1, 1, 1]), (fi, fi)),
                                                                                        shape=(nv, nv))

    # Set boundary condition for big triangle
    x1 = 0.0
    y1 = 0.0
    x2 = 100.0
    y2 = 0.0
    a = vertex[fi[1], :] - vertex[fi[0], :]

    b = vertex[fi[2], :] - vertex[fi[0], :]

    ratio = np.linalg.norm([x1 - x2, y1 - y2]) / np.linalg.norm(a)
    y3 = np.linalg.norm(list([x1 - x2, y1 - y2])) * np.linalg.norm(np.cross(a, b)) / np.power(np.linalg.norm(a), 2)
    x3 = np.sqrt(np.power(np.linalg.norm(b), 2) * ratio * ratio - np.power(y3, 2))

    # Solve matrix equation
    d = np.zeros((nv, 2))
    d[fi[0], :] = (x1, y1)
    d[fi[1], :] = (x2, y2)
    d[fi[2], :] = (x3, y3)

    uv = sp.linalg.spsolve(A, d)
    z = uv[:, 0] + 1j * uv[:, 1]
    z = z - np.mean(z)
    dz2 = np.power(np.abs(z), 2)

    vertex_new = np.zeros((nv, 3))
    vertex_new[:, 0] = 2 * np.real(z) / (1.0 + dz2)
    vertex_new[:, 1] = 2 * np.imag(z) / (1.0 + dz2)
    vertex_new[:, 2] = (-1.0 + dz2) / (1 + dz2)

    # Find optimal big triangle size
    # Reason: the distribution will be the best
    # if the southmost triangle has similar size of the northmost one
    w = np.zeros(nv, np.complex)
    w.real = vertex_new[:, 0] / (1.0 + vertex_new[:, 2])
    w.imag = vertex_new[:, 1] / (1.0 + vertex_new[:, 2])

    index = np.argsort(np.abs(z[face[:, 0]]) + np.abs(z[face[:, 1]]) + np.abs(z[face[:, 2]]))

    ni = index[1]  # since index(1) must be bi, not really
    ns = np.sum(np.abs(z[face[bi, [0, 1, 2]]] - z[face[bi, [1, 2, 0]]])) / 3

    ss = np.sum(abs(w[face[ni, [0, 1, 2]]] - w[face[ni, [1, 2, 0]]])) / 3
    z = z * (np.sqrt(ns * ss)) / ns
    dz2 = np.power(np.abs(z), 2)

    vertex_new[:, 0] = 2 * np.real(z) / (1.0 + dz2)
    vertex_new[:, 1] = 2 * np.imag(z) / (1.0 + dz2)
    vertex_new[:, 2] = (-1.0 + dz2) / (1 + dz2)

    # south pole stereographic projection
    uv[:, 0] = vertex_new[:, 0] / (1.0 + vertex_new[:, 2])
    uv[:, 1] = vertex_new[:, 1] / (1.0 + vertex_new[:, 2])
    mu = compute_bc(face, uv, vertex)

    # find the south pole
    ind = np.argsort(vertex_new[:, 2])
    fixed = ind[0:min([nv, 50])]
    # reconstruct map with given mu and some fixed point
    fuv, fmu = linear_beltrami_solver(face, uv, mu, fixed, uv[fixed, :])

    dfz2 = fuv[:, 0] * fuv[:, 0] + fuv[:, 1] * fuv[:, 1]

    uvw = np.zeros((nv, 3))
    uvw[:, 0] = 2 * fuv[:, 0] / (1.0 + dfz2)
    uvw[:, 1] = 2 * fuv[:, 1] / (1.0 + dfz2)
    uvw[:, 2] = -(-1 + dfz2) / (1 + dfz2)
    return uvw
