from parameterization.disk_harmonic_map import *
from algebra.compute_bc import *
from algebra.compute_bd import *
from algebra.linear_beltrami_solver import *
from algebra.compute_vertex_face_ring import *
import numpy as np

"""
% Compute the disk conformal mapping using the method in [1].
%
% Input:
% vertex: nv x 3 vertex coordinates of a simply-connected open triangle mesh
% face: nf x 3 triangulations of a simply-connected open triangle mesh
% 
% Output:
% map: nv x 2 vertex coordinates of the disk conformal parameterization
% 
% Remark:
% 1. Please make sure that the input mesh does not contain any 
%    unreferenced vertices/non-manifold vertices/non-manifold edges.
% 2. Please remove all valence 1 boundary vertices (i.e. vertices with 
%    only 1 face attached to them) before running the program.
% 
% If you use this code in your own work, please cite the following paper:
% [1] P. T. Choi and L. M. Lui, 
%     "Fast Disk Conformal Parameterization of Simply-Connected Open Surfaces."
%     Journal of Scientific Computing, 65(3), pp. 1065-1090, 2015.
"""

import numpy as np


def f2v(face, vertex, mu):
    vr = compute_vertex_ring(face)
    nv = vertex.shape[0]
    muv = np.zeros(nv, np.complex)
    for i in range(nv):
        vri = vr[i]
        muv[i] = np.sum(mu[vri]) / len(vri)
    return muv


def disk_conformal_map(face, vertex):
    parameter_north = 5
    parameter_south = 100
    parameter_threshold = 0.00001

    map = disk_harmonic_map(face, vertex)

    bdy_index = compute_bd(face)

    z = map[:, 0] + 1j * map[:, 1]
    print("Initialization completed.\n")

    # %% North     Pole    iteration
    mu = compute_bc(face, map, vertex)
    mu_v = f2v(face, vertex, mu)
    bdy_index_temp = np.concatenate((bdy_index[1::], bdy_index[0]))

    least = np.argmin(np.abs(mu_v[bdy_index]) + np.abs(mu_v[bdy_index_temp]))

    mi = np.mod(least+1, bdy_index.shape[0])

    ang1 = np.angle(z[bdy_index[least]])
    ang2 = np.angle(z[bdy_index[mi]])
    z = z*np.exp(-1j * (ang1 + ang2)/ 2.0)
    g = 1j * (1 + z) / (1 - z)
    ind = np.argsort(-np.real(z))

    maxind = np.max(np.round(vertex.shape[0] / parameter_north), np.min(100, z.shape[0])).astype('int')
    indmaxmin = ind[::maxind]
    fixed1 = np.setdiff(indmaxmin ,bdy_index)
    fixed2 = np.where(np.real(g) == np.max(np.real(g)))
    fixed3 = np.where(np.real(g) == np.min(np.real(g)))
    fixed = np.concatenate((fixed1, fixed2, fixed3))

    # TODO : MATLAB 47-48 Line
    P = [real(g),imag(g),ones(length(g),1)];
    mu = beltrami_coefficient(P, face, vertex);

#
# def disk_conformal_map(face, vertex):
#     uv = disk_harmonic_map(face, vertex)
#     uv_new = np.copy(uv)
#     bd = compute_bd(face)
#     nf = face.shape[0]
#     mu = np.zeros(nf, np.complex)
#
#     for i in range(20):
#         uv_new, _ = linear_beltrami_solver(face, uv, mu, bd, uv[bd, :])
#         mu_new = compute_bc(face, uv_new, vertex)
#         dmu = np.max(np.abs(mu - mu_new))
#         print('# %d current change of dmu is %f' % (i, dmu))
#         if dmu < 1e-2:
#             break
#         mumean = np.mean(np.abs(mu_new))
#         mu = mumean * mu_new / np.abs(mu_new)
#
#
#     return uv_new
