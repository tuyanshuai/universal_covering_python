""" disk harmonic map
% Disk harmonic map of a 3D simply-connected surface.
%
%% Syntax
%   uv = disk_harmonic_map(face,vertex)
%
%% Description
%  face  : double array, nf x 3, connectivity of mesh
%  vertex: double array, nv x 3, vertex of mesh
%
%  uv: double array, nv x 2, uv coordinates of vertex on 2D circle domain
"""
from algebra import *
from scipy.linalg import norm
from scipy import sparse
from scipy.sparse.linalg import spsolve
import math

def disk_harmonic_map(face, vertex):
    nv = vertex.shape[0]
    bd = compute_bd(face)
    db = vertex[bd] - vertex[np.append(bd[1::], bd[0])]
    bl = norm(db, axis=1)
    t = np.cumsum(bl) / np.sum(bl) * 2 * math.pi
    t = np.append(t[-1], t[0:-1])
    # use edge length to parameterize boundary
    uvbd = np.array([np.cos(t), np.sin(t)]).transpose()
    uv = np.zeros((nv,2))
    uv[bd,:] = uvbd
    inflag = np.ones((nv,1))
    inflag[bd] = 0
    inflag = np.argwhere(inflag.flatten()).flatten()
    A = laplace_beltrami(face,vertex)
    A = A.tocsc()
    Ain = A[inflag,:][:,inflag]
    rhs = -A[inflag,:][:,bd] @ uvbd
    uvin = spsolve(Ain,rhs)
    uv[inflag,:] = uvin
    return uv