## linear beltrami solver
# Solve Beltrami equation from given \(\mu\) and some landmark constraints.
#  
# We can solve Beltrami equation \( \frac{\partial{f}}{\partial{\bar{z}}} = \mu 
# \frac{\partial{f}}{\partial{z}} \) with given \(\mu\), to get a mapping f. 
# Together with compute_bc, which computes mu from f, we can see that f and 
# mu can mutually determined, or f can be represented by mu.
#
## Syntax
#   [uv_new,mu_new] = linear_beltrami_solver(face,uv,mu)
#   [uv_new,mu_new] = linear_beltrami_solver(face,uv,mu,landmark,target)
#
## Description
#  face: double array, nf x 3, connectivity of mesh
#  uv  : double array, nv x 2, uv coordinate of mesh
#  mu  : complex array, nf x 1, beltrami coefficient on all faces
#  landmark: double array, n x 1, index landmark points to be fixed
#  target  : double array, n x 2, target coordinates of landmark points
#
#  uv_new: double array, nv x 2, uv coordinate of mesh
#  mu_new: complex array, nf x 1, beltrami coefficient of the map uv to
#          uv_new, ideally should be identical to mu, but in discrete case, there
#          may have slight difference. It can be significant with presence
#          of landmark constraints.
#


from algebra.generalized_laplacian import *
from algebra.compute_bc import *

import scipy.sparse as sp
import numpy as np
from dbgtool.dbgtool import *


def linear_beltrami_solver(face, uv, mu, landmark, target):
    A = generalized_laplacian(face, uv, mu)
    b = -A[:, landmark] * target
    b[landmark, :] = target
    A[landmark, :] = 0
    A[:, landmark] = 0

    A = A + sp.csr_matrix((np.ones(landmark.shape[0]), (landmark, landmark)), shape=A.shape)
    uv_new = sp.linalg.spsolve(A, b)
    mu_new = compute_bc(face, uv, uv_new)
    return uv_new, mu_new
