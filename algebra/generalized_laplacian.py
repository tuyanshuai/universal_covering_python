## generalized laplacian
# Discretize generalized laplacian operator as a sparse matrix.
#  
# Generalized laplacian operator has the form 
# 
# \[ \nabla (D \cdot \nabla) \]
# 
# where \( D = [a,bb,c] \) and \( \mu = \rho + \imath \tau \)
# 
# \( a = \frac{(1-\rho)^2+\tau^2}{1-\rho^2-\tau^2} \),
# \( b = \frac{ -2\tau}{1-\rho^2-\tau^2} \),
# \( c = \frac{(1+\rho)^2+\tau^2}{1-\rho^2-\tau^2} \).
# 
# More details on discretization, pleas refer to [1].
# 
# # Lui, L., Lam, K., Yau, S., and Gu, X. Teichmuller Mapping (T-Map) and Its
# Applications to Landmark Matching Registration. SIAM Journal on Imaging 
# Sciences 2014 7:1, 391-426
#
## Syntax
#   A = generalized_laplacian(face,uv,mu)
#
## Description
#  face: double array, nf x 3, connectivity of mesh
#  uv  : double array, nv x 2, uv coordinate of mesh
#  mu  : complex array, nf x 1, beltrami coefficient on all faces
# 
#  A: sparse matrix, nv x nv, discretized generalized laplacian operator
from numpy import *
from algebra.face_area import *
import scipy.sparse as sp
from dbgtool.dbgtool import *

def generalized_laplacian(face, uv, mu):
    af = (1 - 2 * real(mu) + power(abs(mu), 2)) / (1.0 - power(abs(mu), 2))
    bf = -2 * imag(mu) / (1.0 - power(abs(mu),2))
    gf = (1 + 2 * real(mu) + power(abs(mu),2)) / (1.0 - power(abs(mu),2))
    fa = face_area(face, uv)

    f0 = face[:, 0]
    f1 = face[:, 1]
    f2 = face[:, 2]
    uxv0 = uv[f1, 1] - uv[f2, 1]
    uyv0 = uv[f2, 0] - uv[f1, 0]
    uxv1 = uv[f2, 1] - uv[f0, 1]
    uyv1 = uv[f0, 0] - uv[f2, 0]
    uxv2 = uv[f0, 1] - uv[f1, 1]
    uyv2 = uv[f1, 0] - uv[f0, 0]

    v00 = (af * uxv0 * uxv0 + 2 * bf * uxv0 * uyv0 + gf * uyv0 * uyv0) / fa
    v11 = (af * uxv1 * uxv1 + 2 * bf * uxv1 * uyv1 + gf * uyv1 * uyv1) / fa
    v22 = (af * uxv2 * uxv2 + 2 * bf * uxv2 * uyv2 + gf * uyv2 * uyv2) / fa

    v01 = (af * uxv1 * uxv0 + bf * uxv1 * uyv0 + bf * uxv0 * uyv1 + gf * uyv1 * uyv0) / fa
    v12 = (af * uxv2 * uxv1 + bf * uxv2 * uyv1 + bf * uxv1 * uyv2 + gf * uyv2 * uyv1) / fa
    v20 = (af * uxv0 * uxv2 + bf * uxv0 * uyv2 + bf * uxv2 * uyv0 + gf * uyv0 * uyv2) / fa

    v00[fa == 0] = 0
    v11[fa == 0] = 0
    v22[fa == 0] = 0
    v01[fa == 0] = 0
    v12[fa == 0] = 0
    v20[fa == 0] = 0

    I = concatenate((f0, f1, f2, f0, f1, f1, f2, f2, f0))
    J = concatenate((f0, f1, f2, f1, f0, f2, f1, f0, f2))
    V = concatenate((v00, v11, v22, v01, v01, v12, v12, v20, v20))
    A = sp.csr_matrix((-V, (I, J)))

    return A
