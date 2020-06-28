## compute bc
# Compute Beltrami coefficients mu of mapping from uv to vertex, where vertex 
# can be 2D or 3D.
# 
# mu from 2D to 2D is defined by the Beltrami equation:
# 
# \[ \frac{\partial{f}}{\partial{\bar{z}}} = \mu \frac{\partial{f}}{\partial{z}} \]
# 
# mu from 2D to 3D is defined by: 
# 
# \[ \mu = \frac{E-G +2iF}{E+G +2\sqrt{EG-F^2}} \]
# 
# where \( ds^2=Edx^2+2Fdxdy+Gdy^2 \) is metric.
#
## Syntax
#   mu = compute_bc(face,uv,vertex)
#
from algebra.face_area import *


def compute_bc(face, uv, vertex):
    nf = face.shape[0]
    fa = face_area(face, uv)
    Duv = uv[face[:, [2, 0, 1]], :] - uv[face[:, [1, 2, 0]], :]
    Duv[:, 0] = Duv[:, 0] / np.concatenate((fa, fa, fa)) / 2.0
    Duv[:, 1] = Duv[:, 1] / np.concatenate((fa, fa, fa)) / 2.0

    if vertex.shape[1] == 2:
        nv = vertex.shape[0]
        z = np.zeros(nv, np.complex)
        z.real = vertex[:, 0]
        z.imag = vertex[:, 1]

        Dcz = np.sum(reshape((Duv[:, 1] - 1j * Duv[:, 0]) * z[face, :], (nf, 3)), axis=1)
        Dzz = np.sum(reshape((Duv[:, 1] + 1j * Duv[:, 0]) * z[face, :], (nf, 3)), axis=1)
        mu = Dcz / Dzz

    if vertex.shape[1] == 3:
        du = np.zeros((nf, 3))
        du[:, 0] = np.sum(reshape(Duv[:, 1] * vertex[face, 0], (nf, 3)), axis=1)
        du[:, 1] = np.sum(reshape(Duv[:, 1] * vertex[face, 1], (nf, 3)), axis=1)
        du[:, 2] = np.sum(reshape(Duv[:, 1] * vertex[face, 2], (nf, 3)), axis=1)
        dv = np.zeros((nf, 3))
        dv[:, 0] = np.sum(reshape(Duv[:, 0] * vertex[face, 0], (nf, 3)), axis=1)
        dv[:, 1] = np.sum(reshape(Duv[:, 0] * vertex[face, 1], (nf, 3)), axis=1)
        dv[:, 2] = np.sum(reshape(Duv[:, 0] * vertex[face, 2], (nf, 3)), axis=1)

        E = np.sum(du * du, axis=1)
        G = np.sum(dv * dv, axis=1)
        F = -np.sum(du * dv, axis=1)
        mu = (E - G + 2j * F) / (E + G + 2 * np.sqrt(E * G - F ^ 2))

    if vertex.shape[1] != 2 and vertex.shape[1] != 3:
        raise NameError('Dimension of target mesh must be 3 or 2.')

    mu[np.abs(mu) > 1e3] = 1e3
    mu[np.isnan(mu)] = 1

    return mu
