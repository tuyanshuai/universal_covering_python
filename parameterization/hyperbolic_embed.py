# Embed according to the hyperbolic metric
import numpy as np
from algebra import *
from math import *
from parameterization.circle_circle_intersection import circle_circle_intersection

# TODO : check the correctness of the function
def hyper_circle_to_circle(c, r):
    d = abs(c) * abs(c)
    f = (exp(r) - 1.0) / (exp(r) + 1.0)
    f = f * f
    a = 1 - f * d
    b = 2 * (f - 1.0) * c / a
    d = (d - f) / a
    c = -b / 2.0
    r = sqrt(abs(b) * abs(b) / 4.0 - d)
    return c, r


def hyperbolic_embed(face, u):
    nv = u.shape[0]
    z = np.zeros(nv, dtype=np.complex)
    ind = np.zeros(nv, dtype=np.bool)
    indf = np.zeros(face.shape[0], dtype=np.bool)

    def el(i, j):
        return -np.log(np.tanh(-u[i] / 2.0)) - np.log(np.tanh(-u[j] / 2.0))

    def embed_face(fi):
        si = np.sum(ind[fi])
        if si == 0 or si == 1:
            raise NameError("ERROR: wrong order of embedding")
        order = [0, 1, 2]
        if not ind[fi[0]]:
            fi = fi[[1, 2, 0]]
            order = [2, 0, 1]
        else:
            if not ind[fi[1]]:
                fi = fi[[2, 0, 1]]
                order = [1, 2, 0]
        zi = list([0, 0, 0])
        zi[0] = z[fi[0]]
        zi[1] = z[fi[1]]
        r1 = el(fi[0], fi[2])
        r2 = el(fi[1], fi[2])

        c1, r1 = hyper_circle_to_circle(zi[0], r1)
        c2, r2 = hyper_circle_to_circle(zi[1], r2)

        p = circle_circle_intersection(c1, r1, c2, r2)
        if np.isnan(p):
            raise NameError('ERROR: two circles do not intersect, invalid metric')

        zi[2] = p
        zi = np.array(zi)[order]
        ind[fi] = True

        return zi

    fr = compute_face_ring(face)
    z[face[0, 0]] = 0
    r = el(face[0, 0], face[0, 1])
    z[face[0, 1]] = (exp(r) - 1.0) / (exp(r) + 1.0)
    ind[face[0, [0, 1]]] = True
    queue = list([0])
    while len(queue) > 0:
        i = queue[0]
        if indf[i]:
            queue = queue[1::]
            continue
        z[face[i, :]] = embed_face(face[i, :])
        indf[i] = True
        f2 = fr[i, :]
        in_ = f2 > 0
        f2 = f2[in_]

        queue = queue[1::] + (f2-1).tolist()[0]


    return z
