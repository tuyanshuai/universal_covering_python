# Embed according to the hyperbolic metric
from algebra import *
from math import *
from universalcovering.circle_circle_intersection import circle_circle_intersection


def euclidean_embed(face, u):
    nv = u.shape[0]
    z = np.zeros(nv, np.complex)
    ind = np.zeros(nv, dtype=np.bool)
    indf = np.zeros(face.shape[0], dtype=np.bool)

    def el(i, j):
        return np.exp(u[i]) + np.exp(u[j])

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
        zi = np.zeros(3, np.complex)
        zi[0] = z[fi[0]]
        zi[1] = z[fi[1]]
        r1 = el(fi[0], fi[2])
        r2 = el(fi[1], fi[2])

        p = circle_circle_intersection(zi[0], r1, zi[1], r2)
        if np.isnan(p):
            raise NameError('ERROR: two circles do not intersect, invalid metric')

        zi[2] = p
        zi = np.array(zi)[order]
        ind[fi] = True

        return zi

    fr = compute_face_ring(face)
    z[face[0, 0]] = 0
    z[face[0, 1]] = el(face[0, 0], face[0, 1])
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

        queue = queue[1::] + (f2 - 1).tolist()[0]
    return z
