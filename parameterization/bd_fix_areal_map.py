"""
Boundary fixed planar area_preserving mapping

bd_fix_areal_map(face, vertex, bd, uv) for area preserving mapping
bd is the boundary indices
uv is the boundary position

"""


from algebra.face_area import *
from algebra.compute_edge import *
from scipy.sparse import linalg, csc_matrix


def dot(A, bd, axis=0):
    return np.sum(A.conj() * bd, axis=axis)


# # Laplace with area on
def simLap(face, vertex, f):
    edge, eif = compute_edge(face)
    ne = edge.shape[0]
    ew = np.zeros(ne)
    ind = eif[:, 0] > 0
    ev1 = np.sum(face[eif[ind, 0], :], axis=1) - np.sum(edge[ind, :], axis=1)

    fi = f[edge[ind, 0], :]
    fj = f[edge[ind, 1], :]
    fk = f[ev1, :]

    vi = vertex[edge[ind, 0], :]
    vj = vertex[edge[ind, 1], :]
    vk = vertex[ev1, :]

    ew[ind] = ew[ind] + dot(fi - fk, fj - fk, axis=1) / (2.0 * triangle_area(vi, vj, vk))
    ind = eif[:, 1] > 0
    ev2 = np.sum(face[eif[ind, 1], :], axis=1) - np.sum(edge[ind, :], axis=1)

    fi = f[edge[ind, 0], :]
    fj = f[edge[ind, 1], :]
    fk = f[ev2, :]

    vi = vertex[edge[ind, 0], :]
    vj = vertex[edge[ind, 1], :]
    vk = vertex[ev2, :]

    ew[ind] = ew[ind] + dot(fi - fk, fj - fk, axis=1) / (2.0 * triangle_area(vi, vj, vk))

    I = np.concatenate((edge[:, 0], edge[:, 1]))
    J = np.concatenate((edge[:, 1], edge[:, 0]))
    V = np.concatenate((ew / 2.0, ew / 2.0))
    A = csc_matrix((V, (I, J)))
    sA = np.sum(A, axis=1)
    A = A - np.diag(sA)
    return A


def bd_fix_areal_map(face, vertex, bd, uv):
    N = np.array(range(uv.shape[0]))
    I = np.setdiff1d(N, bd)

    uv_new = np.copy(uv)
    L = simLap(face, vertex, uv)

    for k in range(20):

        f1 = np.copy(uv_new)
        fbd = linalg.spsolve(-L[bd, bd], L[bd, I].dot(uv_new[I, :]))
        fbd = fbd - np.mean(fbd, axis=1)

        dl = np.sqrt(np.linalg.norm(fbd, axis=1))
        uv_new[bd, 0] = fbd[:, 0] / dl
        uv_new[bd, 1] = fbd[:, 1] / dl

        fin = linalg.spsolve(-L[I, I], L[I, bd].dot(uv_new[bd, :]))
        uv_new[I, :] = fin
        L = simLap(face, vertex, uv_new)

        max_df = np.max(np.linalg.norm(uv_new - f1, axis=1))
        print('max_df = %f\n' % max_df)
        if max_df < 1e-6:
            break
    return uv_new
