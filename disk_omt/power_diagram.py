""""
power diagram
"""
import numpy as np
from scipy.linalg import norm

from scipy.spatial import ConvexHull
from disk_omt.calculate_face_normal import *
from disk_omt.intersectRayPolygon import *
from algebra import *


def face_dual_uv(p):
    a = p[0, 1] * (p[1, 2] - p[2, 2]) + p[1, 1] * (p[2, 2] - p[0, 2]) + p[2, 1] * (p[0, 2] - p[1, 2])
    b = p[0, 2] * (p[1, 0] - p[2, 0]) + p[1, 2] * (p[2, 0] - p[0, 0]) + p[2, 2] * (p[0, 0] - p[1, 0])
    c = p[0, 0] * (p[1, 1] - p[2, 1]) + p[1, 0] * (p[2, 1] - p[0, 1]) + p[2, 0] * (p[0, 1] - p[1, 1])
    dp = [-a / c / 2, -b / c / 2]
    return dp





def power_diagram(face, uv, h=None, dh=None):
    if h is None:
        h = np.zeros((uv.shape[0], 1))

    if dh is None:
        dh = h * 0

    nf = face.shape[0]
    c = 1

    while True:
        h = h - c * dh
        pl = np.concatenate((uv, np.reshape(np.square(norm(uv, axis=1)), (-1, 1)) - h), axis=1)
        hull = ConvexHull(pl, qhull_options='Qt')
        face = hull.simplices
        # fix ups for the convex hull, as the orientation may inverse
        fn_from_hull = hull.equations[:,2]
        fn = calculate_face_normal(face, pl)
        for i in range(face.shape[0]):
            if fn[i,2] * fn_from_hull[i] < 0 :  # orientation difff
                face[i,:] = face[i,[0, 2, 1]]



        for i in range(face.shape[0]):
            mif = np.argmin(face[i,:])
            face[i, :] = face[i, np.mod(np.arange(mif,mif+3),3)]
        face = face[np.argsort(face[:, 0] * np.max(face) + face[:, 1]), :]
        fn = calculate_face_normal(face, pl)
        ind = fn[:, 2] < 0

        if np.sum(ind) < nf:
            h = h + c * dh
            c = c / 2
        else:
            break

        if np.max(abs(dh)) == 0:
            break

    fn = calculate_face_normal(face, pl)
    ind = fn[:, 2] < 0
    face = face[ind, :]
    pd = dict()
    pd['face'] = face
    vr = compute_vertex_ring(face, uv, ordered=True)
    pd['uv'] = uv
    pd['dp'] = np.zeros((face.shape[0], 2))
    pd['cell'] = [[] for i in range(pl.shape[0])]

    for i in range(face.shape[0]):
        dp = face_dual_uv(pl[face[i,:],:])
        pd['dp'][i,:] = dp

    K =  ConvexHull(uv, qhull_options='Qt').vertices
    ks = np.argmin(K)
    K = np.concatenate((K[ks::],  K[0:ks]), axis=0)
    K = np.append(K,K[0])
    vb = np.zeros((K.shape[0] - 1, 2))
    mindp = np.min(pd["dp"], axis=0) - 1
    maxdp = np.max(pd["dp"], axis=0) + 1
    minx = mindp[0]
    miny = mindp[1]
    maxx = maxdp[0]
    maxy = maxdp[1]
    box = np.array([minx, miny, maxx, miny, maxx, maxy, minx, maxy, minx, miny]).reshape((-1,2))

    for i in range(K.shape[0]- 1):
        i1 = K[i]
        i2 = K[i + 1]
        vec = uv[i2,:] - uv[i1,:]
        vec = np.array([vec[1], -vec[0]])
        mid = (uv[i2,:] + uv[i1,:]) / 2.0
        intersect = intersectRayPolygon(mid, vec, box)
        vb[i,:] = intersect

    pd["dpe"] = np.concatenate((pd["dp"], vb), axis=0)

    vvif, _, _= compute_connectivity(face)

    for i in range(uv.shape[0]):
        vri = vr[i]
        pb = np.argwhere(K==i)
        if pb.size > 0 :
            pb = pb[0][0]
            fr = np.zeros((len(vri) + 1,)).astype(int)
            fr[-1] = face.shape[0] + pb
            if pb == 0:
                fr[0] = face.shape[0] + K.shape[0]-2
            else:
                fr[0] = face.shape[0] + pb - 1
            for j in range(len(vri) - 1):
                fr[j+1] = vvif[i, vri[j]]
        else:
            fr = np.zeros((len(vri),)).astype(int)
            for j in range(len(vri)):
                fr[j] = vvif[i, vri[j]]
        pd["cell"][i] = np.flip(fr)

    return pd, h
