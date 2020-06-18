"""[pd2,h] = discrete_optimal_transport(disk,face,uv,sigma,area);
"""
from matplotlib import path
from numpy import *
from disk_omt.power_diagram import *
import warnings
from disk_omt.intersectRayPolygon import *
from scipy.linalg import norm
from scipy.sparse.linalg import spsolve

def isinpolygon(polygon, xy):
    p = path.Path(polygon)
    flag = p.contains_points(xy)
    return flag


def polyarea(xy):
    n = xy.shape[0]
    # Initialze area
    area = 0.0
    # Calculate value of shoelace formula
    j = n - 1
    for i in range(0, n):
        area += (xy[j,0] + xy[i,0]) * (xy[j,1] - xy[i,1])
        j = i  # j is previous vertex to i
    # Return absolute value
    return abs(area / 2.0)

#
# def polybool(cp, ci):
#     flag = isinpolygon(ci, cp)
#     kpt = argwhere(flag).flatten().astype(int)
#     for i in range(ci.shape[0]):
#         len = sqrt(ci[i, 0] * ci[i, 0] + ci[i, 1] * ci[i, 1])
#         if len > 0.9:
#             ci[i,:] = intersectRayPolygon(array([0.0, 0.0]), ci[i,:], concatenate((cp, cp[0].reshape(1,2))))
#
#
#     newpt = concatenate((cp[kpt,:], ci), axis=0)
#     if (newpt.shape[0] <4):
#         return newpt
#
#     K = ConvexHull(newpt, qhull_options='Qt').vertices
#     K = append(K,K[0])
#     return newpt[K,:]

def polybool(cp, ci):
    flag = isinpolygon(ci, cp)
    kpt = argwhere(flag).flatten().astype(int)
    for i in range(ci.shape[0]):
        len = sqrt(ci[i, 0] * ci[i, 0] + ci[i, 1] * ci[i, 1])
        if len > 0.99:
            ci[i,:] = ci[i,:]/len

    return ci










def calculate_gradient(cp, pd, sigma, more_accurate = False):
    nc = len(pd["cell"])
    D = zeros((nc, ))

    isin = ones((nc, )).astype(bool)
    isin2 = isinpolygon(cp, pd["dpe"])

    for i in range(nc):
        ci = pd["cell"][i]
        if not all(isin2[ci]):
            isin[i] = False
    cp = delete(cp,-1, 0)

    for i in range(nc):
        ci = pd["dpe"][pd["cell"][i],:]

        if isin[i]:
            ci = ci[0:-1]
            mui = (2 * mean(sigma(ci)) + sigma(mean(ci, axis=0).reshape(-1,2))) / 3

            if ci.size ==0:
                D[i] = eps;
                warning('occured empty D');
            else:
                D[i] = polyarea(ci) * mui;

        else:

           # We simplify the problem, as the ci shall inside of radio 1 circle
            xy = polybool(cp, ci)


            mui = (2 * mean(sigma(xy)) + sigma(mean(xy, axis =0 ).reshape(-1,2))) / 3;

            if xy.size ==0 :
                D[i] = eps;
                warnings.warn('occured empty D')
            else:
                D[i] = polyarea(xy) * mui

    return D

def calculate_hessian(cp,pd,sigma):
    nc = len(pd["cell"])
    ne = 0
    for nei in pd["cell"]:
        ne += len(nei)
    ne -= nc
    I = zeros((ne,)).astype(int)
    J = zeros((ne,)).astype(int)
    V = zeros((ne,)).astype(int)
    k = 0
    for i in range(nc):
        ci = pd["cell"][i]
        for j in range(len(ci)-1):
            I[k] = ci[j]
            J[k] = ci[j + 1]
            V[k] = i
            k = k + 1
    C = sparse.coo_matrix((V.flatten()+1, (I.flatten(), J.flatten())))
    C = C.tolil()
    IJ =argwhere(C)
    I = IJ[:, 0]
    J = IJ[:, 1]

    I2 = zeros( (ne, )).astype(int)
    J2 = zeros( (ne, )).astype(int)
    V2 = zeros( (ne, ))
    isin = isinpolygon(cp, pd["dpe"])

    p = sigma(pd["dpe"])
    if len(p) == 1:
        p = ones((pd["dpe"].shape[0]))*p

    k = 0
    for i in range(I.shape[0]):
        I2[k] = C[J[i], I[i]]-1
        J2[k] = C[I[i], J[i]]-1
        # compute         edge        length in convex        polygon
        p1 = pd["dpe"][I[i],:]
        p2 = pd["dpe"][J[i],:]
        in2 = isin[array([I[i], J[i]])]

        if sum(in2) ==2:
            lij = norm(p1 - p2) * (p[I[i]] + p[J[i]]) / 2.0
        if sum(in2) == 1:

            pi1 = intersectRayPolygon(p1, p2 - p1, cp)
            pi2 = intersectRayPolygon(p2, p1 - p2, cp)

            if pi1 is None:
                pi = pi2
            else:
                if pi2 is None:
                    pi = pi1
                else:
                    if norm( (p1+p2)/2.0 - pi1) <   norm( (p1+p2)/2.0 - pi2):
                        pi =pi1
                    else:
                        pi = pi2

            if in2[0]:
                lij = norm(pi - p1) * (sigma(pi) + p[I[i]] ) / 2.0
            else:
                lij = norm(pi - p2) * (sigma(pi) + p[J[i]] ) / 2.0
        if sum(in2) == 0:
            lij =0

        V2[k] = -lij / norm(pd["uv"][I2[k],:]-pd["uv"][J2[k],:])
        k = k + 1

    H = sparse.coo_matrix((V2.flatten(), (I2.flatten(), J2.flatten())))
    H = H.tolil()
    Hs = -sum(H, axis=1)
    for i in range(H.shape[0]):
        H[i,i] += Hs[i,0]
    H = (H + H.transpose())/2.0

    return H




def discrete_optimal_transport(cp, face, uv, sigma, delta, h = None, max_iter =100, eps = 1e-6 ):
    npuv = uv.shape[0]
    if h is None:
        h = zeros((npuv,1))

    pd, _ = power_diagram(face, uv, h)
    k = 1
    while k < max_iter:
        G = calculate_gradient(cp, pd, sigma, 1)
        G = G / sum(G) * sum(delta)
        D = G - delta
        H = calculate_hessian(cp, pd, sigma)

        H[0,0] = H[0,0] +1.0
        dh = spsolve(H, D)
        dh = dh.reshape((-1,1))
        if not all(isfinite(dh)):
            raise ValueError("""ERROR: |dh| goes infinite, most probably due to convexhull 
                   failing, which is due to some cell(s) disappear. Real reason 
                   is mesh quality/measure is too bad""")
        dh = dh - mean(dh)
        dh = dh - mean(dh)

        maxdh = max(abs(dh));
        print('#%02d: max|dh| = %.10f' %( k, maxdh))

        if maxdh < eps:
            break

        pd, h = power_diagram(face, uv, h, dh)
        k = k + 1


    return pd, h, maxdh