# pc_spherical_conformal_map

# Conformally map a genus-0 point cloud to the unit sphere
#
# Usage:
# map = pc_spherical_conformal_map(vertex)
#
# Input:
# vertex: nv x 3 vertex coordinates of the input point cloud
#
# Output:
# map: nv x 3 vertex coordinates of the spherical conformal parameterization
#
# If you use this code in your own work, please cite the following paper:
# [1] G. P.-T. Choi, K. T. Ho and L. M. Lui,
#     "Spherical Conformal Parameterization of Genus-0 Point Clouds for Meshing."
#     SIAM Journal on Imaging Sciences, 9(4), pp. 1582-1618, 2016.
from scipy.spatial import ConvexHull
import numpy as np
from scipy.spatial import distance
def knnsearch(vertex,k):
    D = distance.squareform(distance.pdist(vertex))
    closest = np.argsort(D, axis=1)
    knnInd = closest[:, 1:k + 1]
    return knnInd

def calc_pc_laplacian(vertex, k):

    numofv = vertex.shape[0]
    rings = [list()] * numofv

    if vertex.shape[1]  == 2:
        vertex = np.concatenate( (vertex, np.zeros(numofv)), axis =1)

    knnInd = knnsearch(vertex, k)
    coeff = np.ones(6)
    powerc = np.array([[0,0],[1,0],[0,1],[2,0],[1,1],[0,2]])

    L, rings = None, None
    return L, rings
def pc_spherical_conformal_map(vertex):
    L, rings = calc_pc_laplacian(vertex, 25)

    return None