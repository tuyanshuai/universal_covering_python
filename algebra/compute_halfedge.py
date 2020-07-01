## compute halfedge 
# Halfedge is simply directed edge, each face has three halfedges.
# This function will return all nf x 3 halfedges, as well as a nf x 3
# vector indicate which face the halfedge belongs to.
#
## Syntax
#   [he,heif] = compute_halfedge(face)
#
## Description
#  face: double array, nf x 3, connectivity of mesh
# 
#  he  : double array, (nf x 3) x 2, each row is a halfedge
#  heif: double array, (nf x 3) x 1, face id in which the halfedge lies in

import numpy as np

# TODO check
def compute_halfedge(face):
    nf = face.shape[0]
    he1 = np.concatenate((face[:, 0], face[:, 1]), axis=1)
    he2 = np.concatenate((face[:, 1], face[:, 2]), axis=1)
    he3 = np.concatenate((face[:, 2], face[:, 0]), axis=1)

    he = np.concatenate((he1, he2, he3))

    heif = np.reshape(np.repeat( np.array(range(nf)), (3,1) ), (nf*3))

    return he, heif
