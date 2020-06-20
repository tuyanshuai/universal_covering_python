## compute face ring 
# Face ring of each face. For interior face, there are three faces
# attached: each halfedge (reversed direction) attaches a neighbor face.
# 
# For boundary face (with edge on the boundary), there is no face
# attached with the boundary edge, in such case, -1 is used.
# 
## Syntax
#   fr = compute_face_ring(face)
#
## Description
#  face: double array, nf x 3, connectivity of mesh
# 
#  fr: double array, nf x 3, face ring, each row is three faces attached
#      with three edges of the face, -1 indicates boundary edge
from algebra import compute_adjacency_matrix
import numpy as np
from dbgtool.dbgtool import *

def compute_face_ring(face):
    _, amd = compute_adjacency_matrix(face)
    fr = amd[face[:, [2, 0, 1]], face[:, [1, 2, 0]]]
    fr = fr.todense()
    IJ = np.argwhere(fr == 0)
    I = IJ[:,0]
    J = IJ[:,1]


    fr[I,J] = -1
    return fr
