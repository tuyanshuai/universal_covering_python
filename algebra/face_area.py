## face area 
# Compute area of all face
#
## Syntax
#   fa = face_area(face,vertex)
#
## Description
#  face  : double array, nf x 3, connectivity of mesh
#  vertex: double array, nv x 3, vertex of mesh
# 
#  fa: double array, nf x 1, area of all faces.
#  
from numpy import *

def face_area(face,vertex):
    fi = face[:,0]
    fj = face[:,1]
    fk = face[:,2]
    vij = vertex[fj,:]-vertex[fi,:]
    vjk = vertex[fk,:]-vertex[fj,:]
    vki = vertex[fi,:]-vertex[fk,:]
    a = sqrt(np.linalg.norm(vij,axis=1))
    b = sqrt(np.linalg.norm(vjk,axis=1))
    c = sqrt(np.linalg.norm(vki,axis=1))
    s = (a+b+c)/2.0
    fa = sqrt(s * (s-a)*(s-b)*(s-c))
