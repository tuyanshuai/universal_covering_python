## compute dual graph 
# Dual graph of a triangle mesh, regarded as graph. 
# Each face in original mesh corresponds to a vertex in dual graph, vertex
# position be the centroid of the original face.
#
## Syntax
#   [amf] = compute_dual_graph(face);
#   [amf,dual_vertex] = compute_dual_graph(face,vertex);
#
## Description
#  face  : double array, nf x 3, connectivity of mesh
#  vertex: double array, nv x 3, vertex of mesh
# 
#  amf: sparse matrix, nf x nf, connectivity of dual graph
#  dual_vertex: nf x 3, dual vertex in dual graph, if vertex is not
#               supplied, will return []
from algebra import *
import scipy.sparse as sp
import numpy as np
from dbgtool.dbgtool import  *
def compute_dual_graph(face, vertex = np.array([])):
    edge, eif = compute_edge(face)
    nf = face.shape[0]
    if vertex.size ==0:
        dual_vertex = list([])
    else:
        dual_vertex = (vertex[face[:, 0], :] + vertex[face[:, 1], :] + vertex[face[:, 2], :]) / 3.0
    ind = np.logical_and(eif[:, 0] >= 0,  eif[:, 1] >=0 )
    eif2 = eif[ind,:]

    I = eif2[:, 0]
    J = eif2[:, 1]
    V = np.ones((eif2.shape[0], 1),dtype=np.intc).flatten()
    amf = sp.csr_matrix((V, (I, J)), shape=(nf, nf))
    amf = amf + amf.transpose()

    return amf, dual_vertex