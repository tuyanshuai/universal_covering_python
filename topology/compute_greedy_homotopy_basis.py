#   hb = compute_greedy_homotopy_basis(face,vertex,bi)
##Description
#  face  : double array, nf x 3, connectivity of mesh
#  vertex: double array, nv x 3, vertex of mesh
#  bi    : integer scaler, base point of homotopy group
#
#  hb: cell array, n x 1, a basis of homotopy group, each cell is a closed 
#      loop based at bi. Return empty for genus zero surface.

# Two graph algorithms are needed: minimum spanning tree and shortest
# path, provided in Matlab's bioinformatics toolbox. We supply 
# alternatives for these two functions, implemented purely in Matlab,
# which are a little slower and will be invoked when Maltab's built-in
# functions are not available.
#  
# # Erickson, Jeff, and Kim Whittlesey. "Greedy optimal homotopy and 
#   homology generators." Proceedings of the sixteenth annual ACM-SIAM 
#   symposium on Discrete algorithms. Society for Industrial and Applied 
#   Mathematics, 2005.
from algebra import *
import scipy.sparse as sp
from  dbgtool.dbgtool import *
def compute_greedy_homotopy_basis(face,vertex,bi):
    nf = face.shape[0]
    nv = vertex.shape[0]
    edge, eif = compute_edge(face)
    am, amd = compute_adjacency_matrix(face)

    # G is adjacency matrix for weighted graph
    idij = np.argwhere(am>0)
    I = idij[:,0]
    J = idij[:,1]

    el = np.linalg.norm(vertex[I,:]-vertex[J,:], axis=1)

    G = sp.csr_matrix((el, (I, J)), shape=(nv, nv))
    G = G + G.transpose()

    # TODO : find shortest path to graph

    # [dist,path,pred] = graphshortestpath(G,bi,'METHOD','Dijkstra');
    # % we always use column array
    # dist = dist(:);
    # pred = pred(:);


    # return hb

    print(G)
    return None