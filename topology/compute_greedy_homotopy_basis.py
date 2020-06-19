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

def gen_path(pred, i):
    path = list([i])
    while True:
        if pred[path[0]]<0:
            break
        path.insert(0, pred[path[0]])

    return path

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

    # find shortest path to graph
    dist, pred = sp.csgraph.shortest_path(csgraph=G, directed=False, indices=bi, return_predecessors=True)
    pathes = [list()] * nv
    for di in range(nv): # destination node
        pathes[di] = gen_path(pred, di)

    amf, _ = compute_dual_graph(face)

    I = np.array(range(nv))
    I = np.delete(I,bi)
    J = pred[I]

    I2 = np.array(amd[I,J]).flatten()
    J2 = np.array(amd[J,I]).flatten()
    ind = np.argwhere(np.logical_or(I2 == 0 , J2 == 0))

    I2 = np.delete(I2, ind, 0)
    J2 = np.delete(J2, ind, 0)

    amf[I2-1, J2-1] = 0
    amf[J2-1, I2-1] = 0


    IJ = np.argwhere(amf)
    I = IJ[:, 1]
    J = IJ[:, 0]

    ind = np.argwhere(np.logical_or(eif[:, 0] == -1, eif[:, 1] == -1))
    eif = np.delete(eif, ind,0)
    edge = np.delete(edge, ind,0)


    ti = np.concatenate((eif[:, 0], eif[:, 1]))
    tj = np.concatenate((eif[:, 1], eif[:, 0]))
    tv = np.concatenate((edge[:, 0], edge[:, 1]))
    F2E = sp.csr_matrix((tv, (ti, tj)), shape=(nf, nf))

    ei = np.concatenate((F2E[I, J], F2E[J,I]), axis=0).transpose()

    dvi = vertex[ei[:, 0].flatten(),:]-vertex[ei[:, 1].flatten(),:]
    V = -(dist[ei[:, 0]].flatten()+dist[ei[:, 1]].flatten()+ np.linalg.norm(dvi[0], axis=1))
    amf_w = sp.csr_matrix((V.flatten(), (I, J)), shape=(nf, nf))


    tree = sp.csgraph.minimum_spanning_tree(csgraph=amf_w, overwrite=False)

    # tree = tree + tree.transpose()
    G2 = G.copy()
    I = np.array(range(nv))
    I = np.delete(I, bi, 0)
    J = pred[I]
    G2[I, J] = 0
    G2[J, I] = 0

    idij = np.argwhere(tree )
    I = idij[:, 0]
    J = idij[:, 1]

    ei = np.concatenate((F2E[I, J], F2E[J,I]), axis=0).transpose()
    G2[ei[:, 1], ei[:, 0]] = 0
    G2[ei[:, 0], ei[:, 1]] = 0

    idij = np.argwhere(G2)
    I = idij[:, 0]
    J = idij[:, 1]

    id2del = np.argwhere(I<J)
    I = np.delete(I, id2del)
    J = np.delete(J, id2del)

    argj = np.argsort(J)
    J = J[argj]
    I = I[argj]


    hb = [list([])]*I.size
    for i in range(I.size):
        pi = pathes[I[i]]
        pj = pathes[J[i]]
        hb[i] = pi+pj[::-1]

    if I.size ==0:
       hb = np.array([])
    return hb