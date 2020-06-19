## slice_mesh 
# Slice mesh open along a collection of edges ee, which usually comes from 
# cut_graph(directly), or compute_greedy_homotopy_basis and
# compute_homology_basis (need to form edges from closed loops in basis).
# ee can form a single closed loops or multiple closed loops. 
# 
## Syntax
#   [face_new,vertex_new,father] = slice_mesh(face,vertex,ee)
#
## Description
#  face  : double array, nf x 3, connectivity of mesh
#  vertex: double array, nv x 3, vertex of mesh
#  ee    : double array, n x 2, a collection of edges, each row is an edge on 
#          mesh, may not be in consecutive order. 
# 
#  face_new  : double array, nf x 3, connectivity of new mesh after slice
#  vertex_new: double array, nv' x 3, vertex of new mesh, vertex number is
#              more than original mesh, since slice mesh will separate each
#              vertex on ee to two vertices or more.
#  father    : double array, nv' x 1, father indicates the vertex on original
#              mesh that new vertex comes from.
#
def slice_mesh(face,vertex,ee):


    return None