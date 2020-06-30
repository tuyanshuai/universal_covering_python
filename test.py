from meshio import *
from parameterization import *
from graphics import *
import numpy as np
from topology import *
from universalcovering import *
from dbgtool.dbgtool import *



F, V = read_off('data/maxplanck.nf25k.off')
tic()
uvw = spherical_conformal_map(F, V)
toc()

# plot_mesh(F, uvw)

show_mesh(F, uvw, V)


face, vertex = read_off('data/eightsim.off')

# show_mesh(face, vertex)
u = hyperbolic_ricci_flow(face, vertex)
bi = 1
hb = compute_greedy_homotopy_basis(face, vertex, bi)
face_new, vertex_new, father = slice_mesh(face, vertex, hb)
# show_mesh(face_new, vertex_new)
z = hyperbolic_embed(face_new, u[father])
z = move_mc_to_zero(z)

ucs = compute_ucs_h(face_new, vertex_new, z, hb, father)

# show_mesh(face_new, np.concatenate((z.real[:,None],z.imag[:,None], z.real[:,None]*0),axis=1))
# plot_mesh_Circle(face_new, np.concatenate((z.real[:,None],z.imag[:,None]),axis=1))

plot_ucs(face_new, ucs["pieces"])
#
#
# face, vertex = read_off('data/torus.off')
#
# u = euclidean_ricci_flow(face, vertex)
# bi = 200
# hb = compute_greedy_homotopy_basis(face, vertex, bi)
# face_new, vertex_new, father = slice_mesh(face, vertex, hb)
# z = euclidean_embed(face_new, u[father])
# ucs = compute_ucs(face_new, vertex_new, z, hb, father)
# plot_ucs(face_new, ucs["pieces"])
