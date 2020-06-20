from meshio import *
from parameterization import *
from graphics import *
import numpy as np
from topology import *
from universalcovering import *
from dbgtool.dbgtool import *

face, vertex = read_off('data/eightsim.off')

# show_mesh(face, vertex)
u = hyperbolic_ricci_flow(face, vertex)
bi = 0
hb = compute_greedy_homotopy_basis(face, vertex, bi)
face_new,vertex_new,father = slice_mesh(face, vertex, hb)
# show_mesh(face_new, vertex_new)
z = hyperbolic_embed(face_new,  u[father])
z = move_mc_to_zero(z)

ucs = compute_ucs_h(face_new,vertex_new,z,hb,father)
# show_mesh(face_new, np.concatenate((z.real[:,None],z.imag[:,None], z.real[:,None]*0),axis=1))
# plot_mesh(face_new, np.concatenate((z.real[:,None],z.imag[:,None]),axis=1))
