from meshio import *
from parameterization import *
from graphics import *
import numpy as np
from topology import *

from dbgtool.dbgtool import *

face, vertex = read_off('data/eightsim.off')

# show_mesh(face, vertex)
u = hyperbolic_ricci_flow(face, vertex)
bi = 16
hb = compute_greedy_homotopy_basis(face, vertex, bi)
face_new,vertex_new,father = slice_mesh(face, vertex, hb)

