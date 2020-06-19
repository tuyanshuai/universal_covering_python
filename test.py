from meshio import *
from parameterization import *
from graphics import *
import numpy as np
from topology import *


from dbgtool.dbgtool import *

F, V = read_off('data/eightsim.off')

# show_mesh(F, V)
# compute_dual_graph(F,V)
u = hyperbolic_ricci_flow(F, V)
bi = 16
hb = compute_greedy_homotopy_basis(F, V, bi)
print(hb)



