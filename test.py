from meshio import *
from parameterization import *
from graphics import *
import numpy as np

F, V = read_obj('data/bunny.obj')

plot_mesh(F, V)
