from meshio import *
from parameterization import *
from graphics import *
import numpy as np

F, V = read_off('data/eightsim.off')

plot_mesh(F, V)
