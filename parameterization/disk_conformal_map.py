from parameterization.disk_harmonic_map import *
from algebra.compute_bc import *
from algebra.compute_bd import *
from algebra.linear_beltrami_solver import *
import numpy as np


def disk_conformal_map(face, vertex):
    uv = disk_harmonic_map(face, vertex)
    bd = compute_bd(face)
    mu = compute_bc(face, uv, vertex)
    for i in range(10):
        mu = mu / np.abs(mu) * np.mean(np.abs(mu))
        uv, _u = linear_beltrami_solver(face, uv, mu, bd, uv[bd, :])
        mu = compute_bc(face, uv, vertex)
    return uv
