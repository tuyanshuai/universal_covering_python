import numpy as np
from dbgtool.dbgtool import *

def move_mc_to_zero(z):
    # move the mass center of mesh to zero by Mobius transform on Poincare
    # disk, regard points as complex number
    z_t = np.array(z)
    zc = 1
    while np.absolute(zc) > 1e-6:
        zc = np.sum(z_t)/z_t.shape[0]
        for i in range(z_t.shape[0]):
            z_t[i] = (z_t[i]-zc) / (1-np.conj(zc)*z_t[i])

    return z_t
