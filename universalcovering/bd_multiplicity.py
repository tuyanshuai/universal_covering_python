from algebra import compute_bd
import numpy as np
from dbgtool.dbgtool import *

def bd_multiplicity(face_o, z, hb, father):
    bd = compute_bd(face_o)
    vj = np.zeros(z.shape[0], dtype=np.intc)
    vc = vj.copy()

    for i in range(len(hb)):
        loop = hb[i]
        for loopj in loop:
            bj = loopj == father
            vj[bj] = vj[bj] + 1
            vc[bj] = i+1

    bdn = vj[bd]
    bdc = vc[bd]

    ns = np.max(bdn)
    i = np.argwhere(bdn == ns)[0][0]

    bd =  np.concatenate((bd[i::],  bd[0:i]))
    bdn = np.concatenate((bdn[i::], bdn[0:i]))
    bdc = np.concatenate((bdc[i::], bdc[0:i]))

    return bd, bdn, bdc
