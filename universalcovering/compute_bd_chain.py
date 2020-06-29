import  numpy as np
from dbgtool.dbgtool import *
def compute_bd_chain(bd,bdn,bdc):

    nb = bd.shape[0]  # number of boundaries
    ns = np.max(bdn)
    chain = [list()]*ns

    s = 0
    k = 0
    mc = np.zeros(ns, dtype=np.intc)
    c = np.zeros(ns, dtype=np.intc)

    while k < nb:
        i = k
        while i<= nb and bdn[i] != 1:
            i = i +1

        if i >= nb:
            break

        j = i

        while j < nb and bdn[j] ==1:
            j = j+1

        jn = j
        if jn >= nb:
            jn = jn -nb

        c[s] = bdc[i]
        chain[s] = list(bd[i-1:j]) + list([bd[jn]])
        s = s+1

        if s >= ns:
            break
        k = j

    for s in range(ns):
        i = np.argwhere(c==c[s]).flatten()[0:2]
        mc[i] = np.array([i[1],i[0]])

    return chain, mc
