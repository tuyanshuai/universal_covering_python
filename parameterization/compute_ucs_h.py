from algebra.bd_multiplicity import *
from algebra.compute_bd_chain import *
from algebra.compute_decks_from_chain_h import *
from algebra.compute_bp_from_decks import *
from algebra.segment_pair import *
import numpy as np

def  compute_ucs_h(face,vertex,z,hb,father):

    bd,bdn,bdc = bd_multiplicity(face,z,hb,father)
    chain, mc = compute_bd_chain(bd, bdn, bdc)
    decks = compute_decks_from_chain_h(z, chain, mc)
    bp = compute_bp_from_decks(decks, mc, z(bd[0]))
    sp, ms= segment_pair(bd, father(bd))

    nb = bp.shape[0]
    pieces = [list()]*nb
    ucs_vi = [list()]*nb
    nv = z.shape[0]

    for i in range(nb):
        z_t = z
        bp_t = bp
        ui =  [list()]*(nb-1)
        ui_vi = [list()]*(nb-1)
        ui[0] = z_t
        vi = np.ones(nv, dtype=np.bool)
        vi[bd] = False

        ui_vi[0] = vi
        k = i
        k0 = k - 1
        if k0 < 0:
            k0 = k0 + nb
        s = 0
        while True:
            if sp[s].shape[0] == chain[k0].shape[0]:
                if sp[s] == chain[k0]:
                    break

            s = s + 1
            if s > ms.shape[0]:
                s = s - ms.shape[0]

        s = s + 1

        if s > ms.shape[0]:
            s = s - ms.shape[0]
        # TODO Line 46
        while True:
            0

    return None