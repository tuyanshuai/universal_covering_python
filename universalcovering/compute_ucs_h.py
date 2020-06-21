from universalcovering.bd_multiplicity import *
from universalcovering.compute_bd_chain import *
from universalcovering.compute_decks_from_chain_h import *
from universalcovering.compute_bp_from_decks import *
from universalcovering.segment_pair import *
from universalcovering.compute_decks_from_bp_h import *
import numpy as np


def compute_ucs_h(face, vertex, z, hb, father):
    bd, bdn, bdc = bd_multiplicity(face, z, hb, father)
    chain, mc = compute_bd_chain(bd, bdn, bdc)
    decks = compute_decks_from_chain_h(z, chain, mc)
    bp = compute_bp_from_decks(decks, mc, z[bd[0]])
    sp, ms = segment_pair(bd, father[bd])

    nb = bp.shape[0]
    pieces = [list()] * nb
    ucs_vi = [list()] * nb
    nv = z.shape[0]

    for i in range(nb):
        z_t = z
        bp_t = bp
        ui = [list()] * (nb - 1)
        ui_vi = [list()] * (nb - 1)
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
            if sp[s].shape[0] == len(chain[k0]):
                if np.all(sp[s] == chain[k0]):
                    break

            s = s + 1
            if s >= ms.shape[0]:
                s = s - ms.shape[0]

        s = s + 1

        if s >= ms.shape[0]:
            s = s - ms.shape[0]

        vi0 = np.ones(nv, dtype=np.bool)


        #
        while True:
            if sp[s].shape[0] == len(chain[k]):
                if np.all(sp[s] == chain[k]):
                    break
            vi0[sp[ms[s]][:-1]] = False
            vi0[sp[s][:-1]] = True
            s = s + 1
            if s >= ms.shape[0]:
                s = s - ms.shape[0]

        s = ms[s]
        vi0[sp[s][:-1]] = False
        # TODO: check Sec2
        for j in range(1, nb - 1):
            decks = compute_decks_from_bp_h(bp_t, mc)
            z_t = decks[mc[k]](z_t)
            bp_t = decks[mc[k]](bp_t)
            ui[j] = z_t

            vi = np.ones(nv, dtype=np.bool)

            k = mc[k] + 1
            if k >= nb:
                k = k - nb

            while True:
                if sp[s].shape[0] == len(chain[k]):
                    if np.all(sp[s] == chain[k]):
                        break

                vi[sp[s][:-1]] = 1 - vi0[sp[s][:-1]]
                vi0[sp[s][:-1]] = True
                vi0[sp[ms[s]]][1::-1] = False

                s = s + 1
                if s >= ms.shape[0]:
                    s = s - ms.shape[0]

            s = ms[s]
            vi[chain[k]] = False

            vi0[chain[k][:-1]] = True
            vi0[chain[mc[k]]][:-1] = False
            ui_vi[j] = vi

        pieces[i] = ui
        ucs_vi[i] = ui_vi

    # Sec3
    vi = np.ones(nv, dtype=np.bool)
    vi[bd] = False

    ci = 0
    ci = ci + np.sum(vi)

    for i in range(nb):
        ui_vi = ucs_vi[i]
        for j in range(1, nb - 1):
            vi = ui_vi[j]
            ci = ci + np.sum(vi)

    ucs_z = np.zeros(ci, np.complex)
    ucs_vertex = np.zeros((ci, vertex.shape[1]))

    ci = 0
    vi = np.ones(nv, dtype=np.bool)
    vi[bd] = False
    ucs_z[ci:ci + np.sum(vi)] = z[vi]
    ucs_vertex[ci:ci + np.sum(vi), :] = vertex[vi, :]
    ci = ci + np.sum(vi)

    for i in range(nb):
        ui = pieces[i]
        ui_vi = ucs_vi[i]

        for j in range(1, nb - 1):
            uij = ui[j]
            vi = ui_vi[j]
            ucs_z[ci:ci + np.sum(vi)] = uij[vi]
            ucs_vertex[ci:ci + np.sum(vi), :] = vertex[vi, :]
            ci = ci + np.sum(vi)

    ucs = {"pieces": pieces,
           "z": ucs_z,
           "vertex": ucs_vertex,
           "bp": bp}
    return ucs
