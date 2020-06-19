"""

function vr = compute_vertex_ring(face,vertex,vc,ordered)

"""
import numpy as np
from algebra import *
from algebra import compute_bd

from dbgtool.dbgtool import *
#TODO: modify to support order
def compute_vertex_ring(face, vertex, vc = np.array([]), ordered=None):
    nv = np.max(np.max(face)) + 1
    if vc.size ==0:
        vc = np.array(range(nv))
    if ordered == None:
        ordered = False

    vr = [[] for i in range(vc.size)]
    bd = compute_bd(face)
    isbd = np.zeros((nv, )).astype(int)
    isbd[bd] = 1

    am, _ = compute_adjacency_matrix(face)
    IJ = np.argwhere(am[:, vc])
    I = IJ[:, 1]
    J = IJ[:, 0]

    for i in range(J.shape[0]):
        vr[J[i]].append(I[i])

    if ordered:
        for i in range(vc.size):

            vai = np.append(np.append(np.argwhere((face[:, 0] == i)), np.argwhere((face[:, 1] == i))),
                            np.argwhere((face[:, 2] == i)))
            fai = face[vai]
            if fai.size ==0 :
                vr[i] = list()
                continue

            _, amd = compute_adjacency_matrix(fai)

            vri = vr[i]
            vrinew = list()
            ai  = vri[0]
            vrinew.append(ai)
            for trytime in range(len(vri)):
                q = vrinew[-1]
                for bq in vri:
                    if amd[q, bq]:
                        vrinew.append(bq)
                        break

                if len(vrinew) >= len(vri):
                    break

                q = vrinew[0]
                for bq in vri:
                    if amd[bq, q]:
                        vrinew.insert(0,bq)
                        break

                if len(vrinew) >= len(vri):
                    break
            # Sort the vri


            if not isbd[i]:
                ks = np.argmin(vrinew)
                vrinew = np.concatenate((vrinew[ks::], vrinew[0:ks]), axis=0)
                vrinew = np.append(vrinew, vrinew[0])
            vr[i] = vrinew

            # vai = np.append(np.append(np.argwhere((face[:, 0] == i)), np.argwhere((face[:, 1] == i))),
            #                 np.argwhere((face[:, 2] == i)))
            # fai = face[vai]
            # edge_list = np.append(np.append(fai[:, [0, 1]], fai[:, [1, 2]], axis=0), fai[:, [2, 0]], axis=0)
            #
            # edge_list = np.delete(edge_list, np.argwhere(edge_list[:, 0] == i), axis=0)
            # edge_list = np.delete(edge_list, np.argwhere(edge_list[:, 1] == i), axis=0)
            #
            # # edge_list to list
            # vri = list()
            # if edge_list.shape[0] > 0:  # if empty of face list, none
            #     if not isbd[i, 0]:  # if not boundary index, find a loop
            #         vri.append(edge_list[0, 0])
            #         while True:
            #             new_node = edge_list[np.where(edge_list[:, 0] == vri[-1]), 1]
            #             if new_node.size == 0:
            #                 break
            #             new_node = new_node[0][0]
            #             vri.append(new_node)
            #             if new_node == vri[0]:
            #                 break
            #
            #     else:
            #         vri.append(edge_list[0, 0])
            #         while True:
            #             new_node = edge_list[np.where(edge_list[:, 0] == vri[-1]), 1]
            #             if new_node.size == 0:
            #                 break
            #             vri.append(new_node[0][0])
            #         # back add
            #         while True:
            #             new_node = edge_list[np.where(edge_list[:, 1] == vri[0]), 0]
            #             if new_node.size == 0:
            #                 break
            #             vri.insert(0, new_node[0][0])
            # vr[i] = vri

    return vr
