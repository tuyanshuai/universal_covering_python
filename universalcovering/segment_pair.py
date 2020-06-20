# find the corresponding segment pairs in boundary edges, which have same
# father seqence but opposite direction
import numpy as np
from dbgtool.dbgtool import *


def segment_pair(bd, father):
    father2 = np.append(father, father[0])
    father2_r = np.flipud(father2)
    n = father2.shape[0]
    i = 1
    j = i
    si = 1
    sp = list()
    sd = np.zeros(n)

    while np.argwhere(sd).shape[0] < n and j < n:
        l_ = 1
        j2 = 0
        while l_ + j - i < n and j < n:
            if np.all(father2[i - 1:j + 1] == father2_r[l_ - 1: l_ + j - i + 1]):
                j = j + 1
                j2 = l_
            else:
                l_ = l_ + 1

        if j2 is not 0:

            spi = np.row_stack((np.array(range(i, j + 1)),
                                np.array(range(n + 1 - j2 - j + i, n - j2 + 2))))
            spi = spi.transpose()

            b = spi == n
            spi[b] = 1

            sp.append(bd[spi - 1])
            sd[i - 1:j] = si
            si = si + 1
            i = j
            j = i
        else:
            i = i + 1
            j = i

    ms = np.zeros(len(sp), dtype=np.intc)
    for i in range(ms.shape[0]):
        si = sp[i]
        for j in range(ms.shape[0]):
            if j == i:
                continue
            sj = sp[j]

            if si.shape[0] == sj.shape[0]:
                if np.all(si == np.fliplr(sj)):
                    ms[i] = j

    for i in range(len(sp)):
        sp[i] = sp[i][:, 0]

    return sp, ms
