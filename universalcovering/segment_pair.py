# find the corresponding segment pairs in boundary edges, which have same
# father seqence but opposite direction
import numpy as np
def segment_pair(bd,father):
    father2 = np.concatenate((father, np.array(list(father[0]))))
    father2_r = np.flipud(father2)
    n = father2.shape[0]
    i = 0
    j = i
    si =0
    sp = list()
    sd = np.zeros(n)

    while np.argwhere(sd).shape[0] < n and j<n:
        l_ = 1
        j2 = 0
        while l_+j-i < n and j<n:
            if father2[i:j+1] == father2_r[0:1+j-i]:
                j = j+1
                j2 = l_
            else:
                l_ = l_+1

        if j2 is not 0:
            spi = np.concatenate( (np.array(range(i,j+1)), np.array(range(n+1-j2-j+i, n+1-j2))))
            b = spi ==n
            spi[b]  = 1
            sp[si] =bd[spi]
            sd[i:j+1]=si
            si =si+1
            i = j
            j = i
        else:
            i = i +1
            j = i

    ms = np.zeros(len(sp))
    for i in range(ms.shape[0]):
        si = sp[i]
        for j in range(ms.shape[0]):
            if j==i:
                continue
            sj = sp[j]

            if si.shape[0] == sj.shape[0]:
                if si == np.fliplr(sj):
                    ms[i] =j

    for i in range(len(sp)):
        sp[i] = sp[i][:,0]








    return sp, ms
