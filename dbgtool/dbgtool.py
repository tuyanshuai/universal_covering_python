import scipy.sparse as sp
import numpy as np
def dump(a):
    if sp.issparse(a):
        np.savetxt("var.csv", a.todense(), delimiter=",")
    else:
        np.savetxt("var.csv", a, delimiter=",")


def dump_sparse(a):
    fp = open("var.csv","w")
    for j in range(a.shape[1]):
        for i in range(a.shape[0]):
            if a[i,j]:
                fp.write("(%d,%d) = %f\n" %(i+1,j+1, a[i,j]))

    fp.close()
