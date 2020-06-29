import scipy.sparse as sp
import numpy as np
import time
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


t =0
def tic():
    global t
    t = time.time()
    return None

def toc():
    # do stuff
    global t
    elapsed = time.time() - t
    print("It ellpased %f(s)" % elapsed)
    return elapsed