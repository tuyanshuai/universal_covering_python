import scipy.sparse as sp
import numpy as np
import time

# dump a simple real matrix
def simdump(fn, a):
    np.savetxt(fn, a, delimiter=",",fmt='%60.60f')

def dump(a):

    if type(a)==list:
        if len(a[0]) > 1:  # array of array
            for i in range(len(a)):
                if np.any(np.iscomplex(a[0])):
                    fn = "%d_real.csv" % (i+1)
                    simdump(fn, np.real(a[i]))
                    fn = "%d_imag.csv" % (i+1)
                    simdump(fn, np.imag(a[i]))
                else:
                    fn = "%d.csv" % (i+1)
                    simdump(fn, a[i])
        else:
            simdump("var.csv", a)  # ordinary array
    else:
        if sp.issparse(a):
            simdump("var.csv", a.todense())
        else:
            if np.any(np.iscomplex(a)):
                simdump("real.csv", a.real)
                simdump("imag.csv", a.imag)
            else:
                simdump("var.csv", a)  # ordinary array


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
    global t
    elapsed = time.time() - t
    print("It ellpased %f(s)" % elapsed)
    return elapsed