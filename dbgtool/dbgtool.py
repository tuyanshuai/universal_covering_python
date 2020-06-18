import scipy.sparse as sp
import numpy as np
def dump(a):
    if sp.issparse(a):
        np.savetxt("var.csv", a.todense(), delimiter=",")
    else:
        np.savetxt("var.csv", a, delimiter=",")