# E is extra feature for vertex only
import numpy as np


def write_mfile(fn, F, V, C=np.zeros(0)):
    with open(fn, 'w') as f:
        for i in range(V.shape[0]):
            if C.shape[0] == 0:
                f.write("Vertex %d %.4f %.4f %.4f\n" % (i + 1, V[i, 0], V[i, 1], V[i, 2]))
            else:
                f.write("Vertex %d %.4f %.4f %.4f {rgb=( %.4f %.4f %.4f )}\n" % (
                    i + 1, V[i, 0], V[i, 1], V[i, 2], C[i, 0], C[i, 1], C[i, 2]))
        for i in range(F.shape[0]):
            f.write("Face %d %d %d %d\n" % (i + 1, F[i, 0], F[i, 1], F[i, 2]))
        f.write("\n")
