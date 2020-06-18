# E is extra feature for vertex only

def write_mfile(fn, F, V):

    with open(fn, 'w') as f:
        for i in range(V.shape[0]):
            f.write("Vertex %d %.4f %.4f %.4f\n" % (i+1,V[i,0],V[i,1],V[i,2]))
        for i in range(F.shape[0]):
            f.write("Face %d %d %d %d\n" % (i + 1, F[i, 0], F[i, 1], F[i, 2]))
        f.write("\n")
