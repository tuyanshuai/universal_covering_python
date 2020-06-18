# read off file


import numpy
def read_off(filename):
    file = open(filename)
    if 'OFF' != file.readline().strip():
        raise('Not a valid OFF header')
    n_verts, n_faces, n_e  = tuple([int(s) for s in file.readline().strip().split(' ')])
    points = [[float(s) for s in file.readline().strip().split(' ')] for i_vert in range(n_verts)]
    faces = [[int(s) for s in file.readline().strip().split(' ')][1:] for i_face in range(n_faces)]
    file.close()

    triangle = numpy.array([f for f in faces if len(f) == 3])
    return triangle+1, numpy.array(points)
