'''
convert a mesh obj to points
'''
import numpy as np
from meshio import *
import sys

# compute triangle mesh surface area
def triangle_area(x):
    a = x[:, 0, :] - x[:, 1, :]
    b = x[:, 0, :] - x[:, 2, :]
    cross = np.cross(a, b)
    area = 0.5 * np.linalg.norm(np.cross(a, b), axis=1)
    return area


# compute euclidean distance matrix
def euclidean_distance_matrix(x):
    r = np.sum(x * x, 1)
    r = r.reshape(-1, 1)
    distance_mat = r - 2 * np.dot(x, x.T) + r.T
    # return np.sqrt(distance_mat)
    return distance_mat


# update distance matrix and select the farthest point from set S after a new point is selected
def update_farthest_distance(far_mat, dist_mat, s):
    for i in range(far_mat.shape[0]):
        far_mat[i] = dist_mat[i, s] if far_mat[i] > dist_mat[i, s] else far_mat[i]
    return far_mat, np.argmax(far_mat)


# initialize matrix to keep track of distance from set s
def init_farthest_distance(far_mat, dist_mat, s):
    for i in range(far_mat.shape[0]):
        far_mat[i] = dist_mat[i, s]
    return far_mat


# get sample from farthest point on every iteration
def farthest_point_sampling(off_file, num_samples=1000):
    triangle, vertices= read_off(off_file)  # pymesh is used to load mesh obj files.
    faces = vertices[triangle,:]
    area = triangle_area(faces)
    total_area = np.sum(area)

    set_P = []
    for i in range(faces.shape[0]):
        num_gen = area[i] / total_area * 10000
        for j in range(int(num_gen) + 1):
            r1, r2 = np.random.rand(2)
            d = (1 - np.sqrt(r1)) * faces[i, 0] + np.sqrt(r1) * (1 - r2) * faces[i, 1] + np.sqrt(r1) * r2 * faces[i, 2]
            set_P.append(d)

    set_P = np.array(set_P)
    num_P = set_P.shape[0]

    distance_mat = euclidean_distance_matrix(set_P)

    set_S = []
    s = np.random.randint(num_P)
    far_mat = init_farthest_distance(np.zeros((num_P)), distance_mat, s)

    for i in range(num_samples):
        set_S.append(set_P[s])
        far_mat, s = update_farthest_distance(far_mat, distance_mat, s)

    return np.array(set_S)


def write(fname, pts):
    with open(fname, 'w') as f:
        f.writelines(','.join(str(xyz) for xyz in pt) + '\n' for pt in pts)


def mesh2pts(meshfn, outfn):
    sample_pts = farthest_point_sampling(meshfn)
    write(outfn, sample_pts)

# main
if __name__ == "__main__":
    sample_pts = farthest_point_sampling(sys.argv[1])
    write(sys.argv[2], sample_pts)