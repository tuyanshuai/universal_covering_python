"""
Compute face normal and face area
 [fn,fa] = calculate_face_normal(face,vertex)
"""
import numpy as np
from scipy.linalg import norm


def calculate_face_normal(face, vertex):
    x = vertex[face[:, 1], :] - vertex[face[:, 0], :]
    y = vertex[face[:, 2], :] - vertex[face[:, 0], :]
    fn = np.cross(x, y, axis=1)
    fa = norm(fn,axis=1)
    fn[:, 0] = np.divide(fn[:, 0], fa)
    fn[:, 1] = np.divide(fn[:, 1], fa)
    fn[:, 2] = np.divide(fn[:, 2], fa)

    return fn
