"""
mesh io
"""
import datetime

import numpy


def read_obj(filename):
    with open(filename, "r") as f:
        face, vertex = read_buffer(f)
        return face-1, vertex


def read_buffer(f):
    points = []
    faces = []
    while True:
        line = f.readline()

        if not line:
            # EOF
            break

        strip = line.strip()

        if len(strip) == 0 or strip[0] == "#":
            continue

        split = strip.split()

        if split[0] == "v":
            # vertex
            points.append([numpy.float(item) for item in split[1:]])
        elif split[0] == "vn":
            # skip vertex normals
            pass
        elif split[0] == "s":
            # "s 1" or "s off" controls smooth shading
            pass
        elif split[0] == "f":
            faces.append([int(item.split("//")[0]) for item in split[1:]])
        else:
            # who knows
            pass

    triangle = numpy.array([f for f in faces if len(f) == 3])

    return triangle, numpy.array(points)

#
# def write(filename, mesh):
#     assert (
#         "triangle" in mesh.cells or "quad" in mesh.cells
#     ), "Wavefront .obj files can only contain triangle or quad cells."
#
#     with open(filename, "w") as f:
#         f.write(
#             "# Created by meshio v{}, {}\n".format(
#                 __version__, datetime.datetime.now().isoformat()
#             )
#         )
#         for p in mesh.points:
#             f.write("v {} {} {}\n".format(p[0], p[1], p[2]))
#         if "triangle" in mesh.cells:
#             for c in mesh.cells["triangle"]:
#                 f.write("f {} {} {}\n".format(*(c + 1)))
#         if "quad" in mesh.cells:
#             for c in mesh.cells["quad"]:
#                 f.write("f {} {} {} {}\n".format(*(c + 1)))
#     return
