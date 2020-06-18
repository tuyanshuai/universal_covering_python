import numpy as np

def magnitude(vector):
    return np.sqrt(np.dot(np.array(vector), np.array(vector)))

def normalize(vector):
    return np.array(vector) / magnitude(np.array(vector))


def lineRayIntersectionPoint(rayOrigin, rayDirection, point1, point2):
    # Convert to numpy arrays
    rayOrigin = np.array(rayOrigin, dtype=np.float)
    rayDirection = np.array(normalize(rayDirection), dtype=np.float)
    point1 = np.array(point1, dtype=np.float)
    point2 = np.array(point2, dtype=np.float)


    v1 = rayOrigin - point1
    v2 = point2 - point1
    v3 = np.array([-rayDirection[1], rayDirection[0]])
    t1 = np.cross(v2, v1) / np.dot(v2, v3)
    t2 = np.dot(v1, v3) / np.dot(v2, v3)
    if t1 >= 0.0 and t2 >= 0.0 and t2 <= 1.0:
        return rayOrigin + t1 * rayDirection, t2
    return None,  None


def intersectRayPolygon(origin, direction, polygon):

    nseg = polygon.shape[0] -1
    intersects = list()
    t2s = list()
    for i in range(nseg):
        p0 = polygon[i, :]
        p1 = polygon[i+1, :]
        intersect, t2 = lineRayIntersectionPoint(origin, direction, p0, p1)

        if intersect is None:
            continue
        else:
             return intersect




#
# def intersectRayPolygon(origin, direction, polygon):
#
#     nseg = polygon.shape[0] -1
#     intersects = list()
#     t2s = list()
#     for i in range(nseg):
#         p0 = polygon[i, :]
#         p1 = polygon[i+1, :]
#         intersect, t2 = lineRayIntersectionPoint(origin, direction, p0, p1)
#
#         if intersect is None:
#             continue
#         else:
#             intersects.append(intersect)
#             t2s.append(abs(t2-0.5))
#
#     id = t2s.index(min(t2s))
#     return intersects[id]