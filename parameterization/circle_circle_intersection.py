import numpy as np
from math import *


def circle_circle_intersection(c1, r1, c2, r2):
    p = 0
    dz = c2 - c1
    d = abs(dz)
    a = (r1 * r1 - r2 * r2 + d * d) / d / 2.0
    z0 = c1 + dz * a / d
    h = sqrt(r1 * r1 - a * a)

    rz = 1j * dz * h / d
    p1 = z0 + rz
    p2 = z0 - rz
    e1 = c2 - c1
    e2 = p1 - c1
    ind = (e1.real * e2.imag - e1.imag * e2.real) > 0
    if ind:
        p = p1
    else:
        p = p2
    if d > r1 + r2 or d < abs(r2 - r1):
        p = nan
    return p
