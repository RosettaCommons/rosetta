#!/usr/bin/env python

from numpy import *
from numpy.linalg import *

r_n = array([
    [  97.259, 39.712, 56.443 ],
    [  99.345, 41.629, 56.722 ],
    [ 100.357, 42.410, 59.363 ]])

r_ca = array([
    [  98.448, 39.501, 55.622 ],
    [ 100.335, 42.688, 56.897 ],
    [ 100.816, 42.188, 60.731 ]])

r_c = array([
    [  99.429, 40.636, 55.791 ],
    [ 101.057, 42.540, 58.216 ],
    [  99.665, 42.264, 61.706 ]])

jacobian = zeros((4, 4))

def normalize(array):
    magnitude = sqrt(sum(array**2, axis=1))
    return array / magnitude[:, newaxis]


r_1 = array([r_n[0], r_ca[0], r_n[1], r_ca[1], r_n[2], r_ca[2]])
r_2 = array([r_ca[0], r_c[0], r_ca[1], r_c[1], r_ca[2], r_c[2]])

axis = normalize(r_2 - r_1)
cross_56 = cross(axis[4], axis[5])

for i in range(4):
    distance = r_1[5] - r_1[i]
    direction = axis[i]

    jacobian[0:3, i] = cross(direction, distance)
    jacobian[3, i] = dot(direction, cross_56)

result = 1 / abs(det(jacobian))

print result
