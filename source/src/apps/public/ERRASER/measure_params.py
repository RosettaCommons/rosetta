#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

from math import sqrt, acos, degrees

def norm(a):
    """
    Norm of vector.
    """
    return sqrt( dot(a,a) )

def cross(a, b):
    """
    Cross product.
    """
    c = [a[1]*b[2] - a[2]*b[1],
         a[2]*b[0] - a[0]*b[2],
         a[0]*b[1] - a[1]*b[0]]
    return c

def dot(a, b):
    """
    Dot product.
    """
    c = a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
    return c

def vec_diff(a, b):
    """
    Return difference vector.
    """
    c = [ a[0]-b[0], a[1]-b[1], a[2]-b[2] ]
    return c

def normalize(a):
    """
    Normalize a vector.
    """
    norm_a = norm(a)
    c = [ a[0] / norm_a, a[1] / norm_a, a[2] / norm_a ]
    return c

def compute_squared_dist(coord1, coord2) :
    """
    compute the squared distance between two xyz vector (3D).
    """
    sq_dist = 0
    sq_dist += (coord1[0] - coord2[0]) * (coord1[0] - coord2[0])
    sq_dist += (coord1[1] - coord2[1]) * (coord1[1] - coord2[1])
    sq_dist += (coord1[2] - coord2[2]) * (coord1[2] - coord2[2])
    return sq_dist

def compute_dist(coord1, coord2) :
    """
    compute the distance between two xyz vector (3D).
    """
    sq_dist = compute_squared_dist(coord1, coord2)
    return sqrt(sq_dist)

def compute_angle(v1,v2):
    """
    calculates the angle between two vectors.
    v1 and v2 are array objects.
    returns a float containing the angle in radians.
    """
    length_product = norm(v1) * norm(v2)
    cosine = dot(v1,v2) / length_product
    angle = degrees( acos( cosine ) )
    return angle

def compute_torsion(v1,v2,v3,v4):
    """
    Returns a float value for the dihedral angle between
    the four vectors. They define the bond for which the
    torsion is calculated (~) as:
    V1 - V2 ~ V3 - V4
    """
    # calculate vectors representing bonds
    v12 = vec_diff(v2, v1)
    v23 = vec_diff(v3, v2)
    v34 = vec_diff(v4, v3)

    # calculate vectors perpendicular to the bonds
    normal1 = normalize( cross(v12,v23) )
    normal2 = normalize( cross(v23,v34) )

    # check for linearity
    if norm(normal1) == 0 or norm(normal2)== 0:
        error_exit("Parallel vectors during torsion computation!!!")

    # calculate torsion and convert to degrees
    torsion = compute_angle(normal1,normal2)

    # take into account the determinant
    # (the determinant is a scalar value distinguishing
    # between clockwise and counter-clockwise torsion.
    if dot(normal1,v34) < 0:
        torsion = 360.0 - torsion

    if (torsion >= 180.0) : torsion -= 360.0
    return torsion

