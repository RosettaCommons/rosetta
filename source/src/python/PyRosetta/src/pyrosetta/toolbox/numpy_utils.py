#!/usr/bin/env python
# :noTabs=true:


# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

# @file   numpy_utils.py
# @brief
# @author Brian D. Weitzner (brian.weitzner@gmail.com)

from __future__ import print_function


import numpy as np
from math import acos, degrees


def get_coordinates_from_pose(p, resNo=1):
    """Pack atomic coordinates in a single residue into a numpy matrix.

    Args:
        p (Pose): Pose from coordinates will be extracted
        resNo (int): the pose-numbered residue of interest

    Returns:
        np.ndarray: Set of coordinates for the atoms in the specified residue in the Pose

    """
    assert(resNo > 0 and resNo <= p.size())
    res = p.residues[resNo]
    return np.array(res.xyz(i) for i in range(1, res.natoms() + 1)])


# the following is taken from http://nghiaho.com/?page_id=671
def rigid_transform_3D(A, B):
    """Compute the translation and rotation to optimally superpose two sets
    of points using singular value decomposition.

    Args:
        A (np.ndarray): A set of points
        B (np.ndarray): A set of points (must have a one-to-one correspondence
            to the points in A)

    Returns:
        tuple: Rotation matrix (np.ndarray) and translation vector (np.array) to
            optimally superpose A onto B
    """
    assert(A.shape == B.shape)

    N = A.shape[0]  # total points

    centroid_A = np.mean(A, axis=0)
    centroid_B = np.mean(B, axis=0)

    # center the points
    AA = A - np.tile(centroid_A, (N, 1))
    BB = B - np.tile(centroid_B, (N, 1))

    # construct matrix from coords on which SVD will be performed
    H = np.dot(np.transpose(AA), BB)

    U, S, Vt = np.linalg.svd(H)

    R = np.dot(Vt.T, U.T)

    # special reflection case
    if np.linalg.det(R) < 0:
        print('Reflection detected')
        Vt[2, :] *= -1
        R = Vt.T * U.T

    t = -R * centroid_A.T + centroid_B.T
    return R, t


def calc_dihedral(c):
    """Compute a dihedral angle from four points in space.

    Args:
        c (np.ndarray): Coordinates of four points in space (atoms)

    Returns:
        float: The dihedral angle
    """
    from numpy.linalg import norm
    # make sure there are 4 coordinates in three-space
    assert(np.shape(c) == (4, 3))

    # make vectors pointing from each atom to the next bonded atom
    v_12 = c[1] - c[0]
    v_23 = c[2] - c[1]
    v_34 = c[3] - c[2]

    # calculate the two planes that can be formed by these four atoms
    p_123 = np.cross(v_12, v_23)
    p_234 = np.cross(v_23, v_34)

    # find the angle between the planes
    chi = degrees(acos(np.dot(p_123, p_234.T) / (norm(p_123) *
                  norm(p_234))))

    # the sign of torsion angles is determined by whether the vector formed
    # by the last two points is above (positive) or below (negative) the
    # plane defined by the first three points.
    sign = 1. if np.dot(v_34, p_123.T) > 0. else -1.

    return sign * chi


def translate_pose(p, t):
    """Apply a translation to all of the coordinates in a Pose.

    Args:
        p (Pose): The Pose instance to manipulate
        t (np.array): A vector to add to the Pose coordinates
    """
    assert(t.shape in ((3,), (4,)))
    xform = np.identity(4)
    xform[:t.shape[0], 3] = t
    p.apply_transform(xform)


def rotate_pose(p, R):
    """Apply a rotation matrix to all of the coordinates in a Pose.

    Args:
        p (Pose): The Pose instance to manipulate
        R (np.mat): A rotation matrix to apply to the Pose coordinates
    """
    assert(R.shape in ((3,3), (4,3)))
    xform = np.identity(4)
    xform[:R.shape[0], :3] = R
    p.apply_transform(xform)
