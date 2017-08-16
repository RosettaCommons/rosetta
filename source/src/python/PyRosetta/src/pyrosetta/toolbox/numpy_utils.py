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


from math import acos, degrees
import pyrosetta.rosetta as rosetta

class NotInstalled(object):
    def __init__(self, name):
        self._name = name

    def __getattr__(self, item):
        raise ImportError('The {0} package is required to use this '
                          'feature'.format(self._name))


try:
    import numpy as np
    from numpy.linalg import norm
except ImportError:
    np = NotInstalled(name='numpy')
    norm = NotInstalled(name='numpy.linalg.norm')


def get_coordinates_from_pose(p, resNo=1):
    '''Pack atomic coordinates in a single residue into a numpy matrix.

    Args:
        p (Pose): Pose from coordinates will be extracted
        resNo (int): the pose-numbered residue of interest

    Returns:
        Set of coordinates (np.mat) for the atoms in the Pose

    '''
    assert(resNo > 0 and resNo <= p.size())
    res = p.residue(resNo)

    # only iterate over relevant atoms
    coords = []
    for i in atoms[p.residue(1).name()]:
        coords.append([res.xyz(i).x, res.xyz(i).y, res.xyz(i).z])

    return np.mat(coords)


# the following is taken from http://nghiaho.com/?page_id=671
def rigid_transform_3D(A, B):
    '''Compute the translation and rotation to optimally superpose two sets
    of points using singular value decomposition.

    Args:
        A (np.mat): A set of points
        B (np.mat): A set of points (must have a one-to-one correspondence
            to the points in A)

    Returns:
        Rotation matrix (np.mat) and translation vector (np.array) to
            optimally superpose A onto B

    '''
    assert len(A) == len(B)

    N = A.shape[0]  # total points

    centroid_A = np.mean(A, axis=0)
    centroid_B = np.mean(B, axis=0)

    # center the points
    AA = A - np.tile(centroid_A, (N, 1))
    BB = B - np.tile(centroid_B, (N, 1))

    # construct matrix from coords on which SVD will be performed
    H = np.transpose(AA) * BB

    U, S, Vt = np.linalg.svd(H)

    R = Vt.T * U.T

    # special reflection case
    if np.linalg.det(R) < 0:
        print('Reflection detected')
        Vt[2, :] *= -1
        R = Vt.T * U.T

    t = -R * centroid_A.T + centroid_B.T
    return R, t


def numpy_to_rosetta(np_arr):
    '''Pack values in a numpy data structure into the analogous Rosetta
    data structure.

    Args:
        np_arry (np.mat or np.array): One- (1 x 3) or two-dimensional (3 x 3)
            numpy matrix

    Returns:
        The values in the appropriate Rosetta container
            (numeric.xyzVector_double_t or numeric.xyzMatrix_double_t)

    '''
    # start off assuming a 1D array
    dim = 1

    # ensure that we are in 3D space
    if np.shape(np_arr)[1] != 3:
        dim = -999
    elif np.shape(np_arr)[0] == np.shape(np_arr)[1]:
        dim = 2

    if dim == 1:
        # handle the 1D case
        ros_container = rosetta.numeric.xyzVector_double_t(0.)
        ros_container.x = np_arr[0, 0]
        ros_container.y = np_arr[0, 1]
        ros_container.z = np_arr[0, 2]

        return ros_container

    elif dim == 2:
        # handle the 2D case
        ros_container = rosetta.numeric.xyzMatrix_double_t(0.)
        ros_container.xx = np_arr[0, 0]
        ros_container.xy = np_arr[0, 1]
        ros_container.xz = np_arr[0, 2]

        ros_container.yx = np_arr[1, 0]
        ros_container.yy = np_arr[1, 1]
        ros_container.yz = np_arr[1, 2]

        ros_container.zx = np_arr[2, 0]
        ros_container.zy = np_arr[2, 1]
        ros_container.zz = np_arr[2, 2]

        return ros_container

    # get out of here!
    raise ValueError('Packing {}-dimensional numpy arrays '.format(dim) +
                     'into Rosetta containers is not currently supported')


def calc_dihedral(c):
    '''Compute a dihedral angle from four points in space.

    Args:
        c (np.mat): Coordinates of four points in space (atoms)

    Returns:
        The dihedral angle as a float

    '''
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
    '''Apply a translation to all of the coordinates in a Pose.

    Args:
        p (Pose): The Pose instance to manipulate
        t (np.array): A vector to add to the Pose coordinates

    Returns:
        None. The input Pose is manipulated.

    '''
    # t must be an xyzVector_double_t
    for i in range(1, p.size() + 1):
        for j in range(1, p.residue(i).natoms() + 1):
            p.residue(i).atom(j).xyz(p.residue(i).atom(j).xyz() + t)


def rotate_pose(p, R):
    '''Apply a rotation matrix to all of the coordinates in a Pose.

    Args:
        p (Pose): The Pose instance to manipulate
        R (np.mat): A rotation matrix to apply to the Pose coordinates

    Returns:
        None. The input Pose is manipulated.

    '''
    # t must be an xyzMatrix_Real
    for i in range(1, p.size() + 1):
        for j in range(1, p.residue(i).natoms() + 1):
            v = p.residue(i).atom(j).xyz()

            x = R.xx * v.x + R.xy * v.y + R.xz * v.z
            y = R.yx * v.x + R.yy * v.y + R.yz * v.z
            z = R.zx * v.x + R.zy * v.y + R.zz * v.z

            p.residue(i).atom(j).xyz(rosetta.numeric.xyzVector_double_t(x, y, z))
