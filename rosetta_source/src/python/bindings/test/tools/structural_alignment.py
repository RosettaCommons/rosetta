#!/usr/bin/env python
# :noTabs=true:


# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   structural_alignment.py
## @brief
## @author Evan H. Baugh, Johns Hopkins University
 
##############################################################################
#
# @SUMMARY: -- QKabsch.py.  A python implementation of the optimal superposition
#     of two sets of vectors as proposed by Kabsch 1976 & 1978.
#
# @AUTHOR: Jason Vertrees
# @COPYRIGHT: Jason Vertrees (C), 2005-2007
# @LICENSE: Released under GPL:
# This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
# Street, Fifth Floor, Boston, MA 02110-1301, USA 
#
# DATE  : 2007-01-01
# REV   : 2
# REQUIREMENTS: numpy
#
#############################################################################
# script modified by Evan H. Baugh for usage with PyRosetta
#
# DATA  : 2011-06-20
# REV   : 1
# REQUIREMENTS: numpy, PyRosetta

# using numpy for linear algebra
import numpy

from rosetta import Pose
from rosetta.numeric import xyzVector
 
# other tools
from extract_coords_pose import extract_coordinates_from_pose_1x3

def kabsch_alignment( pose1, pose2 , pose1sel = [], pose2sel =[] ):
    """
    optAlign performs the Kabsch alignment algorithm upon the alpha-carbons of two selections.
  
    By default, this program will optimally align the ALPHA CARBONS of the selections provided.
    To turn off this feature remove the lines between the commented "REMOVE ALPHA CARBONS" below.
 
    This is a hacky modification for usage in PyRosetta,
    The geometry/math is the same but it uses poses instead of PyMOL
 
    @param pose1: First pose
    @param pose2: Second pose
    @param pose1sel: List of residues to align in the first pose
    @param pose2sel: List of residues to align in the second pose
    """	

    # default to all the residues
    if not pose1sel:
        pose1sel = range(1,pose1.total_residue()+1)
    if not pose2sel:
        pose2sel = range(1,pose2.total_residue()+1)

    # obtain a list of the atom positions to align
    stsel1 = extract_coordinates_from_pose_1x3(pose1,pose1sel,['CA'])
    stsel2 = extract_coordinates_from_pose_1x3(pose2,pose2sel,['CA'])
    """
    stsel1 = []
    for r in pose1sel:
        v = pose1.residue(r).atom('CA').xyz()
        stsel1.append([v[0],v[1],v[2]])
    stsel2 = []
    for r in pose2sel:
        v = pose2.residue(r).atom('CA').xyz()
        stsel2.append([v[0],v[1],v[2]])
    """
    # remember all of each pose's atom positions, to apply transformation
    #    later for updating
    molsel1 = extract_coordinates_from_pose_1x3(pose1)
    molsel2 = extract_coordinates_from_pose_1x3(pose2)
    """
    molsel1 = []
    for r in range(pose1.total_residue()):
        r = pose1.residue(r+1)
        for a in range(r.natoms()):
            v = r.xyz(a+1)
            molsel1.append([v[0],v[1],v[2]])
    molsel2 = []
    for r in range(pose2.total_residue()):
        r = pose2.residue(r+1)
        for a in range(r.natoms()):
            v = r.xyz(a+1)
            molsel2.append([v[0],v[1],v[2]])
    """
    
    # check for consistency, same number of selected residues
    assert len(stsel1) == len(stsel2)
    L = len(stsel1)
    assert L > 0

    # must alway center the two proteins to avoid affine transformations
    #     Center the two proteins to their selections
    COM1 = numpy.sum(stsel1,axis=0) / float(L)
    COM2 = numpy.sum(stsel2,axis=0) / float(L)
    stsel1 -= COM1
    stsel2 -= COM2

    # Initial residual, see Kabsch.
    E0 = numpy.sum( numpy.sum(stsel1 * stsel1,axis=0),axis=0) + numpy.sum( numpy.sum(stsel2 * stsel2,axis=0),axis=0)

    # This beautiful step provides the answer.  V and Wt are the orthonormal
    # bases that when multiplied by each other give us the rotation matrix, U.
    # S, (Sigma, from SVD) provides us with the error!  Isn't SVD great!
    V, S, Wt = numpy.linalg.svd( numpy.dot( numpy.transpose(stsel2), stsel1))

    # we already have our solution, in the results from SVD.
    # we just need to check for reflections and then produce
    # the rotation.  V and Wt are orthonormal, so their det's
    # are +/-1.
    reflect = float(str(float(numpy.linalg.det(V) * numpy.linalg.det(Wt))))

    if reflect == -1.0:
        S[-1] = -S[-1]
        V[:,-1] = -V[:,-1]

    RMSD = E0 - (2.0 * sum(S))
    RMSD = numpy.sqrt(abs(RMSD / L))

    #U is simply V*Wt
    U = numpy.dot(V, Wt)

    # rotate and translate the molecule
    stsel2 = numpy.dot((molsel2 - COM2), U)
    stsel2 = stsel2.tolist()
    # center the molecule
    stsel1 = molsel1 - COM1
    stsel1 = stsel1.tolist()

    # apply the changes to both poses
    # first pose
    ind = -1
    #dummy = rosetta.numeric.xyzVector()
    dummy = xyzVector()
    for r in range(pose1.total_residue()):
        res = pose1.residue(r+1)
        for a in range(res.natoms()):
            ind += 1
            v = stsel1[ind]
            dummy.x = v[0]
            dummy.y = v[1]
            dummy.z = v[2]
            pose1.residue(r+1).set_xyz(a+1,dummy)
    # second pose
    ind = -1
    for r in range(pose2.total_residue()):
        res = pose2.residue(r+1)
        for a in range(res.natoms()):
            ind += 1
            v = stsel2[ind]
            dummy.x = v[0]
            dummy.y = v[1]
            dummy.z = v[2]
            pose2.residue(r+1).set_xyz(a+1,dummy)

    print 'RMSD=%f' % RMSD


