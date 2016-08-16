#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   extract_coords_pose.py
## @brief
## @author Evan H. Baugh, Johns Hopkins University

from rosetta import Pose

# returns the x, y, and z coordinates of a pose
#     returns a 1x3 list with elements [x,y,z]
#     optional selection list of desired residues
def extract_coordinates_from_pose_3x1( pose , selection = [] ,
        atom_names = [] , not_atom_names = [] ):
    """
    Returns a list of corresponding x, y, and z coordinates from  <pose>
        for all residues in selection (defaults to all residues)
        with atom names in  <atom_names> and not in  <not_atom_names>
    note: returns a list with three elements, the x, y, and z

    example(s):
        xyz = extract_coordinates_from_pose(pose)
    See Also:
        extract_coordinates_from_pdb
        pose_from_file
        mean
        moment
        vector_moment
    """
    # default to all residues
    if not selection:
        selection = range( 1 , pose.total_residue() + 1 )
    # empty coordinate holders
    coords = [[] for i in range(3)]
    # for each residue
    for resi in selection:
        resi = pose.residue(resi)
        # for each atom
        for atom in range( resi.natoms() ):
            # check if the atom is a type desired
            if ( not atom_names or resi.atom_name( atom + 1 ).strip()
                    in atom_names ) and ( not
                    resi.atom_name( atom + 1 ).strip() in not_atom_names ):
                atom = resi.xyz( atom + 1 )
                # store the coordinates
                coords[0].append( atom[0] )    # x
                coords[1].append( atom[1] )    # y
                coords[2].append( atom[2] )    # z
    return coords

# returns the x, y, and z coordinates of a pose
#     returns a list of elements [x,y,z]
#     optional selection list of desired residues
def extract_coordinates_from_pose_1x3( pose , selection = [] ,
        atom_names = [] , not_atom_names = [] ):
    """
    Returns a list of corresponding x, y, and z coordinates from  <pose>
        for all residues in selection (defaults to all residues)
        with atom names in  <atom_names> and not in  <not_atom_names>
    note: returns a list with three elements, the x, y, and z

    example(s):
        xyz = extract_coordinates_from_pose(pose)
    See Also:
        extract_coordinates_from_pdb
        pose_from_file
        mean
        moment
        vector_moment
    """
    # extract the default, 3 vectors
    xyz = extract_coordinates_from_pose_3x1( pose , selection ,
        atom_names , not_atom_names )
    # change to one vector of 3 vectors
    xyz = [ [xyz[0][i],xyz[1][i],xyz[2][i]] for i in range( len( xyz[0] ) ) ]
    return xyz

# methods to add
# return list of xyzVector objects
# return dictionary list w/ chains or residues
# write coordinates


