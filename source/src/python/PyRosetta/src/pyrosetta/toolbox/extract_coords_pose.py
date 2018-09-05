#!/usr/bin/env python
# :noTabs=true:


# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   extract_coords_pose.py
## @brief
## @author Evan H. Baugh, Johns Hopkins University


def pose_coords_as_rows(pack_or_pose, selection=list(), atom_names=list()):
    """Construct an ndarray of x, y, and z coordinates from all residues in selection
        and atom names in <atom_names> from <pose> as rows and return it.

    Args:
        pack_or_pose (pyrosetta.rosetta.core.pose.Pose OR pyrosetta.distributed.packed_pose.PackedPose):
            the `Pose` instance from which coordinates will be extracted.
        selection (list): a list of residue numbers to use. Defaults to all residues.
        atom_names (list): a list of atom names to include. Defaults to all atoms.

    Returns:
        numpy.ndarray: an Nx3 array of coordinates.
    """
    import pyrosetta.distributed.packed_pose as packed_pose
    import numpy as np

    wpose = packed_pose.to_pose(pack_or_pose)

    # default to all residues
    if not selection:
        selection = range(1, wpose.total_residue() + 1)

    coords = np.array(
        [
            res.xyz(atom)
            for res in [wpose.residues[s] for s in selection]
            for atom in range(1, res.natoms() + 1)
            if (not atom_names or res.atom_name(atom).strip() in atom_names)
        ]
    )
    return coords


def pose_coords_as_cols(pack_or_pose, selection=list(), atom_names=list()):
    """Construct an ndarray of x, y, and z coordinates from all residues in selection
        and atom names in <atom_names> from <pose> as columns and return it.

    Args:
        pack_or_pose (pyrosetta.rosetta.core.pose.Pose OR pyrosetta.distributed.packed_pose.PackedPose):
            the `Pose` instance from which coordinates will be extracted.
        selection (list): a list of residue numbers to use. Defaults to all residues.
        atom_names (list): a list of atom names to include. Defaults to all atoms.

    Returns:
        numpy.ndarray: a 3xN array of coordinates.
    """
    return extract_coordinates_from_pose_as_rows(pack_or_pose, selection, atom_names).T
