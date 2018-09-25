#!/usr/bin/env python
# :noTabs=true:


# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   mutants.py
## @brief
## @author Evan H. Baugh, Johns Hopkins University

# mutate_residue is adapted from an original script by Sid Chaudhury

# UNFINISHED...below!

from __future__ import print_function

import random


def restrict_non_nbrs_from_repacking(pose, res, task, pack_radius):
    """Configure a `PackerTask` to only repack neighboring residues and
    return the task.

    Args:
        pose (pyrosetta.Pose): The `Pose` to opertate on.
        res (int): Pose-numbered residue position to exclude.
        task (pyrosetta.rosetta.core.pack.task.PackerTask): `PackerTask` to modify.
        pack_radius (float): Radius used to define neighboring residues.

    Returns:
        pyrosetta.rosetta.core.pack.task.PackerTask: Configured `PackerTask`.
    """

    def representative_coordinate(resNo):
        return pose.residue(resNo).xyz(pose.residue(resNo).nbr_atom())

    center = representative_coordinate(res)
    for i in range(1, len(pose.residues) + 1):
        # only pack the mutating residue and any within the pack_radius
        if i == res:
            continue

        if center.distance(representative_coordinate(i)) > pack_radius:
            task.nonconst_residue_task(i).prevent_repacking()
        else:
            task.nonconst_residue_task(i).restrict_to_repacking()
    return task


def mutate_residue(
    pack_or_pose, mutant_position, mutant_aa, pack_radius=0., pack_scorefxn=None
):
    """Replace the residue at a single position in a Pose with a new amino acid
        and repack any residues within user-defined radius of selected residue's
        center using.

    Args:
        pack_or_pose (pyrosetta.rosetta.core.pose.Pose OR pyrosetta.distributed.packed_pose.PackedPose):
                    the `Pose` instance to use.
        mutant_position (int): Pose-numbered position of the residue to mutate.
        mutant_aa (str): The single letter name for the desired amino acid.
        pack_radius (float): Radius used to define neighboring residues.
        pack_scorefxn (pyrosetta.ScoreFunction): `ScoreFunction` to use when repacking the `Pose`.
            Defaults to the standard `ScoreFunction`.
    """
    import pyrosetta
    import pyrosetta.distributed.packed_pose as packed_pose

    wpose = packed_pose.to_pose(pack_or_pose)

    if not wpose.is_fullatom():
        raise IOError("mutate_residue only works with fullatom poses")

    # create a standard scorefxn by default
    if not pack_scorefxn:
        pack_scorefxn = pyrosetta.get_score_function()

    # the numbers 1-20 correspond individually to the 20 proteogenic amino acids
    from pyrosetta.rosetta.core.chemical import aa_from_oneletter_code

    mutant_aa = int(aa_from_oneletter_code(mutant_aa))
    aa_bool = pyrosetta.Vector1([aa == mutant_aa for aa in range(1, 21)])

    # mutation is performed by using a PackerTask with only the mutant
    # amino acid available during design
    task = pyrosetta.standard_packer_task(wpose)
    task.nonconst_residue_task(mutant_position).restrict_absent_canonical_aas(aa_bool)

    # prevent residues from packing by setting the per-residue "options" of the PackerTask
    task = restrict_non_nbrs_from_repacking(wpose, mutant_position, task, pack_radius)

    # apply the mutation and pack nearby residues
    from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover

    packer = PackRotamersMover(pack_scorefxn, task)
    packer.apply(wpose)
