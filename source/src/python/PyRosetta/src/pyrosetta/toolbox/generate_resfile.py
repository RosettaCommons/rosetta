#!/usr/bin/env python
# :noTabs=true:


# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   generate_resfile.py
## @brief
## @author Evan H. Baugh, Johns Hopkins University

# adapted from original code by Sid Chaudhury

from __future__ import print_function


def generate_resfile_from_pose(
    pack_or_pose,
    resfilename,
    pack=True,
    design=False,
    input_sc=True,
    freeze=list(),
    specific=dict(),
):
    """Generate a resfile from a pose and write it to a file.

    Args:
        pack_or_pose (pyrosetta.rosetta.core.pose.Pose OR pyrosetta.distributed.packed_pose.PackedPose):
            the `Pose` instance to use to for resfile creation.

        pack (bool): True allows packing, False disallows. Defaults to True.

        design (bool): True allows design using all amino acids, False disables design.
            Defaults to False.

        input_sc (bool): Toggles the inclusion of the original side-chain conformations.
            Defaults to True.

        freeze (list): An optional list of residue numbers to exclude from packing and/or
            design.

        specific (dict): An optional dictionary with residue numbers and resfile keywords
            as key-value pairs.
    """
    # determine the header, default settings
    to_write = (
        (
            "NATAA\n"
            if not design
            else "ALLAA\n# ALLAA will NOT work on bridged Cysteines\n"
        )
        if pack
        else "NATRO\n"
    )
    to_write += ("USE_INPUT_SC\n" if input_sc else "") + "start\n"

    # use frozen list to update dict
    specific.update({resNo: "NATRO" for resNo in freeze})

    # use PDB numbering if possible
    import pyrosetta.distributed.packed_pose as packed_pose

    wpose = packed_pose.to_pose(pack_or_pose)
    info = wpose.pdb_info()
    if info and info.nres():
        res = [
            (info.number(resNo), info.chain(resNo), keyword)
            for resNo, keyword in specific.items()
        ]
    else:
        res = [(resNo, " ", keyword) for resNo, keyword in specific.items()]
    to_write += "\n".join(
        [
            str(num).rjust(4) + str(chain).rjust(3) + "  " + keyword + "  "
            for num, chain, keyword in res
        ]
    )

    with open(resfilename, "w") as f:
        f.write(to_write)


def generate_resfile_from_pdb(
    pdbfilename,
    resfilename,
    pack=True,
    design=False,
    input_sc=True,
    freeze=list(),
    specific=dict(),
):
    """Generate a resfile from a PDB file and write it to a file.

     Args:
         pdbfilename (str): the PDB filename to use to for resfile creation.

         pack (bool): True allows packing, False disallows. Defaults to True.

         design (bool): True allows design using all amino acids, False disables design.
            Defaults to False.

         input_sc (bool): Toggles the inclusion of the original side-chain conformations.
            Defaults to True.

        freeze (list): An optional list of residue numbers to exclude from packing and/or
            design.

        specific (dict): An optional dictionary with residue numbers and resfile keywords
            as key-value pairs.
    """
    from pyrosetta import pose_from_file

    p = pose_from_file(pdbfilename)
    generate_resfile_from_pose(p, resfilename, pack, design, input_sc, freeze, specific)
