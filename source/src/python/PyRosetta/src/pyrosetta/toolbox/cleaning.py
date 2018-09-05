#!/usr/bin/env python
# :noTabs=true:


# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   cleaning.py
## @brief
## @author Evan H. Baugh, Johns Hopkins University

from __future__ import print_function

import os


def cleanATOM(pdb_file, out_file=None, ext=".clean.pdb"):
    """Extract all ATOM and TER records in a PDB file and write them to a new file.

    Args:
        pdb_file (str): Path of the PDB file from which ATOM and TER records
            will be extracted
        out_file (str): Optional argument to specify a particular output filename.
            Defaults to <pdb_file>.clean.pdb.
        ext (str): File extension to use for output file. Defaults to ".clean.pdb"
    """
    # find all ATOM and TER lines
    with open(pdb_file, "r") as fid:
        good = [l for l in fid if l.startswith(("ATOM", "TER"))]

    # default output file to <pdb_file>.clean.pdb
    if out_file is None:
        out_file = os.path.splitext(pdb_file)[0] + ext

    # write the selected records to a new file
    with open(out_file, "w") as fid:
        fid.writelines(good)


def cleanCRYS(pdb_file, olig=2, out_file=None):
    """Extract a monomer from an oligomeric PDB file and write it to a new file.

    Notes:
        This is a simple sequence comparison.

    Args:
        pdb_file (str): Path of the PDB file from which ATOM and TER records
            will be extracted.
        olig (int): Number of subunits in the input complex.
        out_file (str): Optional argument to specify a particular output filename.
            Defaults to <pdb_file>.mono.pdb.
    """
    # load in the PDB as a pose to get the sequence
    from pyrosetta import pose_from_file

    pose = pose_from_file(pdb_file)
    seq = pose.sequence()

    # split sequence into N evenly-sized fragments
    frac = int(round(pose.total_residue() / olig))
    frags = [seq[i : i + frac] for i in range(0, len(seq), frac)]

    # check that there are N identical sequence fragments
    try:
        assert len(frags) == olig and len(set(frags)) == 1
    except AssertionError:
        raise ValueError('"' + pdb_file + '" is not a ' + str(olig) + "-mer")

    [pose.delete_polymer_residue(frac + 1) for _ in range(frac * (olig - 1))]

    # write the new pdb file
    if out_file is None:
        out_file = os.path.splitext(pdb_file)[0] + ".mono.pdb"

    pose.dump_pdb(out_file)
