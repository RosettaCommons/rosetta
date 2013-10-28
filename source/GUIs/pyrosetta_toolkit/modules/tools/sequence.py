#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/modules/tools/sequence.py
## @brief  general sequence tools for the toolkit- To be replaced by region class
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


from rosetta import *
import loops

def get_sequence(pose, loop_string):
    """
    Returns the sequence of the loop string.
    Works on just chain, or N ter and C ter residues.
    """
    
    if pose.total_residue()==0:
        return ""
    residue_array = loops.return_residue_array(pose, loop_string)
    if not residue_array:return

    x = pose.sequence()
    seq = ""
    for num in residue_array:
        seq = seq+x[num-1]
    if seq=="I":
        print "Region does not exist in PDB"
        return
    if seq=="":
        seq = pose.sequence()
    return seq
