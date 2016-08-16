#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   /GUIs/pyrosetta_toolkit/modules/tools/protocols.py
## @brief  basic protocol functions
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


from rosetta import *
import tkFileDialog
import tkMessageBox
import tkSimpleDialog
from Tkinter import Listbox
import os
import re



def toCen(p):
    print "Switching to Centroid"
    to_cen = SwitchResidueTypeSetMover('centroid')
    to_cen.apply(p)
    return p

def toFA(start, p):
    print "Switching to Full Atom"
    recover_sidechains = ReturnSidechainMover(start)
    #task_pack = TaskFactory.create_packer_task(start)
    #task_pack.restrict_to_repacking()
    #task_pack.or_include_current( True )
    #score = create_score_function_ws_patch('standard', 'score12')
    #pack = PackRotamersMover(score, task_pack )
    to_FA = SwitchResidueTypeSetMover('fa_standard')
    to_FA.apply(p)
    recover_sidechains.apply(p)
    #pack.apply(p)

    return p


def bbMover(p, res, chain, phi, psi, omega, score):
    """
    Moves individual Phi and Psi angles
    """

    res = p.pdb_info().pdb2pose(chain, res)
    p.set_phi(res, phi)
    print score(p)
    p.set_psi(res, psi)
    print score(p)
    p.set_omega(res, omega)
    print score(p)

    return p





