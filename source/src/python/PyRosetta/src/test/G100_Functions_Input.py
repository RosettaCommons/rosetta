#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/tests/Test_Input_Functions
## @brief  PyRosetta Toolkit Test script
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#Python Imports
import os
import sys


#Must have Toolkit added to path!
#Toolkit Imports
from app.pyrosetta_toolkit.modules.tools import input as input_tools

#Rosetta Imports
from rosetta import *
rosetta.init()

p = pose_from_file(os.path.dirname(os.path.abspath(__file__))+"/data/gui/2j88.pdb")
score = get_score_function()


#print "Testing Directory Search."

print "Testing Adding constraints to pose and scorefunction."

cst_file = os.path.dirname(os.path.abspath(__file__))+"/data/gui/L1-11-1.txt"
input_tools.add_constraints_to_pose_and_scorefunction(p, score, 1.0, [dihedral_constraint], cst_file)

#print "Testing loading a ResidueTypeSet from path array."

#print "Complete."
