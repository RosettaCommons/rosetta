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
rel = "../"
sys.path.append(os.path.abspath(rel))

#Toolkit Imports
from modules.tools import input

#Rosetta Imports
from rosetta import *
rosetta.init()

p = pose_from_pdb(os.path.dirname(os.path.abspath(__file__))+"/data/2j88.pdb")


print "Testing Directory Search."

print "Testing Adding constraints to pose and scorefunction."

print "Testing loading a ResidueTypeSet from path array."

print "Complete."
