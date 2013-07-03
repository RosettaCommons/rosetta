#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/tests/Test_Sequence_Functions
## @brief  PyRosetta Toolkit Test script
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#Python Imports
import os
import sys
rel = "../"
sys.path.append(os.path.abspath(rel))

#Toolkit Imports
from modules.tools import sequence

from rosetta import *
rosetta.init()

p = pose_from_pdb(os.path.dirname(os.path.abspath(__file__))+"/data/2j88.pdb")

print "Testing alternative way to get sequence:"
print sequence.get_sequence(p, "24:42:L")