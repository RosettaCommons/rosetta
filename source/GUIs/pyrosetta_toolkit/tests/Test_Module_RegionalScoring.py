#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/tests/Test_Modules
## @brief  PyRosetta Toolkit Test script
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#Python Imports
import os
import sys
rel = "../"
sys.path.append(os.path.abspath(rel))

#Toolkit Imports
from modules.RegionalScoring import RegionalScoring

#Rosetta Imports
from rosetta import *
rosetta.init()

p = pose_from_pdb(os.path.dirname(os.path.abspath(__file__))+"/data/2j88.pdb")
scorefxn = create_score_function("talaris2013")


print "Testing RegionalScoring Class \n\n"
reg_score = RegionalScoring(p, scorefxn)
score = reg_score.ret_energy_string()
print score

print "Weighted energy of residue 30: "+repr(reg_score.ret_total_weighted_residue_energy(30))
print "Unweighted energy of residue 30: "+repr(reg_score.ret_total_unweighted_residue_energy(30))