#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   /GUIs/tests/Test_Analysis_Functions
## @brief  PyRosetta Toolkit Test script
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

from __future__ import print_function

#Python Imports
import os
import sys


#Must have Toolkit added to path!
#Toolkit Imports
from app.pyrosetta_toolkit.modules.tools import analysis

#Rosetta Imports
from rosetta import *
rosetta.init()


p = pose_from_file(os.path.dirname(os.path.abspath(__file__))+"/data/gui/2j88.pdb")

print "Printing all Phi Psi:"
analysis.print_all_phi_psi(p)

print "Printing all energies:"
analysis.print_all_residue_energies(p)

print "Running PackStat:"
analysis.analyze_packing(p)

print "Energy of FaDun Residue 5:"
print repr(analysis.return_energy(p, 5))

print "Rotamer Prob of Residue 5:"
print repr(analysis.return_probability(p, 5))

print "Running LoopAnalyzeMover:"
loop_string = "24:42:L"
loops_as_strings = []
loops_as_strings.append(loop_string)
analysis.analyze_loops(p, loops_as_strings)

print "!! Cannot test interface analyzer"
print "!! Not testing VIP"

print "\n Complete \n"


