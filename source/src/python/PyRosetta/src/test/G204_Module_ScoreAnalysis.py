#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/tests/Test_Module_Region.py
## @brief  PyRosetta Toolkit Test script
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#Python Imports
import os
import sys

#Must have Toolkit added to path!
#Toolkit Imports
from app.pyrosetta_toolkit.modules.ScoreAnalysis import ScoreAnalysis


score_analysis = ScoreAnalysis()
score_analysis.read_scores(os.path.dirname(os.path.abspath(__file__))+"/data/gui/SCORED_PDBLIST.txt")


print "Getting top score model:"
top = score_analysis.get_top_scoring(False)
print top

print "Getting top score by percent:"
tops = score_analysis.get_top_scoring_by_percent(20.0, False)

print "Getting top score by number"
tops = score_analysis.get_top_scoring_by_number(3, False)

#No testing of RMSD as I don't have 10 PDBs to throw in trunk.  Nor should I.  Figure out another way to test this.
