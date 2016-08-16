#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   /GUIs/tests/Test_Output_Functions
## @brief  PyRosetta Toolkit Test script
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

from __future__ import print_function

#Python Imports
import os
import sys

                
#Must have Toolkit added to path!
#Toolkit Imports
from app.pyrosetta_toolkit.modules.tools import output
from app.pyrosetta_toolkit.modules.tools import general_tools

#Rosetta Imports
from rosetta import *
rosetta.init()


p = pose_from_file(os.path.dirname(os.path.abspath(__file__))+"/data/gui/2j88.pdb")

os.chdir(".test.output/")
scorefxn = create_score_function("talaris2013")

loop_string = "24:42:L"
loops_as_strings = []; loops_as_strings.append(loop_string); loops_as_strings.append("24:42:H")
regions = general_tools.loops_as_strings_to_regions(loops_as_strings)

print "Testing main PDB Dump:"
output.dumpPDB(p, p, "gui_output_pdb_test.pdb", scorefxn, True)

print "Testing loop output:"
output.save_loop_file(p, regions, False, "gui_output_loop_test.loop")

print "Testing SCWRL seq output:"
output.saveSeqFile(p, "gui_output_seq_test.txt", loops_as_strings)

print "Testing basic Resfile output:"
output.save_basic_resfile(p, "gui_output_basic_resfile_test.resfile")

print "Testing basic Blueprint output:"
output.save_basic_blueprint(p, True, "gui_output_basic_blueprint_test.blueprint")

#Need to test with a designdic.  Actually, need to refactor the designdic.  Its horrible.

print "Testing full pose fasta output:"
output.save_FASTA(p, "2J88", "gui_output_full_pose_fasta_test.fasta")

print "Testing region based fasta output:"
output.save_FASTA(p, "2J88", "gui_output_regions_fasta_test.fasta", regions)

#make_pdblist

#make_pdblist_recursively

#convert pdblist to sqlite3

#extract pdb from sqlite3

#rescore single pdb

#score PDBList

#Save FASTPdbList

#Export PDBScore

#Output molfile to params

#Save param path list?

                          
