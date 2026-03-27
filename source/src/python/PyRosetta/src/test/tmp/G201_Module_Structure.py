#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   /GUIs/tests/Test_Module_Structure.py
## @brief  PyRosetta Toolkit Test script
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

from __future__ import print_function

#Python Imports
import os
import sys

#Must have Toolkit added to path!
from app.pyrosetta_toolkit.modules.Structure import Antibody_Structure

print "Testing simple Antibody Structure class- Use AntibodyInfo if using Rosetta."

ab_struct = Antibody_Structure()
cdrs = ab_struct.get_CDRs()
for cdr in cdrs:
    print str(cdr)
    
l1 = ab_struct.get_CDR("L1")
print l1.get_pdb_chain()+" "+str(l1.get_pdb_start())+" "+str(l1.get_pdb_end())
