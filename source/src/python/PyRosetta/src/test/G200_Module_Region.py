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
rel = "../"
sys.path.append(os.path.abspath(rel))

#Toolkit Imports
from app.pyrosetta_toolkit.modules.Region import Region
from app.pyrosetta_toolkit.modules.Region import Regions

#Rosetta Imports
from rosetta import *
rosetta.init()

p = pose_from_file(os.path.dirname(os.path.abspath(__file__))+"/data/gui/2j88.pdb")


print "Testing Region Class\n\n"
reg = Region("L", 24, 42)
print str(reg)

print " Testing Basic Region functions"
assert (reg.get_chain() == "L")
assert (reg.get_start() == 24)
assert (reg.get_end() == 42)
assert reg.get_length(p) == 11
assert reg.get_region_string() == "24:42:L"
assert reg.get_region_type() == "loop"

print " Testing Region movemap functions"
mm = reg.get_movemap(p)
mm.show()

print " Testing Region loop functions"
reg.set_Loop_for_region(p, 30)
loop = reg.get_loop()
repr(loop)

print " Testing Region rosetta num and sequence functions"
print "Rosetta num: start "+repr(reg.get_rosetta_start(p))+" end "+repr(reg.get_rosetta_end(p))
print reg.get_sequence(p)


print "Testing Regions Class\n\n"
regs = Regions()
regs.add_region(reg)
region2 = Region("H", 24, 42)
region2.set_Loop_for_region(p, 133)
regs.add_region(region2)
print str(regs)
for reg in regs:
    print str(reg)

print "Testing Regions tf"
tf = regs.get_basic_tf(p)
task = tf.create_task_and_apply_taskoperations(p)
task.show()

print "Testing Regions loops"
loops = regs.get_Loops_from_regions()
mm = regs.get_movemap(p)
mm.show()

print "Testing Regions packer"
task = regs.get_packer_task(p)
task.show()

