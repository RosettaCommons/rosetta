#!/usr/bin/env python
# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @author Sergey Lyskov
## @brief  Demo for PyRosetta sub classing

from __future__ import print_function

import sys

import pyrosetta
import pyrosetta.rosetta as rosetta

import os; os.chdir('.test.output')

pyrosetta.init(extra_options = "-constant_seed")  # WARNING: option '-constant_seed' is for testing only! MAKE SURE TO REMOVE IT IN PRODUCTION RUNS!!!!!

pose = pyrosetta.pose_from_file("../test/data/test_in.pdb")


some_storage = []
# rosetta.core.scoring.methods.ContextIndependentOneBodyEnergy sub-classing -----------------------------------
@pyrosetta.EnergyMethod(version=2)  # version is optional here
class MyCI1B_Method(rosetta.core.scoring.methods.ContextIndependentOneBodyEnergy):
    def __init__(self):
        rosetta.core.scoring.methods.ContextIndependentOneBodyEnergy.__init__(self, self.creator() )

    def residue_energy(self, rsd, pose, emap):
        emap.set( self.scoreType, 2.0)

    # you can define this functions by hand. And if you don't @EnergyMethod will do that for you
    #def clone(self):
    #    rosetta._mem_EnergyMethods_.append( self.__class__() )
    #    return rosetta._mem_EnergyMethods_[-1]
    #def version(self): return 141
    #def indicate_required_context_graphs(self, v): pass


sf_new = rosetta.core.scoring.ScoreFunction()
sf_new.set_weight(MyCI1B_Method.scoreType, 1)
print( '---------------------------------------------' )
print( 'MyCI1B_Method Score:', sf_new.score(pose) )
assert sf_new.score(pose) == 232.0

kT = 1.0
mc = pyrosetta.MonteCarlo(pose, sf_new, kT)
print( mc )
