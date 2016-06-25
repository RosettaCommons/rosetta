#!/usr/bin/env python
# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @author Sergey Lyskov
## @brief  Demo for PyRosetta sub classing 2

from __future__ import print_function

# RosettaCon 2011 demo
import pyrosetta
from rosetta import *
from pyrosetta import *
from rosetta.core.scoring.methods import *

init(extra_options = "-constant_seed")  # WARNING: option '-constant_seed' is for testing only! MAKE SURE TO REMOVE IT IN PRODUCTION RUNS!!!!!
import os; os.chdir('.test.output')

pose = pose_from_file("../test/data/test_in.pdb")


@pyrosetta.EnergyMethod()
class CI1B_Method(ContextIndependentOneBodyEnergy):
    def __init__(self):
        ContextIndependentOneBodyEnergy.__init__(self, self.creator() )

    def residue_energy(self, rsd, pose, emap):
        emap.set( self.scoreType, 2.0)

    # you can define this functions by hand. And if you don't @EnergyMethod will do that for you
    #def clone(self): return MyNewCI1B();
    #def version(self): return 141
    #def indicate_required_context_graphs(self, v): pass


sf_new = ScoreFunction()
sf_new.set_weight(CI1B_Method.scoreType, 1)
print( '---------------------------------------------' )
print( 'CI1B_Method Score:', sf_new.score(pose) )
assert sf_new.score(pose) == 232.0


@pyrosetta.EnergyMethod()
class CI2B_Method(ContextIndependentTwoBodyEnergy):
    def __init__(self):
        ContextIndependentTwoBodyEnergy.__init__(self, self.creator() )

    def residue_pair_energy(self, rsd1, rsd2, pose, sfxn, emap):
        emap.set( self.scoreType, 1.0)

    def atomic_interaction_cutoff(self): return 0.0

    def defines_intrares_energy(self, weights): return True;

    def eval_intrares_energy(self, rsd, pose, sfxn, emap): pass


sf_new = ScoreFunction()
sf_new.set_weight(CI2B_Method.scoreType, 1)
print( '---------------------------------------------' )
print( 'CI2B_Method Score:', sf_new.score(pose) )
assert sf_new.score(pose) == 525.0


@pyrosetta.EnergyMethod()
class CD2B_Method(ContextDependentTwoBodyEnergy):
    def __init__(self):
        ContextDependentTwoBodyEnergy.__init__(self, self.creator() )

    def residue_pair_energy(self, rsd1, rsd2, pose, sfxn, emap):
        emap.set( self.scoreType, 1.0)

    def atomic_interaction_cutoff(self): return 0.0

    def defines_intrares_energy(self, weights): return True;

    def eval_intrares_energy(self, rsd, pose, sfxn, emap): pass


sf_new = ScoreFunction()
sf_new.set_weight(CD2B_Method.scoreType, 1)
print( '---------------------------------------------' )
print( 'CD2B_Method Score:', sf_new.score(pose) )
assert sf_new.score(pose) == 525.0

print('\n\n')
