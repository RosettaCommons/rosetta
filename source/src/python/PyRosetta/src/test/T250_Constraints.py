from __future__ import print_function

from rosetta import *

import rosetta.core.scoring.func

rosetta.init(extra_options = "-constant_seed")  # WARNING: option '-constant_seed' is for testing only! MAKE SURE TO REMOVE IT IN PRODUCTION RUNS!!!!!
import os; os.chdir('.test.output')

print('constraints ----------------------------------------------')

pose = core.import_pose.pose_from_file("../test/data/test_in.pdb")

scorefxn = get_fa_scorefxn() #  create_score_function('standard')
scorefxn.set_weight( atom_pair_constraint, 10 )
scorefxn(pose)

print()
print( 'Score before constraints applied:', scorefxn(pose) )

hterm2 = AtomID(1, 55)
gterm2 = AtomID(1, 5)

GF = rosetta.core.scoring.func.GaussianFunc( 4.0, 2.0 )
apc = constraints.AtomPairConstraint( hterm2, gterm2, GF )
pose.add_constraint( apc )

print('Score after constraints applied:', scorefxn(pose))
