
from rosetta import *
rosetta.init()


print 'constraints ----------------------------------------------'

pose = pose_from_pdb("test/data/test_in.pdb")

scorefxn = get_fa_scorefxn() #  create_score_function('standard')
scorefxn.set_weight( atom_pair_constraint, 10 )
scorefxn(pose)

print
print 'Score before constraints applied:', scorefxn(pose)

hterm2 = AtomID(1, 55)
gterm2 = AtomID(1, 5)

GF = constraints.GaussianFunc( 4.0, 2.0 )
apc = constraints.AtomPairConstraint( hterm2, gterm2, GF )
pose.add_constraint( apc )

print 'Score after constraints applied:', scorefxn(pose)
