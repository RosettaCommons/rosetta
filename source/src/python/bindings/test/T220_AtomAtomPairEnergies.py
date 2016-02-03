import sys


import rosetta
from rosetta import *

init(extra_options = "-constant_seed")  # WARNING: option '-constant_seed' is for testing only! MAKE SURE TO REMOVE IT IN PRODUCTION RUNS!!!!!
import os; os.chdir('.test.output')

test_pose = pose_from_file('../test/data/test_dock.pdb')
residue_1 = test_pose.residue(275)
residue_2 = test_pose.residue(55)

#sfxn =create_score_function('score12')
sfxn = create_score_function('talaris2013')

#calculating atom-atom pairwise interactions and summing
#to get total energy (which should match, residue-residue energy!!!)
atr_total, rep_total, solv_total = 0.0, 0.0, 0.0


for i in range(residue_1.natoms()):
	for j in range(residue_2.natoms()):
		atom1 = residue_1.atom(i+1)
		atom2 = residue_2.atom(j+1)

		atr, rep ,solv = etable_atom_pair_energies(atom1, atom2, sfxn)

		atr_total  += atr
		rep_total  += rep
		solv_total += solv


emap = EMapVector()
sfxn.eval_ci_2b(residue_1, residue_2, test_pose,emap)

print '\n\n'
print 'Printing individual energies:'

print 'res-res atr score emap:      ', emap[fa_atr]
print 'res-res atr score pairwise:  ', atr_total
print 'res-res rep score emap:      ', emap[fa_rep]
print 'res-res rep score pairwise:  ', rep_total
print 'res-res solv score emap:     ', emap[fa_sol]
print 'res-res solv score pairwise: ', solv_total

print 'Checking if scores match...'

if ( abs(emap[fa_atr] - atr_total) +
     abs(emap[fa_rep] - rep_total) +
     abs(emap[fa_sol] - solv_total) ) > 1.0e-10:
    print 'Score did not match, exiting!!!'
    sys.exit(1)
else:
    print 'Score did match, yay...'
