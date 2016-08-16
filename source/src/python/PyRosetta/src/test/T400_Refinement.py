# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @author Sergey Lyskov

from __future__ import print_function

import sys
#if sys.platform == "darwin": sys.exit(0)  # skipping this test on Mac OS due to memory error (mover: ClassicAbinitio ERROR)

from rosetta import *
from pyrosetta import *

init(extra_options = "-constant_seed")  # WARNING: option '-constant_seed' is for testing only! MAKE SURE TO REMOVE IT IN PRODUCTION RUNS!!!!!
import os; os.chdir('.test.output')

print('Refinement ----------------------------------------------')

kT=1.0
n_moves=10

pose = core.import_pose.pose_from_file("../test/data/test_fragments.pdb")
pose_frag = core.import_pose.pose_from_file("../test/data/test_fragments.pdb")


print('setting up a move map')
movemap = MoveMap()
print('setting all backbone movement to true, and residue 10 to false')
movemap.set_bb(True)
movemap.set_bb(10,False)
print('outputting movemap')
movemap.show(pose.total_residue())


fragset3mer = core.fragment.ConstantLengthFragSet(3, "../test/data/test3_fragments")# "aatestA03_05.200_v1_3")
fragset9mer = core.fragment.ConstantLengthFragSet(9, "../test/data/test9_fragments")# "aatestA09_05.200_v1_3")
print('mover: ClassicFragmentMover, 3mer')
movemap.set_bb(True)
mover_3mer = protocols.simple_moves.ClassicFragmentMover(fragset3mer, movemap)
mover_3mer.apply(pose_frag)


print('Creating standard score function with patch and scoring')
scorefxn = get_score_function()
scorefxn(pose)


print('mover: SmallMover')
movemap.set_bb(True)
smallmover = protocols.simple_moves.SmallMover(movemap,kT,n_moves)
smallmover = protocols.simple_moves.SmallMover()
smallmover.angle_max('L',50)
smallmover.apply(pose_frag)
print(smallmover)
# TODO: fix printout.  show mover params (temp, nmoves, anglemax) and last output status

# TODO: mover base class printout with name, last output status, brief description(?)

print('mover: ShearMover')
shearmover = protocols.simple_moves.ShearMover(movemap, kT, n_moves)
shearmover = protocols.simple_moves.ShearMover()
shearmover.angle_max('E',50)
shearmover.apply(pose_frag)
print(shearmover)
#TODO: same as above

print('mover: MinMover')
minmover = protocols.simple_moves.MinMover()
minmover.movemap(movemap)
minmover.score_function(scorefxn)
minmover.apply(pose_frag)

# currently crash!
minmover = protocols.simple_moves.MinMover(movemap,scorefxn,'linmin',0.5, True) #nblist mover was recent undefaulted
# the above crashes unless False changes to True -- what is this about?
# TODO: prefered solution is to overload the last argument to a default, valid value
minmover.apply(pose_frag)
minmover.min_type('dfpmin')  # set_mintype is no longer in C++!
minmover.apply(pose_frag)

#NEEDED: AbInitio [ERROR - OP], Relax [no bindings], Idealize[no bindings]

print('mover: MonteCarlo')
mc = MonteCarlo(pose, scorefxn, kT)
mc.boltzmann(pose)

mc.reset(pose)
mc.boltzmann(pose)
smallmover.apply(pose)
mc.boltzmann(pose)
smallmover.apply(pose)
mc.boltzmann(pose)
print(mc) # TODO: this should tell the pose, scorefunction and temperature
mc.show_scores()
mc.show_counters()
mc.show_state()

#print 'mover: ClassicAbinitio ERROR'
#''' Disabled until Abinitio error is fixed '''
ab = protocols.abinitio.ClassicAbinitio(fragset3mer, fragset9mer, movemap)
#ab.set_mc( MonteCarlo(pose, scorefxn, kT) )
ab.init(pose)  # needed to be se or ClassicAbinitio will crash
#ab.apply(pose)
# ^^^ uncommenting lead to an error:
#protocols.simple_moves.FragmentMover: BEGIN: 89 SIZE: 9 TOTAL_RES: 10
#protocols.simple_moves.FragmentMover: Are the fragments compatible with the fasta or the input PDB used to extract the folding sequence ?
#protocols.simple_moves.FragmentMover: It appears that the fragments go up to residue 97 while the pose only has 10 residues!


print('mover: SequenceMover')
seqmover = protocols.moves.SequenceMover()
seqmover.add_mover(minmover)
seqmover.apply(pose)

print('mover: TrialMover')
trialmover = TrialMover(mover_3mer, mc)
trialmover.apply(pose)
trialmover.num_accepts()
trialmover.acceptance_rate()
mc.show_state() # this should have breakdown by mover type now
