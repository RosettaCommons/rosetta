# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @author Sergey Lyskov

from rosetta import *

import rosetta.protocols.loops.loop_mover.refine

rosetta.init(extra_options = "-constant_seed")  # WARNING: option '-constant_seed' is for testing only! MAKE SURE TO REMOVE IT IN PRODUCTION RUNS!!!!!
import os; os.chdir('.test.output')

print 'Loop building --------------------------------------------------'

#loop functions
loop_p =  core.import_pose.pose_from_file("../test/data/test_in.pdb")

loop = Loop(70,80,75)
loops = Loops()
loops.add_loop(loop)

set_single_loop_fold_tree(loop_p, loop)

print loop_p.fold_tree()

#recover_sidechain = protocols.simple_moves.ReturnSidechainMover(loop_p)
#to_centroid.apply(loop_p)

#fragset3mer = ConstantLengthFragSet(3, "test_in3_fragments")
#scorefxn = create_score_function_ws_patch('cen_std', 'score4L')
#scorefxn(loop_p)
#loop_perturb = LoopMover_Perturb_CCD(loops, scorefxn, fragset3mer)
#loop_perturb.apply(p)
#loop_perturb.model_loop(loop_p, loop)

#recover_sidechain.apply(loop_p)

scorefxn = get_fa_scorefxn() #  create_score_function_ws_patch('standard', 'score12')
scorefxn(loop_p)
loop_refine = rosetta.protocols.loops.loop_mover.refine.LoopMover_Refine_CCD( loops, scorefxn )
loop_refine.max_inner_cycles(10)
loop_refine.apply(loop_p)
