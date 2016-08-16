# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @author Sergey Lyskov

from __future__ import print_function

import sys
#if sys.platform == "darwin": sys.exit(0)  # skipping this test on Mac OS due to memory error

from rosetta import *
from pyrosetta import *

init(extra_options = "-constant_seed")  # WARNING: option '-constant_seed' is for testing only! MAKE SURE TO REMOVE IT IN PRODUCTION RUNS!!!!!
import os; os.chdir('.test.output')


print( 'testing ClassicRelax' )
relax_p = Pose()
core.import_pose.pose_from_file(relax_p, "../test/data/test_in.pdb")
scorefxn = get_fa_scorefxn() #  create_score_function_ws_patch('standard', 'score12')
relax = protocols.relax.ClassicRelax(scorefxn)
relax.set_lj_ramp_cycles(3)
relax.set_lj_ramp_inner_cycles(3)
relax.set_stage2_cycles(10)
relax.set_stage3_cycles(10)
relax.apply(relax_p)
