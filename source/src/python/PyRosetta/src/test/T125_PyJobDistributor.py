# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @author Sergey Lyskov

from __future__ import print_function

import glob

import pyrosetta

pyrosetta.init(extra_options = "-constant_seed")  # WARNING: option '-constant_seed' is for testing only! MAKE SURE TO REMOVE IT IN PRODUCTION RUNS!!!!!
import os; os.chdir('.test.output')

sf = pyrosetta.create_score_function('ref2015')
jd = pyrosetta.PyJobDistributor("_jd_test", 10, sf)

while not jd.job_complete:
    pp = pyrosetta.pose_from_sequence("DSEEKFLRRIGRFGYGYGPYE")
    jd.output_decoy(pp)

assert len( glob.glob('_jd_test*') ) == 10+1  # 10 decoys + score file


sf = pyrosetta.create_score_function('ref2015')
jd = pyrosetta.PyJobDistributor("_jd_at_test", 10, sf, compress=True)

while not jd.job_complete:
    pp = pyrosetta.pose_from_sequence("DSEEKFLRRIGRFGYGYGPYE")
    jd.output_decoy(pp)

assert len( glob.glob('_jd_at_test*') ) == 10+1  # 10 decoys + score file
