# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

# test to check that core::pack::rotamer_set::RotamerSets is convertable to core::pack_basic::RotamerSetsBase in Python

## @author James Lucas, Sergey Lyskov

from __future__ import print_function

import pyrosetta
import pyrosetta.rosetta as rosetta

pyrosetta.init()
print( pyrosetta.version() )

match_pose = pyrosetta.pose_from_sequence('ASDF')
sfxn = rosetta.core.scoring.get_score_function()
rotamer_sets = rosetta.core.pack.rotamer_set.RotamerSetsFactory.create_rotamer_sets(match_pose)
sfxn.setup_for_packing_with_rotsets(match_pose, rotamer_sets)
