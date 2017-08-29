# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @author Sergey Lyskov

from __future__ import print_function

import pyrosetta

pyrosetta.init(extra_options = "-constant_seed")  # WARNING: option '-constant_seed' is for testing only! MAKE SURE TO REMOVE IT IN PRODUCTION RUNS!!!!!
import os; os.chdir('.test.output')

# Testing bind_dunder_method
@pyrosetta.util.bind_method(pyrosetta.rosetta.core.pose.Pose)
def __test__(self):
    pass

# Testing vector1_indices
assert(list(range(*pyrosetta.util.vector1_indices(slice(None, None, None), 10))) == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
assert(list(range(*pyrosetta.util.vector1_indices(slice(None, None, 2), 10))) == [1, 3, 5, 7, 9])
assert(list(range(*pyrosetta.util.vector1_indices(slice(None, None, 3), 10))) == [1, 4, 7, 10])
assert(list(range(*pyrosetta.util.vector1_indices(slice(5, None, None), 10))) == [5, 6, 7, 8, 9, 10])
assert(list(range(*pyrosetta.util.vector1_indices(slice(5, None, 2), 10))) == [5, 7, 9])
assert(list(range(*pyrosetta.util.vector1_indices(slice(None, 5, None), 10))) == [1, 2, 3, 4])
assert(list(range(*pyrosetta.util.vector1_indices(slice(None, 5, 2), 10))) == [1, 3])
assert(list(range(*pyrosetta.util.vector1_indices(slice(1, 11, None), 10))) == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
assert(list(range(*pyrosetta.util.vector1_indices(slice(None, -3, None), 10))) == [1, 2, 3, 4, 5, 6, 7])
