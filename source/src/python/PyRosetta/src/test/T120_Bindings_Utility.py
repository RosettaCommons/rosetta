# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

from __future__ import print_function

import pyrosetta

pyrosetta.init(extra_options = "-constant_seed")  # WARNING: option '-constant_seed' is for testing only! MAKE SURE TO REMOVE IT IN PRODUCTION RUNS!!!!!
import os; os.chdir('.test.output')


# Testing bind_method
@pyrosetta.bindings.utility.bind_method(pyrosetta.rosetta.core.pose.Pose)
def test_method(self):
    return True
pose = pyrosetta.rosetta.core.pose.Pose()
assert pose.test_method() == True, 'bind_method failed'


# Testing bind_property
@pyrosetta.bindings.utility.bind_property(pyrosetta.rosetta.core.pose.Pose)
def test_property(self):
    return True
pose = pyrosetta.rosetta.core.pose.Pose()
assert pose.test_property == True, 'bind_property failed'


# Testing vector1_indices
from pyrosetta.bindings.utility import slice_1base_indicies
try:
    slice_1base_indicies(slice(0, None, None), 10)
    assert False, "Did not raise ValueError"
except ValueError:
    pass

assert slice_1base_indicies(slice(None, None, None), 10) == (1, 11, 1)
assert slice_1base_indicies(slice(None, None, 2), 10) == (1, 11, 2)
assert slice_1base_indicies(slice(5, None, None), 10) == (5, 11, 1)
assert slice_1base_indicies(slice(5, None, 2), 10) == (5, 11, 2)
assert slice_1base_indicies(slice(None, 10, None), 10) == (1, 10, 1)
assert slice_1base_indicies(slice(None, 10, 2), 10) == (1, 10, 2)
assert slice_1base_indicies(slice(None, 5, None), 10) == (1, 5, 1)
assert slice_1base_indicies(slice(None, 5, 2), 10) == (1, 5, 2)
assert slice_1base_indicies(slice(1, 11, 1), 10) == (1, 11, 1)
assert slice_1base_indicies(slice(None, -3, None), 10) == (1, 8, 1)
