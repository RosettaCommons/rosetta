#!/usr/bin/env python
# :noTabs=true:
# -*- coding: utf-8 -*-

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   T040_Types.py
## @brief  Tests for bindings for various types
## @author Sergey Lyskov

from __future__ import print_function

import pyrosetta
import pyrosetta.rosetta as rosetta
from pyrosetta import *

pyrosetta.init(extra_options = "-constant_seed")  # WARNING: option '-constant_seed' is for testing only! MAKE SURE TO REMOVE IT IN PRODUCTION RUNS!!!!!
import os; os.chdir('.test.output')


cpp_set = pyrosetta.rosetta.std.set_int_t()
py_set = set()

values = [1, 1, 3, 4, 5, 6, 7, 1, 3, 4, 5, 6, 7, ]

for v in values: cpp_set.add(v)

assert {v for v in cpp_set} == set(values)
