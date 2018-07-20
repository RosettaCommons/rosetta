#!/usr/bin/env python
# :noTabs=true:
# -*- coding: utf-8 -*-
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   T201_Scoring_pre_talaris.py
## @brief  Test to trigger creation of database/rotamer/bbdep02.May.sortlib.Dunbrack02.lib.bin
## @brief        binaries on windows.
## @author Sergey Lyskov

from __future__ import print_function

import pyrosetta
from pyrosetta import *

pyrosetta.init(extra_options = "-restore_pre_talaris_2013_behavior")  # trigger generation of database/rotamer/bbdep02.May.sortlib.Dunbrack02.lib.bin
pose2 = pyrosetta.pose_from_sequence("ARNDCEQGHILKMFPSTWYV", 'fa_standard')

scorefxn = rosetta.core.scoring.get_score_function()
scorefxn(pose2)
