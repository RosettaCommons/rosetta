#!/usr/bin/env python
# :noTabs=true:
# -*- coding: utf-8 -*-
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   T201_Scoring_pre_talaris.py
## @brief  Test to trigger creation of database/rotamer/bbdep02.May.sortlib.Dunbrack02.lib.bin
## @brief        binaries on windows.
## @author Sergey Lyskov

from __future__ import print_function

import rosetta, pyrosetta
from pyrosetta import *

pyrosetta.init(extra_options = "-restore_pre_talaris_2013_behavior")  # trigger generation of database/rotamer/bbdep02.May.sortlib.Dunbrack02.lib.bin
pose2 = rosetta.core.pose.Pose()
pyrosetta.make_pose_from_sequence(pose2, "ARNDCEQGHILKMFPSTWYV", 'fa_standard')

scorefxn = rosetta.core.scoring.get_score_function()
scorefxn(pose2)
