#!/usr/bin/env python
# :noTabs=true:
# -*- coding: utf-8 -*-

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   LoadPDB.py
## @brief  Serve as example on how to import rosetta, load pdb and as script that create database
## @brief        binaries on windows.
## @author Sergey Lyskov

from rosetta import *

rosetta.init()
print version()

pose = pose_from_pdb("test/data/test_in.pdb")

scorefxn = create_score_function('standard')
scorefxn(pose)


pose2 = Pose()
make_pose_from_sequence(pose2, "ARNDCEQGHILKMFPSTWYV", 'fa_standard')

scorefxn = create_score_function_ws_patch("standard", "score12")
scorefxn(pose2)


pose3 = Pose()
make_pose_from_sequence(pose3, "DSEEKFLRRIGRFGYGYGPYE",'centroid')

# Creating standard centroid score function and scoring
scorefxn = create_score_function('score3')
scorefxn(pose3)


