#!/usr/bin/env python
# :noTabs=true:
# -*- coding: utf-8 -*-

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   T150_core.misc.lkball.py
## @brief  Tests for bindings in core::scoring::lkball namespaces
## @author Jack Maguire and Sergey Lyskov

from __future__ import print_function

import pyrosetta
import pyrosetta.rosetta as rosetta

import pyrosetta
from pyrosetta import *
from pyrosetta.rosetta import core

pyrosetta.init()

sfxn = get_score_function()

pose = pose_from_sequence("NNN")
sfxn.score( pose )

restype = pose.residue(2).type()

singleton = core.scoring.lkball.LKBallDatabase.get_instance()

water_builder_4_restype = singleton.get_water_builder_for_restype(restype)

water_builders_list = water_builder_4_restype.builders()

for i in range(1,len(water_builders_list)+1):
    print( i, len(water_builders_list[i]) )
