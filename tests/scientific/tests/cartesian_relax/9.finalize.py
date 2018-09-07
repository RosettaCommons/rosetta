#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  cartesian_relax/9.finalize.py
## @brief this script is part of cartesian_relax scientific test
## @author Sergey Lyskov


import json
import benchmark

from benchmark import *

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals

with open(_multi_step_result_, 'w') as f:
    r = {
        _StateKey_  : _S_passed_,
        _ResultsKey_ : {},
        _LogKey_ : 'Done!',
    }

    json.dump(r, f, sort_keys=True, indent=2)


benchmark.save_variables()  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
