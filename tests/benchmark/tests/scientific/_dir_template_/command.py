#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   _tempate_.py
## @brief  scientific/_template_dir_/command.py
## Tempalte Python test script showing how to implement test located in separate dir
## For example of how to implement test which does not require extra files (beside access to 'data' checkout) please see rosetta/benchmark/tests/benchmark/tests/scientific/_template_.py
## @author Sergey Lyskov

import os, json

import imp
imp.load_source(__name__, '/'.join(__file__.split('/')[:-1]) +  '/../../__init__.py')  # A bit of Python magic here, what we trying to say is this: from ../../__init__ import *, but init is calculated from file location

_api_version_ = '1.0'  # api version

def run(test, rosetta_dir, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    return {_StateKey_ : _S_failed_,  _ResultsKey_ : {},  _LogKey_ : 'Put test log here...'}
