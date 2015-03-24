#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   _tempate_.py
## @brief  scientific/_tempate_.py
## Tempalte Python script for implementing Scientific test. This is template for tests that does not require extra file to run.
## For example of how to implement test which require extra files (beside access to ‘data’ checkout) please see rosetta/benchmark/tests/benchmark/tests/scientific/_template_dir_
## @author Sergey Lyskov

import os, json

import imp
imp.load_source(__name__, '/'.join(__file__.split('/')[:-1]) +  '/../__init__.py')  # A bit of Python magic here, what we trying to say is this: from ../__init__ import *, but init is calculated from file location

_api_version_ = '1.0'  # api version


def run(test, rosetta_dir, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):

    jobs = config['cpu_count']

    # Building Rosetta binaries
    res, output, build_command_line = build_rosetta(rosetta_dir, platform, jobs, debug=debug)
    if res: return { _StateKey_ : _S_build_failed_,  _ResultsKey_ : {},  _LogKey_ : 'Building rosetta failed!\n{}\n{}\n'.format(build_command_line, output) }
    else:
        extension = calculate_extension(platform)

        res, output = execute('Scoring PDB...', 'cd  {working_dir} && {rosetta_dir}/source/bin/score.{extension} -s {rosetta_dir}/source/test/core/io/test_in.pdb'.format(**vars()), return_='tuple')

        state = _S_failed_ if res else _S_finished_

        return {_StateKey_ : state,  _ResultsKey_ : {},  _LogKey_ : output }
