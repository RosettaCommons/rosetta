#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   _tempate_.py
## @brief  scientific/_tempate_.py
## Tempalte Python script for implementing Scientific test. This is template for tests that does not require extra file to run.
## For example of how to implement test which require extra files (beside access to ‘data’ checkout) please see rosetta/benchmark/tests/benchmark/tests/scientific/_template_dir_
## @author Sergey Lyskov

import os, json

import imp
imp.load_source(__name__, '/'.join(__file__.split('/')[:-1]) +  '/../__init__.py')  # A bit of Python magic here, what we trying to say is this: from ../__init__ import *, but init path is calculated relatively to this location


_api_version_ = '1.0'  # api version


def run(test, rosetta_dir, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):

    jobs = config['cpu_count']

    # Building Rosetta binaries
    res, output, build_command_line = build_rosetta(rosetta_dir, platform, config)
    if res: return { _StateKey_ : _S_build_failed_,  _ResultsKey_ : {},
                     _LogKey_ : 'Building rosetta failed!\n{}\n{}\n'.format(build_command_line, output) }
    else:
        extension = calculate_extension(platform)

        # running Rosetta executable on benchmark machine itself
        res, output = execute('Scoring PDB...',
                              'cd  {working_dir} && {rosetta_dir}/source/bin/score.{extension} -s {rosetta_dir}/source/test/core/io/test_in.pdb'.format(**vars()), return_='tuple')


        # running Rosetta executable on benchmark machine in parallel and using multiple CPU's (=jobs)
        scoring_jobs = dict(scoring1='cd {working_dir} && {rosetta_dir}/source/bin/score.{extension} -s {rosetta_dir}/source/test/core/io/test_in.pdb -out:file:scorefile parallel.score.1.sc'.format(**vars()),
                            scoring2='cd {working_dir} && {rosetta_dir}/source/bin/score.{extension} -s {rosetta_dir}/source/test/core/io/test_in.pdb -out:file:scorefile parallel.score.2.sc'.format(**vars()),
                            scoring3='cd {working_dir} && {rosetta_dir}/source/bin/score.{extension} -s {rosetta_dir}/source/test/core/io/test_in.pdb -out:file:scorefile parallel.score.3.sc'.format(**vars()),
                            scoring4='cd {working_dir} && {rosetta_dir}/source/bin/score.{extension} -s {rosetta_dir}/source/test/core/io/test_in.pdb -out:file:scorefile parallel.score.4.sc'.format(**vars()), )

        parallel_execute('parallel_scoring', scoring_jobs, rosetta_dir, working_dir, cpu_count=jobs, time=1)


        # running executable on HPC cluster
        hpc_driver.execute(executable='{rosetta_dir}/source/bin/score.{extension}'.format(**vars()),
                           arguments='-s {rosetta_dir}/source/test/core/io/test_in.pdb -out:file:scorefile hpc_score.sc'.format(**vars()),
                           working_dir=working_dir, name='scoring')

        # insert result analyzing code here
        if not ( os.path.isfile(working_dir+'/default.sc') and
                 os.path.isfile(working_dir+'/hpc_score.sc') ): res = 1 # No score files found - something is wrong, terminating with error!!!

        state = _S_failed_ if res else _S_passed_

        # Simple form when our test does not have any sub-tests:
        # return {_StateKey_ : state,  _ResultsKey_ : {},  _LogKey_ : output }

        # More complex form when we have a sub-tests which we want Testing Server to track
        # _SummaryKey_ is dict to provide optional information specific to this test
        results = { _SummaryKey_ : { 'arbitrary_key_1' : 'arbitrary_info_1' },

                    _TestsKey_ : {'sub-test-1' : dict(state=_S_passed_, log='sub-test-1 log...'),
                                  'sub-test-2' : dict(state=_S_failed_, log='sub-test-2 log...') } }

        return {_StateKey_ : state,  _ResultsKey_ : results,  _LogKey_ : output }
