#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   build.py
## @brief  Rosetta/PyRosetta build tests
## @author Sergey Lyskov

import os, json, commands

import imp
imp.load_source(__name__, '/'.join(__file__.split('/')[:-1]) +  '/__init__.py')  # A bit of Python magic here, what we trying to say is this: from __init__ import *, but init is calculated from file location

tests = dict(
    debug   = './scons.py bin cxx={compiler} -j{jobs}',
    release = './scons.py bin cxx={compiler} mode=release -j{jobs}',
    static  = './scons.py bin cxx={compiler} mode=release extras=static -j{jobs}',
    header  = 'cd src && python ./../../../tools/python_cc_reader/test_all_headers_compile_w_fork.py -n {jobs}',
)

_TestSuite_ = False  # Set to True for TestSuite-like tests (Unit, Integration, Sfxn_fingerprint) and False other wise

def set_up():
    pass


def tear_down():
    pass


def rollover():
    pass


def get_tests():
    return tests.keys()



def run_test(test, rosetta_dir, working_dir, platform, jobs=1, hpc_driver=None, verbose=False, debug=False):
    ''' Run single test.
        Platform is a dict-like object, mandatory fields: {os='Mac', compiler='gcc'}
    '''

    TR = Tracer(verbose)

    TR('Running test: "{test}" at working_dir={working_dir!r} with rosetta_dir={rosetta_dir}, platform={platform}, jobs={jobs}, hpc_driver={hpc_driver}...'.format( **vars() ) )

    if debug: res, output = 0, 'build.py: debug is enabled, skippig build phase...\n'
    else: res, output = execute('Compiling...', 'cd {}/source && {}'.format(rosetta_dir,
                                                                            tests[test].format(compiler=platform['compiler'], jobs=jobs)),
                                return_='tuple')

    file(working_dir+'/build-log.txt', 'w').write(output)

    res_code = _S_failed_ if res else _S_finished_

    #if res: return _Failed_,   output  # We do not use '_BuildFailed_' because build failed for us actually mean test failure
    #else:   return _Finished_, output

    r = {_StateKey_ : res_code,  _ResultsKey_ : {},  _LogKey_ : output }

    # makeing sure that results could be serialize in to json, but ommiting logs because they could take too much space
    json.dump({_ResultsKey_:r[_ResultsKey_], _StateKey_:r[_StateKey_]}, file(working_dir+'/output.json', 'w'), sort_keys=True, indent=2)

    return r


def run_test_suite(rosetta_dir, working_dir, platform, jobs=1, hpc_driver=None, verbose=False, debug=False):
    TR = Tracer(verbose)
    TR('Build script does not support TestSuite-like run!')
    raise BenchmarkError()


'''
def run_test_suite(rosetta_dir, working_dir, jobs=1, hpc_driver=None, verbose=False):
    TR = Tracer(verbose)

    TR('Running test_suite: "{}" at working_dir={working_dir!r} with rosetta_dir={rosetta_dir} jobs={jobs}, hpc_driver={hpc_driver}...'.format(__name__, **vars() ) )
    results = {}
    for t in tests:
        test_working_dir = working_dir + '/' + t
        if os.path.isdir(test_working_dir): shutil.rmtree(test_working_dir);  #print('Removing old job dir %s...' % working_dir)  # remove old dir if any
        os.makedirs(test_working_dir)

        results[t] = run_test(test=t, rosetta_dir=rosetta_dir, working_dir=test_working_dir, jobs=jobs, hpc_driver=hpc_driver, verbose=verbose)

    return results
'''
