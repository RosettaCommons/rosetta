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

class BenchmarkBuildError(Exception): pass


tests = dict(
    debug   = './scons.py bin cxx={compiler} -j{jobs}',
    release = './scons.py bin cxx={compiler} mode=release -j{jobs}',
    static  = './scons.py bin cxx={compiler} mode=release extras=static -j{jobs}',
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



def run_test(test, rosetta_dir, working_dir, platform, jobs=1, hpc_driver=None, verbose=False):
    ''' Run single test.
        Platform is a dict-like object, mandatory fields: {os='Mac', compiler='gcc'}
    '''

    TR = Tracer(verbose)

    TR('Running test: "{test}" at working_dir={working_dir!r} with rosetta_dir={rosetta_dir}, platform={platform}, jobs={jobs}, hpc_driver={hpc_driver}...'.format( **vars() ) )

    res, output = execute('Compiling...', 'cd {}/source && {}'.format(rosetta_dir, tests[test].format(compiler=platform['compiler'], jobs=jobs)), return_='tuple')

    file(working_dir+'/build-log.txt', 'w').write(output)

    res_code = _S_Failed_ if res else _S_Finished_

    #if res: return _Failed_,   output  # We do not use '_BuildFailed_' because build failed for us actually mean test failure
    #else:   return _Finished_, output

    r = {_StateKey_ : res_code,  _ResultsKey_ : {},  _LogKey_ : output }

    # makeing sure that results could be serialize in to json, but ommiting logs because they could take too much space
    json.dump({_ResultsKey_:r[_ResultsKey_], _StateKey_:r[_StateKey_]}, file(working_dir+'/output.json', 'w'), sort_keys=True, indent=2)

    return r


def run_test_suite(rosetta_dir, working_dir, platform, jobs=1, hpc_driver=None, verbose=False):
    TR = Tracer(verbose)
    TR('Build script does not support TestSuite-like run!')
    raise BenchmarkBuildError()


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


# Standard funtions and constants below ---------------------------------------------------------------------------------
# Do not change this wording, they have to stay in sync with upstream (up to benchmark-model).
# Copied from benchmark-model, standard state code's for tests results.
_S_Draft_         = 'draft'
_S_Queued_        = 'queued'
_S_Running_       = 'running'
_S_Finished_      = 'finished'
_S_Failed_        = 'failed'
_S_BuildFailed_   = 'build failed'
_S_ScriptFailed_  = 'script failed'

_S_Values_ = [_S_Draft_, _S_Queued_, _S_Running_, _S_Finished_, _S_Failed_, _S_BuildFailed_, _S_ScriptFailed_]

_StateKey_    = 'state'
_ResultsKey_  = 'results'
_LogKey_      = 'log'
# Keys below will be used only in up-stream, do to set them here
#_StartedKey_  = 'started'
#_FinishedKey_ = 'finished'


def Tracer(verbose=False):
    def print_(x): print x
    return print_ if verbose else lambda x: None



def execute(message, commandline, return_=False, untilSuccesses=False):
    TR = Tracer()
    TR(message);  TR(commandline)
    while True:
        (res, output) = commands.getstatusoutput(commandline)
        TR(output)

        if res and untilSuccesses: pass  # Thats right - redability COUNT!
        else: break

        print "Error while executing %s: %s\n" % (message, output)
        print "Sleeping 60s... then I will retry..."
        time.sleep(60)

    if return_ == 'tuple': return(res, output)

    if res:
        TR("\nEncounter error while executing: " + commandline )
        if return_==True: return True
        else: raise BenchmarkBuildError()

    if return_ == 'output': return output
    else: return False
