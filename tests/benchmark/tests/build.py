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

import os, commands

class BenchmarkBuildError(Exception): pass


tests = dict(
    debug   = './scons.py bin -j{jobs}',
    release = './scons.py bin mode=release -j{jobs}',
)


def set_up():
    pass


def tear_down():
    pass


def rollover():
    pass


def get_tests():
    return tests.keys()


def run_test(test, rosetta_dir, working_dir, jobs=1, hpc_driver=None, verbose=False):
    TR = Tracer(verbose)

    TR('Running test: "{test}" at working_dir={working_dir!r} with rosetta_dir={rosetta_dir} jobs={jobs}, hpc_driver={hpc_driver}...'.format( **vars() ) )

    res, output = execute('Compiling...', 'cd {}/source && {}'.format(rosetta_dir, tests[test].format(jobs=jobs)), return_='tuple')

    file(working_dir+'/build-log.txt', 'w').write(output)

    res_code = _Failed_ if res else _Finished_

    #if res: return _Failed_,   output  # We do not use '_BuildFailed_' because build failed for us actually mean test failure
    #else:   return _Finished_, output

    return {_StateKey_ : res_code,  _ResultsKey_ : {},  _LogKey_ : output }



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



# Standard funtions and constants below ---------------------------------------------------------------------------------
# Do not change this wording, they have to stay in sync with upstream (up to benchmark-model).
_Finished_     = '_Finished_'
_Failed_       = '_Failed_'
_BuildFailed_  = '_BuildFailed_'
_ScriptFailed_ = '_ScriptFailed_'

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
