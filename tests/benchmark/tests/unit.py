#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   unit.py
## @brief  Rosetta/PyRosetta unit tests
## @author Sergey Lyskov

import os, json, commands

class BenchmarkIntegrationError(Exception): pass


#tests = ['i']

_TestSuite_ = True  # Set to True for TestSuite-like tests (Unit, Integration, Sfxn_fingerprint) and False other wise


def set_up():  pass


def tear_down(): pass


def get_tests():
    TR = Tracer(verbose=True)
    TR('Unit Test script does not support get_tests! Use run_test_suite instead!')
    raise BenchmarkIntegrationError()


def run_test(test, rosetta_dir, working_dir, platform, jobs=1, hpc_driver=None, verbose=False):
    TR = Tracer(verbose)
    TR('Unit Test script does not support run_test! Use run_test_suite instead!')
    raise BenchmarkIntegrationError()


def run_test_suite(rosetta_dir, working_dir, platform, jobs=1, hpc_driver=None, verbose=False):
    ''' Run TestSuite.
        Platform is a dict-like object, mandatory fields: {os='Mac', compiler='gcc'}
    '''

    TR = Tracer(verbose)
    full_log = ''

    TR('Running test_suite: "{}" at working_dir={working_dir!r} with rosetta_dir={rosetta_dir}, platform={platform}, jobs={jobs}, hpc_driver={hpc_driver}...'.format(__name__, **vars() ) )

    results = {}

    TR('Compiling...')
    res, output = execute('Compiling...', 'cd {}/source && ./scons.py -j{jobs} && ./scons.py cat=test -j{jobs}'.format(rosetta_dir, jobs=jobs), return_='tuple')
    #res, output = 0, 'debug... compiling...\n'

    full_log += output  #file(working_dir+'/build-log.txt', 'w').write(output)

    if res:
        results[_StateKey_] = _S_BuildFailed_
        results[_LogKey_]   = full_log
        return results

    else:
        json_results_file = rosetta_dir+'/source/.unit_test_results.yaml'
        if os.path.isfile(json_results_file): os.remove(json_results_file)

        command_line = 'cd {}/source && test/run.py -j{jobs}'.format(rosetta_dir, jobs=jobs)
        TR( 'Running unit test script: {}'.format(command_line) )
        res, output = execute('Running unit test script...', command_line, return_='tuple')
        full_log += output

        if res:
            results[_StateKey_] = _S_ScriptFailed_
            results[_LogKey_]   = output  # ommiting compilation log and only including integration.py output
            return results

    json_results = json.load( file(json_results_file) )
    r = {}

    for lib in json_results:
        key = lib[:-5]  # core.test → core
        for t in json_results[lib]['ALL_TESTS']: r[ key + '.' + t.replace(':', '.')] = _S_Failed_ if t in json_results[lib]['FAILED_TESTS'] else _S_Finished_

    results[_StateKey_]   = reduce(lambda a, b: _S_Finished_ if a==_S_Finished_ and b==_S_Finished_ else _S_Failed_, r.values())
    results[_LogKey_]     = output  # ommiting compilation log and only including integration.py output
    results[_ResultsKey_] = r
    return results



# do not change this wording, they have to stay in sync with upstream (up to benchmark-model).
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
