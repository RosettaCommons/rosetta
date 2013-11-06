#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   integration.py
## @brief  Rosetta/PyRosetta integrtion tests
## @author Sergey Lyskov

import os, shutil, commands

class BenchmarkIntegrationError(Exception): pass


#tests = ['i']

_TestSuite_ = True  # Set to True for TestSuite-like tests (Unit, Integration, Sfxn_fingerprint) and False other wise


def set_up(): pass


def tear_down(): pass


#def rollover():
#    a_commands.getoutput("cd main/tests/integration && ./accept-changes.sh" )


def get_tests():
    TR = Tracer(verbose=True)
    TR('Integration Test script does not support get_tests! Use run_test_suite instead!')
    raise BenchmarkIntegrationError()


def run_test(test, rosetta_dir, working_dir, platform, jobs=1, hpc_driver=None, verbose=False, debug=False):
    TR = Tracer(verbose)
    TR('Integration Test script does not support run_test! Use run_test_suite instead!')
    raise BenchmarkIntegrationError()


def run_test_suite(rosetta_dir, working_dir, platform, jobs=1, hpc_driver=None, verbose=False, debug=False):
    ''' Run TestSuite.
        Platform is a dict-like object, mandatory fields: {os='Mac', compiler='gcc'}
    '''

    TR = Tracer(verbose)
    full_log = ''

    TR('Running test_suite: "{}" at working_dir={working_dir!r} with rosetta_dir={rosetta_dir}, platform={platform}, jobs={jobs}, hpc_driver={hpc_driver}...'.format(__name__, **vars() ) )

    results = {}

    TR('Compiling...')

    if debug:
        res, output = 0, 'integration.py: debug is enabled, skippig build...\n'
    else:
        res, output = execute('Compiling...', 'cd {}/source && ./scons.py bin mode=release -j{jobs}'.format(rosetta_dir, jobs=jobs), return_='tuple')

    full_log += output  #file(working_dir+'/build-log.txt', 'w').write(output)

    if res:
        results[_StateKey_] = _S_BuildFailed_
        results[_LogKey_]   = full_log
        return results

    else:
        ref_files_location = rosetta_dir+'/tests/integration/ref/'
        if os.path.isdir(ref_files_location): TR( 'Removing old ref dir {}...'.format(ref_files_location) );  shutil.rmtree(ref_files_location)
        TR('Creating a dummy ref dir {}...'.format(ref_files_location));
        os.mkdir(ref_files_location)

        files_location = rosetta_dir+'/tests/integration/new/'
        #if os.path.isdir(files_location): TR('Removing old ref dir %s...' % files_location);  shutil.rmtree(files_location)  # remove old dir if any

        #output_json = working_dir + '/output.json'  , output_json=output_json   --yaml={output_json}
        command_line = 'cd {}/tests/integration && ./integration.py --timeout=480 -j{jobs}'.format(rosetta_dir, jobs=jobs)
        TR( 'Running integration script: {}'.format(command_line) )

        if debug: res, output = 0, 'integration.py: debug is enabled, skippig integration script run...\n'
        else:  res, output = execute('Running integration script...', command_line, return_='tuple')

        full_log += output

        if res:
            results[_StateKey_] = _S_ScriptFailed_
            results[_LogKey_]   = output  # ommiting compilation log and only including integration.py output
            return results

    for d in os.listdir(files_location):
        if os.path.isdir(files_location + '/' + d):
            #print 'linking: %s <-- %s' % (root + d, working_dir + d)
            #os.symlink( os.path.abspath(files_location + '/' + d), working_dir + '/' + d)
            shutil.copytree(os.path.abspath(files_location + '/' + d), working_dir + '/' + d)

            command_sh = working_dir + '/' + d + '/command.sh '
            if os.path.isfile(command_sh): os.remove(command_sh)  # deleting non-tempalte command.sh files to avoid stroing absolute paths in database

    results[_StateKey_] = _S_QueuedForComparison_
    results[_LogKey_]   = output  # ommiting compilation log and only including integration.py output
    return results



# ⚔ do not change this wording, they have to stay in sync with upstream (up to benchmark-model).
# Copied from benchmark-model, standard state code's for tests results.
_S_Draft_               = 'draft'
_S_Queued_              = 'queued'
_S_Running_             = 'running'
_S_Finished_            = 'finished'
_S_Failed_              = 'failed'
_S_BuildFailed_         = 'build failed'
_S_ScriptFailed_        = 'script failed'
_S_QueuedForComparison_ = 'queued for comparison'

_S_Values_ = [_S_Draft_, _S_Queued_, _S_Running_, _S_Finished_, _S_Failed_, _S_BuildFailed_, _S_ScriptFailed_, _S_QueuedForComparison_]

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
