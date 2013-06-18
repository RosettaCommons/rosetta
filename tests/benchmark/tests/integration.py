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

class BenchmarkIntegrationError(Exception): pass


tests = ['i']


def set_up():
    pass


def tear_down():
    pass


def rollover():
    a_commands.getoutput("cd main/tests/integration && ./accept-changes.sh" )


def get_tests():
    TR = Tracer(verbose=True)
    TR('Integration Test script does not support get_tests! Use run_test_suite instead!')
    raise BenchmarkIntegrationError()


def run_test(test, rosetta_dir, working_dir, jobs=1, hpc_driver=None, verbose=False):
    TR = Tracer(verbose)

    TR('Integration Test script does not support run_test! Use run_test_suite instead!')
    raise BenchmarkIntegrationError()


def run_test_suite(rosetta_dir, working_dir, jobs=1, hpc_driver=None, verbose=False):
    TR = Tracer(verbose)
    full_log = ''

    TR('Running test_suite: "{}" at working_dir={working_dir!r} with rosetta_dir={rosetta_dir} jobs={jobs}, hpc_driver={hpc_driver}...'.format(__name__, **vars() ) )

    results = {}

    TR('Compiling...')
    #res, output = execute('Compiling...', 'cd {}/source && ./scons.py bin mode=release -j{jobs}'.format(rosetta_dir, jobs=jobs), return_='tuple')
    res, output = 0, 'debug... compiling...\n'

    full_log += output  #file(working_dir+'/build-log.txt', 'w').write(output)

    if res:
        results['suite_result'] = _BuildFailed_
        results['suite_log']    = full_log
        return results

    else:
        TR('Running integration script...')
        #res, output = execute('Running integration script...', 'cd {}/tests/integration && ./integration.py -j{jobs}'.format(rosetta_dir, jobs=jobs), return_='tuple')
        res, output = 0, 'debug... integration.py...\n'
        full_log += output

        if res:
            results['suite_result'] = _ScriptFailed_
            results['suite_log']    = full_log
            return results

        files_location = rosetta_dir+'/tests/integration/new/'
        for root, dirs, _ in os.walk(files_location):
            if root == files_location:  # we only want to iterate over first layer
                for d in dirs:
                    #print 'linking: %s <-- %s' % (root + d, working_dir + d)
                    os.symlink(root + '/' + d, working_dir + '/' + d)

        results['suite_result'] = _Finished_
        results['suite_log']    = full_log
        return results


# do not change this wording, they have to stay in sync with upstream (up to benchmark-model).
_Finished_     = '_Finished_'
_Failed_       = '_Failed_'
_BuildFailed_  = '_BuildFailed_'
_ScriptFailed_ = '_ScriptFailed_'


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
