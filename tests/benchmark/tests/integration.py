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

import imp
imp.load_source(__name__, '/'.join(__file__.split('/')[:-1]) +  '/__init__.py')  # A bit of Python magic here, what we trying to say is this: from __init__ import *, but init is calculated from file location


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
        res, output = execute('Compiling...', 'cd {}/source && ./scons.py bin mode=release cxx={compiler} -j{jobs}'.format(rosetta_dir, jobs=jobs, compiler=platform['compiler']), return_='tuple')

    full_log += output  #file(working_dir+'/build-log.txt', 'w').write(output)

    if res:
        results[_StateKey_] = _S_build_failed_
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
        command_line = 'cd {}/tests/integration && ./integration.py --compiler={platform[compiler]} --timeout=480 -j{jobs}'.format(rosetta_dir, jobs=jobs, platform=platform)
        TR( 'Running integration script: {}'.format(command_line) )

        if debug: res, output = 0, 'integration.py: debug is enabled, skippig integration script run...\n'
        else:  res, output = execute('Running integration script...', command_line, return_='tuple')

        full_log += output

        if res:
            results[_StateKey_] = _S_script_failed_
            results[_LogKey_]   = output  # ommiting compilation log and only including integration.py output
            return results

    for d in os.listdir(files_location):
        if os.path.isdir(files_location + '/' + d):
            #print 'linking: %s <-- %s' % (root + d, working_dir + d)
            #os.symlink( os.path.abspath(files_location + '/' + d), working_dir + '/' + d)
            shutil.copytree(os.path.abspath(files_location + '/' + d), working_dir + '/' + d)

            command_sh = working_dir + '/' + d + '/command.sh '
            if os.path.isfile(command_sh): os.remove(command_sh)  # deleting non-tempalte command.sh files to avoid stroing absolute paths in database

    results[_StateKey_] = _S_queued_for_comparison_
    results[_LogKey_]   = output  # ommiting compilation log and only including integration.py output
    return results



# compare results of two tests run (new vs. previous)
# take two dict and two paths
# must return standard dict
def compare_tests(results, files_path, previous_results, previous_files_path):
    TR = Tracer(verbose=True)
    TR('Integration Test script does not support compare_tests! Use compare_test_suites instead!')
    raise BenchmarkError()


# compare results of two test suites run (new vs. previous)
# take two dict and two paths
# must return standard dict
def compare_test_suites(results, files_path, previous_results, previous_files_path):
    results = dict(tests={}, summary=dict(total=0, failed=0, failed_tests=[]), config={})

    for test in os.listdir(files_path):
        if os.path.isdir(files_path + '/' + test):
            res, brief_diff = execute('Comparing {}...'.format(test), 'diff -rq {0}/{test} {1}/{test}'.format(files_path, previous_files_path, test=test), return_='tuple')
            res, full_diff  = execute('Comparing {}...'.format(test), 'diff -r  {0}/{test} {1}/{test}'.format(files_path, previous_files_path, test=test), return_='tuple')
            results['tests'][test] = {_StateKey_: _S_failed_ if res else _S_finished_, _LogKey_: brief_diff + '\n\n' + full_diff[:16384] if res else ''}

            results['summary']['total'] += 1
            if res: results['summary']['failed'] += 1; results['summary']['failed_tests'].append(test)


    return {_StateKey_: _S_failed_ if results['summary']['failed'] else _S_finished_, _LogKey_: '', _ResultsKey_: results}
