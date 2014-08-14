#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   integration.py
## @brief  Rosetta integrtion tests
## @author Sergey Lyskov

import os, shutil, commands

import imp
imp.load_source(__name__, '/'.join(__file__.split('/')[:-1]) +  '/__init__.py')  # A bit of Python magic here, what we trying to say is this: from __init__ import *, but init is calculated from file location

ignore_files = 'command.sh observers'.split()

#tests = ['i']

_TestSuite_ = True  # Set to True for TestSuite-like tests (Unit, Integration, Sfxn_fingerprint) and False other wise


def set_up(): pass


def tear_down(): pass


#def rollover():
#    a_commands.getoutput("cd main/tests/integration && ./accept-changes.sh" )


def get_tests():
    raise BenchmarkError('Integration Test script does not support get_tests!')


def run_test(test, rosetta_dir, working_dir, platform, jobs=1, hpc_driver=None, verbose=False, debug=False):
    raise BenchmarkError('Integration Test script does not support run_test! Use run_test_suite instead!')


def run_itegration_tests(mode, rosetta_dir, working_dir, platform, jobs=1, hpc_driver=None, verbose=False, debug=False):
    ''' Run TestSuite.
        Platform is a dict-like object, mandatory fields: {os='Mac', compiler='gcc'}
    '''

    TR = Tracer(verbose)
    full_log = ''

    TR('Running test_suite: "{}" at working_dir={working_dir!r} with rosetta_dir={rosetta_dir}, platform={platform}, jobs={jobs}, hpc_driver={hpc_driver}...'.format(__name__, **vars() ) )

    results = {}

    TR('Compiling...')
    compiler = platform['compiler']
    extras   = ','.join(platform['extras'])

    # removing all symlinks from bin/ and then building binaries...
    build_command_line = 'find bin -type l ! -name ".*" -exec rm {{}} \\; && ./scons.py bin mode={mode} cxx={compiler} extras={extras} -j{jobs}'.format(jobs=jobs, mode=mode, compiler=compiler, extras=extras)

    if debug:
        res, output = 0, 'integration.py: debug is enabled, skippig build...\n'
    else:
        res, output = execute('Compiling...', 'cd {}/source && {}'.format(rosetta_dir, build_command_line), return_='tuple')

    full_log += output  #file(working_dir+'/build-log.txt', 'w').write(output)

    if res:
        results[_StateKey_] = _S_build_failed_
        results[_LogKey_]   = 'Compiling: {}\n'.format(build_command_line) + full_log
        return results

    else:
        ref_files_location = rosetta_dir+'/tests/integration/ref/'
        if os.path.isdir(ref_files_location): TR( 'Removing old ref dir {}...'.format(ref_files_location) );  shutil.rmtree(ref_files_location)
        TR('Creating a dummy ref dir {}...'.format(ref_files_location));
        os.mkdir(ref_files_location)

        files_location = rosetta_dir+'/tests/integration/new/'
        #if os.path.isdir(files_location): TR('Removing old ref dir %s...' % files_location);  shutil.rmtree(files_location)  # remove old dir if any

        #output_json = working_dir + '/output.json'  , output_json=output_json   --yaml={output_json}
        timeout = 60*8 if mode == 'release' else 60*60  # setting timeout to 8min on release and one hour on debug

        command_line = 'cd {}/tests/integration && ./integration.py --skip-comparison --mode={mode} --compiler={compiler} --extras={extras} --timeout={timeout} -j{jobs}'.format(rosetta_dir, jobs=jobs, mode=mode, compiler=compiler, extras=extras, timeout=timeout)
        TR( 'Running integration script: {}'.format(command_line) )

        if debug: res, output = 0, 'integration.py: debug is enabled, skippig integration script run...\n'
        else:  res, output = execute('Running integration script...', command_line, return_='tuple')

        full_log += output

        if res:
            results[_StateKey_] = _S_script_failed_
            results[_LogKey_]   = 'Compiling: {}\nRunning: {}\n'.format(build_command_line, command_line) + output  # ommiting compilation log and only including integration.py output
            return results

    ignore = []
    for d in os.listdir(files_location):
        if os.path.isdir(files_location + '/' + d):
            #print 'linking: %s <-- %s' % (root + d, working_dir + d)
            #os.symlink( os.path.abspath(files_location + '/' + d), working_dir + '/' + d)
            shutil.copytree(os.path.abspath(files_location + '/' + d), working_dir + '/' + d)

            #command_sh = working_dir + '/' + d + '/command.sh'
            #if os.path.isfile(command_sh): os.remove(command_sh)  # deleting non-tempalte command.sh files to avoid stroing absolute paths in database

            # for f in ignore_files:
            #     fl = working_dir + '/' + d + '/' + f
            #     if os.path.isfile(fl): ignore.append(d + '/' + f)  #os.remove(command_sh)  # deleting non-tempalte command.sh files to avoid stroing absolute paths in database


    results[_StateKey_]  = _S_queued_for_comparison_
    results[_LogKey_]    =  'Compiling: {}\nRunning: {}\n'.format(build_command_line, command_line) + output  # ommiting compilation log and only including integration.py output
    results[_IgnoreKey_] = ignore
    return results


def run(test, rosetta_dir, working_dir, platform, jobs=1, hpc_driver=None, verbose=False, debug=False):
    if not test:          return run_itegration_tests('release', rosetta_dir, working_dir, platform, jobs=jobs, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    elif test == 'debug': return run_itegration_tests('debug',   rosetta_dir, working_dir, platform, jobs=jobs, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    elif test == 'release_debug': return run_itegration_tests('release_debug',   rosetta_dir, working_dir, platform, jobs=jobs, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    else: raise BenchmarkError('Integration Test script does not support run with test="{}"!'.format(test))

    #if test: return run_test(test, rosetta_dir, working_dir, platform, jobs=jobs, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    #else: return run_test_suite(rosetta_dir, working_dir, platform, jobs=jobs, hpc_driver=hpc_driver, verbose=verbose, debug=debug)


# compare results of two tests run (new vs. previous)
# take two dict and two paths
# must return standard dict with results
def compare(test, results, files_path, previous_results, previous_files_path):
    #if test: raise BenchmarkError('Integration Test script does not support compare function for {} test!'.format(test))

    results = dict(tests={}, summary=dict(total=0, failed=0, failed_tests=[]), config={})

    if previous_files_path:
        for test in os.listdir(files_path):
            if os.path.isdir(files_path + '/' + test):
                exclude = ''.join([' --exclude="{}"'.format(f) for f in ignore_files] )
                res, brief_diff = execute('Comparing {}...'.format(test), 'diff -rq {exclude} {0}/{test} {1}/{test}'.format(previous_files_path, files_path, test=test, exclude=exclude), return_='tuple')
                res, full_diff  = execute('Comparing {}...'.format(test), 'diff -r  {exclude} {0}/{test} {1}/{test}'.format(previous_files_path, files_path, test=test, exclude=exclude), return_='tuple')
                diff = 'Brief Diff:\n' + brief_diff + ( ('\n\nFull Diff:\n' + full_diff[:1024*1024*16]) if full_diff != brief_diff else '' )
                results['tests'][test] = {_StateKey_: _S_failed_ if res else _S_finished_, _LogKey_: diff if res else ''}

                results['summary']['total'] += 1
                if res: results['summary']['failed'] += 1; results['summary']['failed_tests'].append(test)

    else: # no previous test failed, returning finished for all sub_tests
        for test in os.listdir(files_path):
            if os.path.isdir(files_path + '/' + test):
                results['tests'][test] = {_StateKey_: _S_finished_, _LogKey_: 'First run, no previous results available. Skipping comparison...\n'}
                results['summary']['total'] += 1

    return {_StateKey_: _S_failed_ if results['summary']['failed'] else _S_finished_, _LogKey_: '', _ResultsKey_: results}
