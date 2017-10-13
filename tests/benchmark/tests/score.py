#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   score.py
## @brief  Rosetta score function fingerprint tests
## @author Sergey Lyskov

import os, shutil, commands

import imp
imp.load_source(__name__, '/'.join(__file__.split('/')[:-1]) +  '/__init__.py')  # A bit of Python magic here, what we trying to say is this: from __init__ import *, but init is calculated from file location

_api_version_ = '1.0'

ignore_files = 'command.sh observers'.split()

def run_score_tests(mode, rosetta_dir, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    ''' Run TestSuite.
        Platform is a dict-like object, mandatory fields: {os='Mac', compiler='gcc'}
    '''

    jobs = config['cpu_count']

    TR = Tracer(verbose)
    full_log = ''

    TR('Running test_suite: "{}" at working_dir={working_dir!r} with rosetta_dir={rosetta_dir}, platform={platform}, jobs={jobs}, hpc_driver={hpc_driver}...'.format(__name__, **vars() ) )

    results = {}

    TR('Compiling...')
    res, output, build_command_line = build_rosetta(rosetta_dir, platform, config, mode=mode, verbose=verbose)
    full_log += output

    if res:
        results[_StateKey_] = _S_build_failed_
        results[_LogKey_]   = 'Compiling: {}\n'.format(build_command_line) + full_log
        return results

    else:
        ref_files_location = rosetta_dir+'/tests/sfxn_fingerprint/ref/'
        if os.path.isdir(ref_files_location): TR( 'Removing old ref dir {}...'.format(ref_files_location) );  shutil.rmtree(ref_files_location)
        TR('Creating a dummy ref dir {}...'.format(ref_files_location));
        os.mkdir(ref_files_location)

        files_location = rosetta_dir+'/tests/sfxn_fingerprint/new/'
        #if os.path.isdir(files_location): TR('Removing old ref dir %s...' % files_location);  shutil.rmtree(files_location)  # remove old dir if any

        #output_json = working_dir + '/output.json'  , output_json=output_json   --yaml={output_json}
        timeout = 60*8 if mode == 'release' else 60*60  # setting timeout to 8min on release and one hour on debug

        compiler = platform['compiler']
        extras   = ','.join(platform['extras'])

        command_line = 'cd {}/tests/sfxn_fingerprint && ./sfxn_fingerprint.py --skip-comparison --database=./../../database --mode={mode} --compiler={compiler} --extras={extras} --timeout={timeout} -j{jobs}'.format(rosetta_dir, jobs=jobs, mode=mode, compiler=compiler, extras=extras, timeout=timeout)
        TR( 'Running sfxn_fingerprint script: {}'.format(command_line) )

        if debug: res, output = 0, 'score.py: debug is enabled, skippig score script run...\n'
        else:  res, output = execute('Running sfxn_fingerprint script...', command_line, return_='tuple')

        full_log += output

        if res:
            results[_StateKey_] = _S_script_failed_
            results[_LogKey_]   = 'Compiling: {}\nRunning: {}\n'.format(build_command_line, command_line) + output  # ommiting compilation log and only including sfxn_fingerprint.py output
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
    results[_LogKey_]    =  'Compiling: {}\nRunning: {}\n'.format(build_command_line, command_line) + output  # ommiting compilation log and only including sfxn_fingerprint.py output
    results[_IgnoreKey_] = ignore
    return results


def run(test, rosetta_dir, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    if not test:          return run_score_tests('release', rosetta_dir, working_dir, platform, config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    elif test == 'debug': return run_score_tests('debug',   rosetta_dir, working_dir, platform, config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    else: raise BenchmarkError('Score Test script does not support run with test="{}"!'.format(test))

    #if test: return run_test(test, rosetta_dir, working_dir, platform, jobs=jobs, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    #else: return run_test_suite(rosetta_dir, working_dir, platform, jobs=jobs, hpc_driver=hpc_driver, verbose=verbose, debug=debug)


# compare results of two tests run (new vs. previous)
# take two dict and two paths
# must return standard dict with results
def compare(test, results, files_path, previous_results, previous_files_path):
    #if test: raise BenchmarkError('Score Test script does not support compare function for {} test!'.format(test))

    results = dict(tests={}, summary=dict(total=0, failed=0, failed_tests=[]), config={})
    has_failed_scripts = False

    if previous_files_path:
        for test in os.listdir(files_path):
            if os.path.isdir(files_path + '/' + test):
                exclude = ''.join([' --exclude="{}"'.format(f) for f in ignore_files] ) + ' --exclude="*.ignore"'
                res, brief_diff = execute('Comparing {}...'.format(test), 'diff -rq {exclude} {0}/{test} {1}/{test}'.format(previous_files_path, files_path, test=test, exclude=exclude), return_='tuple')
                res, full_diff  = execute('Comparing {}...'.format(test), 'diff -r  {exclude} {0}/{test} {1}/{test}'.format(previous_files_path, files_path, test=test, exclude=exclude), return_='tuple')
                diff = 'Brief Diff:\n' + brief_diff + ( ('\n\nFull Diff:\n' + full_diff[:1024*1024*1]) if full_diff != brief_diff else '' )

                if os.path.isfile(files_path+'/'+test+'/.test_did_not_run.log')  or  os.path.isfile(files_path+'/'+test+'/.test_got_timeout_kill.log'): state = _S_script_failed_;  has_failed_scripts=True
                else: state = _S_failed_ if res else _S_passed_
                results['tests'][test] = {_StateKey_: state, _LogKey_: diff if state != _S_passed_ else ''}

                results['summary']['total'] += 1
                if res: results['summary']['failed'] += 1; results['summary']['failed_tests'].append(test)

    else: # no previous test failed, returning 'passed' for all sub_tests
        for test in os.listdir(files_path):
            if os.path.isdir(files_path + '/' + test):
                results['tests'][test] = {_StateKey_: _S_passed_, _LogKey_: 'First run, no previous results available. Skipping comparison...\n'}
                results['summary']['total'] += 1

    #if has_failed_scripts: state = _S_script_failed_
    #else: state = _S_failed_ if results['summary']['failed'] else _S_passed_
    state = _S_failed_ if results['summary']['failed'] else _S_passed_

    return {_StateKey_: state, _LogKey_: '', _ResultsKey_: results}
