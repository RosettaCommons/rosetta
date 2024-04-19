#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   integration.py
## @brief  Rosetta integrtion tests
## @author Sergey Lyskov
## @author Rocco Moretti (Valgrind additions)
## @author Jared Adolf-Bryfogle( Demos/Tutorials additions)

import os, os.path, shutil, re, fnmatch
from collections import defaultdict
import json
import codecs

# A bit of Python magic here, what we trying to say is this: from __init__ import *, but init is calculated from file location
import importlib.util, sys
importlib.util.spec_from_file_location(__name__, '/'.join(__file__.split('/')[:-1]) +  '/__init__.py').loader.exec_module(sys.modules[__name__])

_api_version_ = '1.1'

ignore_files = 'command.sh command.mpi.sh observers'.split()


def run_test(test, rosetta_dir, working_dir, platform, jobs=1, hpc_driver=None, verbose=False, debug=False):
    raise BenchmarkError('Integration Test script does not support run_test! Use run_test_suite instead!')


def copy_files(files_location, working_dir):
    """
    Copy files to working directory for comparison or checking failures.
    :param files_location: str
    :param working_dir: str
    """
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

def run_general(mode, rosetta_dir, platform, jobs, TR, debug, full_log, build_command_line, additional_flags, results, config):
    """
    Run the typical integration tests for integration and demos

    Returns results and full_log so far.

    :param platform: dict-like object, mandatory fields: {os='Mac', compiler='gcc'}
    :param TR: Tracer
    :param debug: bool
    :param full_log: str
    :param build_command_line: str
    :param results: dict
    :rtype: str, str, dict, str
    """

    if re.search("--demos", additional_flags) or re.search("--tutorials", additional_flags):
        new_ref_files_location = os.path.realpath(rosetta_dir+'/demos')
    else:
        new_ref_files_location = rosetta_dir+'/tests/integration'

    for suffix in ['new', 'ref']:
        res_dir = new_ref_files_location + '/' + suffix

        if os.path.isdir(res_dir): TR(f'Removing old dir {res_dir}...');  shutil.rmtree(res_dir)
        TR(f'Creating an empty dir {res_dir}...')
        os.mkdir(res_dir)


    # Make a dummy ref directory and skip comparison if we're not in valgrind mode.
    if not re.search("--valgrind", additional_flags):
        additional_flags = "--skip-comparison " + additional_flags

    compiler = platform['compiler']
    extras   = ','.join(platform['extras'])
    #output_json = working_dir + '/output.json'  , output_json=output_json   --yaml={output_json}
    timeout = 60*8 if mode == 'release' else 60*60  # setting timeout to 8min on release and one hour on debug

    if re.search("--demos", additional_flags) or re.search("--tutorials", additional_flags):
        timeout = 60*64

    #Special case valgrind as it takes a *long* time - be very generous with the timeout.
    if re.search("--valgrind", additional_flags):
        timeout = 24*60*60  # If we've spent a full day on it, and it's still running, we're probably hosed.

    python_environment = local_python_install(platform, config)
    python_virtual_environment = setup_persistent_python_virtual_environment(python_environment, 'numpy scipy mmtf-python pubmed-lookup networkx mrcfile dask-jobqueue')

    command_line = f'{python_virtual_environment.activate} && cd {rosetta_dir}/tests/integration && python ./integration.py --mode={mode} --compiler={compiler} --extras={extras} --timeout={timeout} -j{jobs} {additional_flags}'

    TR( 'Running integration script: {}'.format(command_line) )

    #JAB - Why is this not run in debug mode, isn't this half the point?
    # RM - the "debug" variable is debug of the benchmark machinery. Debug mode is mode=="debug"
    if debug: res, output = 0, 'integration.py: debug is enabled, skipping integration script run...\n'
    else: res, output = execute('Running integration script...', command_line, return_='tuple', add_message_and_command_line_to_output=True)
    full_log += output

    if res:
        results[_StateKey_] = _S_script_failed_
        results[_LogKey_]   = 'Compiling: {}\nRunning: {}\n'.format(build_command_line, command_line) + output  # omitting compilation log and only including integration.py output

    return command_line, output, results, full_log


######################################################
def run_integration_tests(mode, rosetta_dir, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False, additional_flags=''):
    ''' Run TestSuite.
        Platform is a dict-like object, mandatory fields: {os='Mac', compiler='gcc'}
    '''

    jobs = config['cpu_count']

    TR = Tracer(verbose)
    full_log = ''

    TR('Running test_suite: "{}" at working_dir={working_dir!r} with rosetta_dir={rosetta_dir}, platform={platform}, jobs={jobs}, hpc_driver={hpc_driver}...'.format(__name__, **vars() ) )

    results = defaultdict()

    TR('Compiling...')
    res, output, build_command_line = build_rosetta(rosetta_dir, platform, config, mode=mode, verbose=verbose)

    full_log += output  #file(working_dir+'/build-log.txt', 'w').write(output)

    if res:
        results[_StateKey_] = _S_build_failed_
        results[_LogKey_]   = 'Compiling: {}\n'.format(build_command_line) + full_log
        return results

    # For trajectory consistency from run-to-run we want to have the Dunbrack libraries generated with the particular mode (release/debug) that we're running with.
    os_spec = platform['os']
    compiler = platform['compiler']
    extras = ''.join( sorted(platform['extras']) )
    if extras:
        extras += '.'
    database_cache_path = f"{rosetta_dir}/.database-binaries/{extras}{os_spec}{compiler}{mode}"
    if not os.path.exists(database_cache_path):
        os.makedirs(database_cache_path)
    additional_flags += f' --additional_flags "-in:path:database_cache_dir {database_cache_path}" ' # Quoted to keep it together

    #if os.path.isdir(files_location): TR('Removing old ref dir %s...' % files_location);  shutil.rmtree(files_location)  # remove old dir if any
    command_line, output, results, full_log = run_general(mode, rosetta_dir, platform, jobs, TR, debug, full_log, build_command_line, additional_flags, results, config=config)

    ignore = []
    files_location = rosetta_dir+'/tests/integration/new/'
    copy_files(files_location, working_dir)

    with open(working_dir + '/full-log.txt', 'wb') as f: f.write( to_bytes(full_log) )

    if _StateKey_ not in results: # it possible that failure was already reported by setting state to _S_failed_ or even _S_script_failed_
        if mode in ['ubsan','addsan']:
            # Don't queue for comparison - sanitizers are hard pass/fail.
            results =  gather_pass_fail_results(files_location)
        else:
            results[_StateKey_]  = _S_queued_for_comparison_ #This queues the comparison for the test server, which will run compare function in this file.

        results[_LogKey_]    =  'Compiling: {}\nRunning: {}\n'.format(build_command_line, command_line) + output  # omitting compilation log and only including integration.py output
        results[_IgnoreKey_] = ignore
    return results


def run_demo_tests(mode, rosetta_dir, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False, additional_flags='', diff_causes_failure=False):
    ''' Run TestSuite for the Demos and Tutorials
        Platform is a dict-like object, mandatory fields: {os='Mac', compiler='gcc'}

        IF  diff_causes_failure=True, this indicates to use the DIFF for failure.  Here, we check to make sure it all runs without error.
        Integration tests proper should check the diffs.  If a demo or tutorial is not covered sufficiently in integration and unit tests, so be it.
    '''

    jobs = config['cpu_count']

    TR = Tracer(verbose)
    full_log = ''

    TR('Running test_suite: "{}" at working_dir={working_dir!r} with rosetta_dir={rosetta_dir}, platform={platform}, jobs={jobs}, hpc_driver={hpc_driver}...'.format(__name__, **vars() ) )

    results = defaultdict()

    TR('Compiling...')
    res, output, build_command_line = build_rosetta(rosetta_dir, platform, config, mode=mode, verbose=verbose)
    full_log += output

    if res:
        results[_StateKey_] = _S_build_failed_
        results[_LogKey_]   = 'Compiling: {}\n'.format(build_command_line) + full_log
        return results

    command_line, output, results, full_log = run_general(mode, rosetta_dir, platform, jobs, TR, debug, full_log, build_command_line, additional_flags, results, config=config)

    ignore = []

    files_location = os.path.realpath(rosetta_dir+'/demos/new')
    copy_files(files_location, working_dir)

    with open(working_dir + '/full-log.txt', 'wb') as f: f.write( to_bytes(full_log) )

    #Do a DIFF to indicate pass or f
    if diff_causes_failure:
        results[_StateKey_]  = _S_queued_for_comparison_
        results[_LogKey_]    =  'Compiling: {}\nRunning: {}\n'.format(build_command_line, command_line) + output  # omitting compilation log and only including integration.py output
        results[_IgnoreKey_] = ignore
        return results

    else:
        results =  gather_pass_fail_results(files_location)
        results[_LogKey_]    =  'Compiling: {}\nRunning: {}\n'.format(build_command_line, command_line) + output
        return results


def run_valgrind_tests(mode, rosetta_dir, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False, additional_flags=''):
    ''' Run TestSuite under valgrind
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

    files_location = rosetta_dir+'/tests/integration/valgrind/'
    json_results_file = os.path.abspath( files_location +"/valgrind_results.yaml" )
    command_line, output, results, full_log = run_general(mode, rosetta_dir, platform, jobs, TR, debug, full_log, build_command_line, additional_flags, results, config=config)

    ignore = []
    # Copy just the valgrind output to the archive directory - we don't need the other output
    for d, dirnames, filenames in os.walk(files_location):
        for filename in filenames:
            keep_file = filename.endswith("valgrind.out")
            keep_file = keep_file or filename.startswith(".test") # Also keep .test_did_not_run.log or .test_got_timeout_kill.log

            # For debugging purposes, keep the expanded command shell and the log files
            keep_file = keep_file or filename == 'command.sh' or 'log' in filename

            if keep_file:
                relpath = os.path.relpath(d, files_location)
                try:
                    os.makedirs(working_dir + '/' + relpath)
                except OSError:
                    if not os.path.isdir(working_dir + '/' + relpath):
                        raise
                shutil.copy(os.path.abspath(files_location + '/' + relpath + '/' + filename), working_dir + '/' + relpath + '/' + filename)

    json_file_results = json.load( open( json_results_file ) ) #JSON handles unicode internally
    if json_file_results[ "failed" ] > 0:
        results[_StateKey_] = _S_failed_
    else:
        results[_StateKey_] = _S_passed_

    json_results = dict(tests={}, summary=dict(total=json_file_results[ "total" ], failed=json_file_results[ "failed" ], failed_tests=[]))
    for test, nfailures in json_file_results[ "details" ].items():
        log = ''
        if nfailures > 0 or os.path.isfile(files_location+'/'+test+'/.test_did_not_run.log')  or  os.path.isfile(files_location+'/'+test+'/.test_got_timeout_kill.log'):
            if nfailures > 0:
                state = _S_failed_
                log = "Found {} Valgrind error(s).\n\n".format(nfailures)
            else:
                state = _S_script_failed_
                log = "Test script did not run correctly.\n\n"
            json_results['summary']['failed_tests'].append(test)
            if os.path.isfile(files_location+'/'+test+'/valgrind.out'):
                log += codecs.open(files_location+'/'+test+'/valgrind.out', 'r', encoding='utf-8', errors='backslashreplace').read()
        else:
            state = _S_passed_
        json_results['tests'][test] = {_StateKey_: state, _LogKey_: log }

    with open(working_dir + '/full-log.txt', 'wb') as f: f.write( to_bytes(full_log) )

    results[_LogKey_]    =  'Compiling: {}\nRunning: {}\n'.format(build_command_line, command_line) + output  # omitting compilation log and only including integration.py output
    results[_IgnoreKey_] = ignore
    results[_ResultsKey_] = json_results
    return results


def compare(test, results, files_path, previous_results, previous_files_path):
    """
    Compare the results of two tests run (new vs. previous) for integration test
    Take two dict and two paths
    Must return standard dict with results

    :param test: str
    :param results: dict
    :param files_path: str
    :param previous_results: dict
    :param previous_files_path: str
    :rtype: dict
    """
    #if test: raise BenchmarkError('Integration Test script does not support compare function for {} test!'.format(test))

    results = dict(tests={}, summary=dict(total=0, failed=0, failed_tests=[]))  # , config={}

    if previous_files_path:
        for test in os.listdir(files_path):
            if os.path.isdir(files_path + '/' + test):
                exclude = ''.join([' --exclude="{}"'.format(f) for f in ignore_files] ) + ' --exclude="*.ignore"'
                res, brief_diff = execute('Comparing {}...'.format(test), 'diff -rq {exclude} {0}/{test} {1}/{test}'.format(previous_files_path, files_path, test=test, exclude=exclude), return_='tuple')
                res, full_diff  = execute('Comparing {}...'.format(test), 'diff -r  {exclude} {0}/{test} {1}/{test}'.format(previous_files_path, files_path, test=test, exclude=exclude), return_='tuple')
                diff = 'Brief Diff:\n' + brief_diff + ( ('\n\nFull Diff:\n' + full_diff[:1024*1024*1]) if full_diff != brief_diff else '' )

                state = _S_failed_ if res else _S_passed_
                results['tests'][test] = {_StateKey_: state, _LogKey_: diff if state != _S_passed_ else ''}

                results['summary']['total'] += 1
                if res: results['summary']['failed'] += 1; results['summary']['failed_tests'].append(test)

    else: # no previous tests case, returning 'passed' for all sub_tests
        for test in os.listdir(files_path):
            if os.path.isdir(files_path + '/' + test):
                results['tests'][test] = {_StateKey_: _S_passed_, _LogKey_: 'First run, no previous results available. Skipping comparison...\n'}
                results['summary']['total'] += 1

    for test in os.listdir(files_path):
        if os.path.isdir(files_path + '/' + test):
            if os.path.isfile(files_path+'/'+test+'/.test_did_not_run.log')  or  os.path.isfile(files_path+'/'+test+'/.test_got_timeout_kill.log'):
                results['tests'][test][_StateKey_] = _S_script_failed_
                results['tests'][test][_LogKey_] += '\nCompare(...): Marking as "Script failed" due to presense of .test_did_not_run.log or .test_got_timeout_kill.log file!\n'
                if test not in results['summary']['failed_tests']:
                    results['summary']['failed'] += 1
                    results['summary']['failed_tests'].append(test)

    state = _S_failed_ if results['summary']['failed'] else _S_passed_

    return {_StateKey_: state, _LogKey_: '', _ResultsKey_: results}

def gather_pass_fail_results(files_path):

    results = dict(tests={}, summary=dict(total=0, failed=0, failed_tests=[]))  # , config={}

    for test in os.listdir(files_path):
        if os.path.isdir(files_path + '/' + test):

            if os.path.isfile(files_path+'/'+test+'/.test_did_not_run.log')  or  \
               os.path.isfile(files_path+'/'+test+'/.test_got_timeout_kill.log'): state = _S_script_failed_;  has_failed_scripts=True
            else: state = _S_passed_

            # Get all log files recursively.
            log_files = []
            for root, dirnames, filenames in os.walk(os.path.join(files_path, test)):
                for filename in fnmatch.filter(filenames, 'log*'):
                    log_files.append(os.path.join(root, filename))

            log_files.sort(key=sort_natural)
            out_log = ""
            for log in log_files:
                FILE = codecs.open(log, 'r', encoding='utf-8', errors='backslashreplace') # Be robust to unicode in the log file
                lines = FILE.readlines()
                if lines:
                    out_log += "\n\nLogFile "+os.path.basename(log)+"\n"
                    out_log += "".join(lines)
                FILE.close()



            results['tests'][test] = {_StateKey_: state, _LogKey_: out_log if state != _S_passed_ else ''}

            results['summary']['total'] += 1
            if state == _S_script_failed_:
                results['summary']['failed'] += 1
                results['summary']['failed_tests'].append(test)




    #if has_failed_scripts: state = _S_script_failed_
    #else: state = _S_failed_ if results['summary']['failed'] else _S_passed_
    state = _S_failed_ if results['summary']['failed'] else _S_passed_

    return {_StateKey_: state, _LogKey_: '', _ResultsKey_: results}




def atoi(text):
    return int(text) if text.isdigit() else text

def sort_natural(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]



def run(test, repository_root, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    if not test:                             return run_integration_tests('release',                  repository_root, working_dir, platform, config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    elif test == 'debug':                    return run_integration_tests('debug',                    repository_root, working_dir, platform, config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    elif test == 'release_debug':            return run_integration_tests('release_debug',            repository_root, working_dir, platform, config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    elif test == 'release_debug_no_symbols': return run_integration_tests('release_debug_no_symbols', repository_root, working_dir, platform, config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    elif test == 'mpi':                      return run_integration_tests('release_debug',            repository_root, working_dir, platform, config, hpc_driver=hpc_driver, verbose=verbose, debug=debug, additional_flags='--suffix mpi')
    elif test == 'torch':                    return run_integration_tests('debug',                    repository_root, working_dir, platform, config, hpc_driver=hpc_driver, verbose=verbose, debug=debug, additional_flags='--suffix torch')
    elif test == 'tensorflow':               return run_integration_tests('debug',                    repository_root, working_dir, platform, config, hpc_driver=hpc_driver, verbose=verbose, debug=debug, additional_flags='--suffix tensorflow')
    elif test == 'thread':                   return run_integration_tests('debug',                    repository_root, working_dir, platform, config, hpc_driver=hpc_driver, verbose=verbose, debug=debug, additional_flags='--suffix thread')
    elif test == 'bcl':                      return run_integration_tests('debug',                    repository_root, working_dir, platform, config, hpc_driver=hpc_driver, verbose=verbose, debug=debug, additional_flags='--suffix bcl')
    elif test == 'addsan':                   return run_integration_tests('addsan',                   repository_root, working_dir, platform, config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    elif test == 'ubsan':
        os.environ["UBSAN_OPTIONS"]="print_stacktrace=1" # Get the backtrace in the log when running ubsan
        return run_integration_tests('ubsan',   repository_root, working_dir, platform, config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    elif test == 'valgrind':          return run_valgrind_tests('release_symbols', repository_root, working_dir, platform, config, hpc_driver=hpc_driver, verbose=verbose, debug=debug, additional_flags='--valgrind --yaml valgrind/valgrind_results.yaml') # 'release_symbols' for line # information
    elif test == 'valgrind_detailed': return run_valgrind_tests('release_symbols', repository_root, working_dir, platform, config, hpc_driver=hpc_driver, verbose=verbose, debug=debug, additional_flags='--valgrind --yaml valgrind/valgrind_results.yaml --trackorigins') # 'release_symbols' for line # information
    elif test == 'demos':                    return run_demo_tests('release', repository_root, working_dir, platform, config, hpc_driver=hpc_driver, verbose=verbose, debug=debug, additional_flags='--demos')
    elif test == 'tutorials':                return run_demo_tests('release', repository_root, working_dir, platform, config, hpc_driver=hpc_driver, verbose=verbose, debug=debug, additional_flags='--tutorials')
    else: raise BenchmarkError('Integration Test script does not support run with test="{}"!'.format(test))

    #if test: return run_test(test, repository_root, working_dir, platform, jobs=jobs, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    #else: return run_test_suite(repository_root, working_dir, platform, jobs=jobs, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
