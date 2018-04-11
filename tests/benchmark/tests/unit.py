#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   unit.py
## @brief  Rosetta/PyRosetta unit tests
## @author Sergey Lyskov

import os, json, functools
import codecs

import imp
imp.load_source(__name__, '/'.join(__file__.split('/')[:-1]) +  '/__init__.py')  # A bit of Python magic here, what we trying to say is this: from __init__ import *, but init is calculated from file location

_api_version_ = '1.0'

_TestSuite_ = True  # Set to True for TestSuite-like tests (Unit, Integration, Sfxn_fingerprint) and False other wise


def set_up():  pass


def tear_down(): pass


def get_tests():
    TR = Tracer(verbose=True)
    TR('Unit Test script does not support get_tests! Use run_test_suite instead!')
    raise BenchmarkError()


def run_test(test, rosetta_dir, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    TR = Tracer(verbose)
    TR('Unit Test script does not support run_test! Use run_test_suite instead!')
    raise BenchmarkError()


def run_test_suite(rosetta_dir, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False, additional_flags="", mode="debug"):
    ''' Run TestSuite.
        Platform is a dict-like object, mandatory fields: {os='Mac', compiler='gcc'}
    '''

    jobs = config['cpu_count']

    TR = Tracer(verbose)
    full_log = ''

    TR('Running test_suite: "{}" at working_dir={working_dir!r} with rosetta_dir={rosetta_dir}, platform={platform}, mode={mode}, jobs={jobs}, hpc_driver={hpc_driver}...'.format(__name__, **vars() ) )

    results = {}

    TR('Compiling...')
    res, output, build_command_line = build_rosetta(rosetta_dir, platform, config, mode=mode, build_unit=True, verbose=verbose)

    full_log += 'Compiling: {}\n'.format(build_command_line) + output

    with codecs.open(working_dir+'/build-log.txt', 'w', encoding='utf-8', errors='replace') as f: f.write(full_log)

    if res:
        results[_StateKey_] = _S_build_failed_
        results[_LogKey_]   = full_log
        return results

    else:
        json_results_file = rosetta_dir+'/source/.unit_test_results.json'
        if (not debug) and  os.path.isfile(json_results_file): os.remove(json_results_file)

        compiler = platform['compiler']
        extras   = ','.join(platform['extras'])

        command_line = 'cd {}/source && test/run.py --compiler={compiler} --mode={mode} --extras={extras} -j{jobs} --mute all {additional_flags}'.format(rosetta_dir, jobs=jobs, compiler=compiler, extras=extras, mode=mode, additional_flags=additional_flags)
        TR( 'Running unit test script: {}'.format(command_line) )

        if debug: res, output = 0, 'unit.py: debug is enabled, skipping unit-tests script run...\n'
        else: res, output = execute('Running unit test script...', command_line, return_='tuple')
        full_log += output

        if res:
            results[_StateKey_] = _S_script_failed_
            results[_LogKey_]   = 'Compiling: {}\nRunning: {}\n'.format(build_command_line, command_line) + output  # ommiting compilation log and only including run.py output
            return results

    with open(json_results_file) as f: json_results = json.load(f)  # JSON handles unicode internally

    #r = {}
    # for lib in json_results:
    #     key = lib[:-5]  # core.test → core
    #     # u'∙'
    #     for t in json_results[lib]['ALL_TESTS']: r[ key.replace('.', '_') + '_' + t.replace(':', '_')] = _S_failed_ if t in json_results[lib]['FAILED_TESTS'] else _S_passed_

    results[_StateKey_]   = functools.reduce(lambda a, b: _S_passed_ if a==_S_passed_ and b==_S_passed_ else _S_failed_, [ json_results['tests'][t][_StateKey_] for t in json_results['tests'] ] )
    results[_LogKey_]     = output  # ommiting compilation log and only including unit tests output
    results[_ResultsKey_] = json_results

    with open(working_dir+'/unit.json', 'w') as f: json.dump(json_results, f, sort_keys=True, indent=2)

    return results


def run(test, rosetta_dir, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    if test == "valgrind": return run_test_suite(rosetta_dir, working_dir, platform, config, hpc_driver=hpc_driver, verbose=verbose, debug=debug, additional_flags="--valgrind --timeout 1440") # 24 hour time limit
    if test == "valgrind_detailed": return run_test_suite(rosetta_dir, working_dir, platform, config, hpc_driver=hpc_driver, verbose=verbose, debug=debug, additional_flags="--valgrind --trackorigins --timeout 1440") # 24 hour time limit
    elif test == "addsan":  return run_test_suite(rosetta_dir, working_dir, platform, config, hpc_driver=hpc_driver, verbose=verbose, debug=debug, mode="addsan")
    elif test == "ubsan":
        os.environ["UBSAN_OPTIONS"]="print_stacktrace=1" # Get the backtrace in the log when running ubsan
        return run_test_suite(rosetta_dir, working_dir, platform, config, hpc_driver=hpc_driver, verbose=verbose, debug=debug, mode="ubsan")
    elif test == 'release': return run_test_suite(rosetta_dir, working_dir, platform, config, hpc_driver=hpc_driver, verbose=verbose, debug=debug, mode='release')
    #elif test: return run_test(test, rosetta_dir, working_dir, platform, config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    elif test == '': return run_test_suite(rosetta_dir, working_dir, platform, config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    else: raise BenchmarkError('Unknow scripts test: {}!'.format(test))
