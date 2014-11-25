#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   tests/performance.py
## @brief  Performace benchmark test
## @author Sergey Lyskov

import os, os.path, json, commands, shutil

import imp
imp.load_source(__name__, '/'.join(__file__.split('/')[:-1]) +  '/__init__.py')  # A bit of Python magic here, what we trying to say is this: from __init__ import *, but init is calculated from file location

_api_version_ = '1.0'  # api version


def run_performance_tests(rosetta_dir, working_dir, platform, config, hpc_driver, verbose=False, debug=False):
    memory = config['memory']
    jobs = config['cpu_count']
    mode = 'release'

    TR = Tracer(verbose)

    TR('Running Performance benchmark at working_dir={working_dir!r} with rosetta_dir={rosetta_dir}, platform={platform}, jobs={jobs}, memory={memory}GB, hpc_driver={hpc_driver}...'.format( **vars() ) )

    compiler = platform['compiler']
    extras   = ','.join(platform['extras'])

    build_command_line = './scons.py bin cxx={compiler} extras={extras} mode={mode} -j{jobs}'.format( **vars() )

    if debug: res, output = 0, 'build.py: debug is enabled, skippig build phase...\n'
    else: res, output = execute('Compiling...', 'cd {}/source && {}'.format(rosetta_dir, build_command_line), return_='tuple')

    if res:  res, output = execute('Compiling...', 'cd {}/source && {}'.format(rosetta_dir, build_command_line), return_='tuple')

    file(working_dir+'/build-log.txt', 'w').write( 'Running: {}\n{}\n'.format(build_command_line, output) )

    ext = calculate_extension(platform, mode)

    json_results_file = '{rosetta_dir}/source/_performance_'.format(**vars())
    command_line = '{rosetta_dir}/source/bin/performance_benchmark.{ext} -database {rosetta_dir}/database -mute core protocols -in:file:extra_res_path extra_params'.format(**vars())

    if not debug:
        if os.path.isfile(json_results_file): os.remove(json_results_file)
        hpc_driver.execute(command_line, '{rosetta_dir}/source/src/apps/benchmark/performance'.format(**vars()), 'performance_benchmark')

    json_results = json.load( file(json_results_file) )

    results = {}
    results[_StateKey_] = _S_queued_for_comparison_
    results[_LogKey_]   = 'Compiling: {}\nRunning: {}\n'.format(build_command_line, command_line) + test_output

    return results



def run(test, rosetta_dir, working_dir, platform, config, hpc_driver, verbose=False, debug=False):
    ''' Run single test.
        Platform is a dict-like object, mandatory fields: {os='Mac', compiler='gcc'}
    '''
    if   test =='': return run_performance_tests(rosetta_dir, working_dir, platform, config=config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    else: raise BenchmarkError('Unknow performance test: {}!'.format(test))
