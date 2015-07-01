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

import os, os.path, json, commands, shutil, stat

import imp
imp.load_source(__name__, '/'.join(__file__.split('/')[:-1]) +  '/__init__.py')  # A bit of Python magic here, what we trying to say is this: from __init__ import *, but init is calculated from file location

_api_version_ = '1.0'  # api version


_failure_threshold_pct_ = 5  # specify how much execution time could deviate (percent) from previous value without raising the alarm


def run_performance_tests(rosetta_dir, working_dir, platform, config, hpc_driver, verbose=False, debug=False):
    memory = config['memory']
    jobs = config['cpu_count']
    mode = 'release'

    TR = Tracer(verbose)

    TR('Running Performance benchmark at working_dir={working_dir!r} with rosetta_dir={rosetta_dir}, platform={platform}, jobs={jobs}, memory={memory}GB, hpc_driver={hpc_driver}...'.format( **vars() ) )

    compiler = platform['compiler']
    extras   = ','.join(platform['extras'])

    build_command_line = './scons.py bin cxx={compiler} extras={extras} mode={mode} -j{jobs}'.format( **vars() )

    if debug: res, output = 0, 'run_performance_tests: debug is enabled, skippig build phase...\n'
    else: res, output = execute('Compiling...', 'cd {}/source && {}'.format(rosetta_dir, build_command_line), return_='tuple')

    file(working_dir+'/build-log.txt', 'w').write( 'Running: {}\n{}\n'.format(build_command_line, output) )

    results = {}

    if res:
        results[_StateKey_] = _S_build_failed_
        results[_LogKey_]   = 'Compiling: {}\n'.format(build_command_line) + output
        return results

    else:
        ext = calculate_extension(platform, mode)

        json_results_file = '{rosetta_dir}/source/_performance_'.format(**vars())
        output_log_file = '{working_dir}/_performance_.log'.format(**vars())
        command_line = '{rosetta_dir}/source/bin/performance_benchmark.{ext} -database {rosetta_dir}/database -benchmark_scale 120 -mute core protocols -in:file:extra_res_path extra_params 2>&1 >{output_log_file}'.format(**vars())

        #performance_benchmark_sh = os.path.abspath(working_dir + '/performance_benchmark.sh')
        #with file(performance_benchmark_sh, 'w') as f: f.write('#!/bin/bash\n{}\n'.format(command_line));  os.fchmod(f.fileno(), stat.S_IEXEC | stat.S_IREAD | stat.S_IWRITE)

        if debug and False: res, output = 0, 'run_performance_tests: debug is enabled, skipping actual run...\n'
        else:
            if os.path.isfile(json_results_file): os.remove(json_results_file)
            hpc_driver.execute(command_line, '{rosetta_dir}/source/src/apps/benchmark/performance'.format(**vars()), 'performance_benchmark')

        output = file(output_log_file).read()
        results[_LogKey_]   = 'Compiling: {}\nRunning: {}\n'.format(build_command_line, command_line) + output

        try:
            json_results = json.load( file(json_results_file) )
            results[_ResultsKey_] = { _TestsKey_:{} }
            for t in json_results: results[_ResultsKey_][_TestsKey_][t] = {_StateKey_: _S_queued_for_comparison_, _LogKey_: '', 'num_executions': json_results[t][0], 'run_time': json_results[t][1] }

            results[_ResultsKey_][_PlotsKey_] = [ dict(type='sub_test:revision_value', data=[ dict(y='num_executions', legend='num_executions, s', color='#66f') ] ) ]

            results[_StateKey_] = _S_queued_for_comparison_

        except IOError:
            results[_ResultsKey_] = { _TestsKey_:{} }
            results[_StateKey_] = _S_script_failed_

        return results



def run(test, rosetta_dir, working_dir, platform, config, hpc_driver, verbose=False, debug=False):
    ''' Run single test.
        Platform is a dict-like object, mandatory fields: {os='Mac', compiler='gcc'}
    '''
    if test =='': return run_performance_tests(rosetta_dir, working_dir, platform, config=config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    else: raise BenchmarkError('Unknow performance test: {}!'.format(test))



# compare results of two tests run (new vs. previous)
# take two dict and two paths
# must return standard dict with results
def compare(test, results, files_path, previous_results, previous_files_path):
    cr = {_TestsKey_:{}, _SummaryKey_:{_TotalKey_:0, _FailedKey_:0, _FailedTestsKey_:[]} }

    if previous_results[_ResultsKey_]  and  _TestsKey_ in previous_results[_ResultsKey_]:
        for test in results[_ResultsKey_][_TestsKey_]:
            num_executions = results[_ResultsKey_][_TestsKey_][test]['num_executions']

            if test in previous_results[_ResultsKey_][_TestsKey_]  and  'num_executions' in previous_results[_ResultsKey_][_TestsKey_][test]: previous_num_executions = previous_results[_ResultsKey_][_TestsKey_][test]['num_executions']
            else: previous_num_executions = None

            cr[_SummaryKey_][_TotalKey_] += 1
            cr[_TestsKey_][test] = {_StateKey_: _S_finished_, 'num_executions':num_executions, 'previous_num_executions':previous_num_executions,
                                    _LogKey_: 'previous_num_executions={}\nnum_executions={}\n'.format(previous_num_executions, num_executions) }

            if previous_num_executions  and  2.0 * abs(previous_num_executions - num_executions) / abs(previous_num_executions + num_executions + 1.0e-200) > _failure_threshold_pct_/100.0:  # mark test as failed if there is more then 5% difference in number of executions
                cr[_TestsKey_][test][_StateKey_] = _S_failed_
                cr[_SummaryKey_][_FailedKey_] += 1
                cr[_SummaryKey_][_FailedTestsKey_].append(test)

    else: # no previous tests case, returning 'finished' for all sub_tests
        for test in results[_ResultsKey_][_TestsKey_]:
            cr[_TestsKey_][test] = {_StateKey_: _S_finished_, 'run_time':results[_ResultsKey_][_TestsKey_][test]['num_executions'], 'previous_num_executions':None, _LogKey_: 'First run, no previous results available. Skipping comparison...\n'}
            cr[_SummaryKey_][_TotalKey_] += 1

    state = _S_failed_ if cr[_SummaryKey_][_FailedKey_] else _S_finished_

    cr[_PlotsKey_] = [ dict(type='sub_test:revision_value', data=[ dict(y='num_executions', legend='num_executions, s', color='#66f') ] ) ]

    return {_StateKey_: state, _LogKey_: '', _ResultsKey_: cr}
