#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   profile.py
## @brief  profile.py
## Benchmark script for running Rosetta profile tests
## @author Sergey Lyskov

import json, os, os.path, shutil
import codecs

import imp
imp.load_source(__name__, '/'.join(__file__.split('/')[:-1]) +  '/__init__.py')  # A bit of Python magic here, what we trying to say is this: from ../__init__ import *, but init path is calculated relatively to this location


_api_version_ = '1.0'  # api version


_failure_threshold_min_execution_time_   = 128
_failure_threshold_min_alloacted_memory_ = 256

_failure_threshold_execution_time_pct_       = 25
_failure_threshold_max_memory_allocated_pct_ = 4


def run(test, rosetta_dir, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):

    jobs = config['cpu_count']

    # Building Rosetta binaries
    res, output, build_command_line = build_rosetta(rosetta_dir, platform, config)
    if res: return { _StateKey_ : _S_build_failed_,  _ResultsKey_ : {},
                     _LogKey_ : u'Building rosetta failed!\n{}\n{}\n'.format(build_command_line, output) }
    else:
        extension = calculate_extension(platform)
        compiler  = platform['compiler']

        json_results_file = '{rosetta_dir}/tests/profile/.profile_test_results.json'.format( **vars() )
        if os.path.isfile(json_results_file): os.remove(json_results_file)  # removing old results file if it is present

        output_log_file = '{working_dir}/profile_py.log'.format( **vars() )

        # running profile script on HPC cluster
        # hpc_driver.execute(executable='cd'.format(**vars()), arguments='{rosetta_dir}/tests/profile && ./profile.py --compiler={compiler} --daemon 2>&1 >{output_log_file}'.format(output_log_file=output_log_file, rosetta_dir=rosetta_dir, compiler=compiler),  # relax_native ligand_dock_script fixbb jd2test
        #                    working_dir=working_dir, name='profile', shell_wrapper=True)

        execute('Running Profile tests...', 'cd {rosetta_dir}/tests/profile && ./profile.py --compiler={compiler} --daemon 2>&1 >{output_log_file}'.format(output_log_file=output_log_file, rosetta_dir=rosetta_dir, compiler=compiler) )

        files_location = '{rosetta_dir}//tests/profile/tests'.format( **vars() )

        for d in os.listdir(files_location):
            if os.path.isdir(files_location + '/' + d)  and  os.path.isdir(files_location + '/' + d + '/output'):
                shutil.copytree(os.path.abspath(files_location + '/' + d + '/output'), working_dir + '/' + d)
            # if os.path.isdir(files_location + '/' + d):
            #     shutil.copytree(os.path.abspath(files_location + '/' + d), working_dir + '/' + d)

        with open(json_results_file) as f: tests = json.load(f) #JSON handles unicode internally

        for t in tests: tests[t][_LogKey_]   = ''  #tests[t][_StateKey_] = _S_queued_for_comparison_


        results = { _TestsKey_ : tests }

        return {_StateKey_ : _S_queued_for_comparison_,  _ResultsKey_ : results,  _LogKey_ : codecs.open(output_log_file, encoding='utf-8', errors='replace').read() }



# compare results of two tests run (new vs. previous)
# take two dict and two paths
# must return standard dict with results
def compare(test, results, files_path, previous_results, previous_files_path):
    cr = {_TestsKey_:{}, _SummaryKey_:{_TotalKey_:0, _FailedKey_:0, _FailedTestsKey_:[]} }

    if previous_results  and  _TestsKey_ in previous_results:
        for test in results[_TestsKey_]:
            execution_time       = results[_TestsKey_][test]['execution_time']
            max_memory_allocated = results[_TestsKey_][test]['max_memory_allocated']

            if test in previous_results[_TestsKey_]: previous_values = { k:previous_results[_TestsKey_][test][k] for k in previous_results[_TestsKey_][test] if k not in [_StateKey_, _LogKey_, 'previous_values'] }
            else: previous_values = None

            if previous_values  and 'execution_time' in previous_values  and 'max_memory_allocated' in previous_values:
                previous_execution_time = previous_values['execution_time']
                previous_max_memory_allocated = previous_values['max_memory_allocated']
            else:
                previous_execution_time = None
                previous_max_memory_allocated = None


            log = 'execution_time={}       previous_execution_time={}\nmax_memory_allocated={} previous_max_memory_allocated={}\n'.format(execution_time, previous_execution_time, max_memory_allocated, previous_max_memory_allocated)


            if _StateKey_ in results[_TestsKey_][test]:
                state = results[_TestsKey_][test][_StateKey_]

                if state == _S_failed_:
                    state = _S_script_failed_
                    log = 'Sub test script terminated with non-zero exit status, marking test as script-failed!\n\n' + log

                elif execution_time < _failure_threshold_min_execution_time_:
                    state = _S_failed_
                    log = 'Execution time is below {}s, marking test as failed!\n\n'.format(_failure_threshold_min_execution_time_) + log

                elif max_memory_allocated < _failure_threshold_min_alloacted_memory_:
                    state = _S_failed_
                    log = 'Max memory allocated is below {}Mb, marking test as failed!\n\n'.format(_failure_threshold_min_alloacted_memory_) + log

                elif previous_execution_time  and 2.0 * abs(previous_execution_time - execution_time) / abs(previous_execution_time + execution_time + 1.0e-200) > _failure_threshold_execution_time_pct_/100.0:
                    state = _S_failed_
                    log = 'Execution time is deviate more then {}% from previous resuls, marking test as failed!\n\n'.format(_failure_threshold_execution_time_pct_) + log

                elif previous_max_memory_allocated  and 2.0 * abs(previous_max_memory_allocated - max_memory_allocated) / abs(previous_max_memory_allocated + max_memory_allocated + 1.0e-200) > _failure_threshold_max_memory_allocated_pct_/100.0:
                    state = _S_failed_
                    log = 'Max memory allocated is deviate more then {}% from previous resuls, marking test as failed!\n\n'.format(_failure_threshold_max_memory_allocated_pct_) + log

            else:
                state = _S_script_failed_
                log = 'Sub test result does not have "state" key, marking test as script-failed!\n\n' + log


            cr[_SummaryKey_][_TotalKey_] += 1

            cr[_TestsKey_][test] = dict(results[_TestsKey_][test])
            cr[_TestsKey_][test].update( {_StateKey_: state, 'previous_values': previous_values, _LogKey_: log } )

            if state != _S_passed_:
                cr[_SummaryKey_][_FailedKey_] += 1
                cr[_SummaryKey_][_FailedTestsKey_].append(test)

    else: # no previous tests case, returning 'passed' for all sub_tests
        for test in results[_TestsKey_]:
            cr[_TestsKey_][test] = dict(results[_TestsKey_][test])

            if results[_TestsKey_][test][_StateKey_] != _S_passed_:
                cr[_TestsKey_][test].update( {_StateKey_: _S_script_failed_, 'previous_values':None, _LogKey_: 'Test bash script termiated with error! Skipping comparison...\n'} )
                cr[_SummaryKey_][_FailedKey_] += 1
                cr[_SummaryKey_][_FailedTestsKey_].append(test)

            else:
                cr[_TestsKey_][test].update( {_StateKey_: _S_passed_, 'previous_values':None, _LogKey_: 'First run, no previous results available. Skipping comparison...\n'} )

            cr[_SummaryKey_][_TotalKey_] += 1

    state = _S_failed_ if cr[_SummaryKey_][_FailedKey_] else _S_passed_

    # Two values on same plot:
    #cr[_PlotsKey_] = [ dict(type='sub_test:revision_value', data=[ dict(y='execution_time',       legend='execution_time',       color='#66f'),
    #                                                               dict(y='max_memory_allocated', legend='max_memory_allocated', color='#6f6') ] ) ]

    cr[_PlotsKey_] = [ dict(type='sub_test:revision_value', data=[ dict(y='execution_time',       legend='execution_time(s)/revision',       color='#66f') ] ),
                       dict(type='sub_test:revision_value', data=[ dict(y='max_memory_allocated', legend='max_memory_allocated(MiB)/revision', color='#6f6') ] ),
                       dict(type='sub_test:value_value', data=[ dict(array='memory_usage', x='time', y='memory', legend='memory(MiB)/time(s)', color='#c33') ] ) ]


    return {_StateKey_: state, _LogKey_: '', _ResultsKey_: cr}


    #failed_tests = [ t for t in tests if tests[t][_StateKey_] == _S_failed_ ]
