#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   serialization.py
## @brief  scientific/serialization.py
## Benchmark script for running serialization test
## @author Sergey Lyskov

import os, json

import imp
imp.load_source(__name__, '/'.join(__file__.split('/')[:-1]) +  '/__init__.py')  # A bit of Python magic here, what we trying to say is this: from ../__init__ import *, but init path is calculated relatively to this location


_api_version_ = '1.0'  # api version


def run_serialization_test(rosetta_dir, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    jobs = config['cpu_count']

    execute('Updating options, ResidueTypes and version info...', 'cd {}/source && ./update_options.sh && ./update_ResidueType_enum_files.sh && python version.py'.format(rosetta_dir) )

    executable = install_llvm_tool('serialization_validator', source_location='{}/../tools/clang_ast_transform/rosetta-refactor-tool'.format(rosetta_dir), config=config)
    #executable_path = executable.rpartition('/')[0]

    tools_dir = os.path.abspath(rosetta_dir+'/../tools')

    json_output_dir = working_dir + '/serialization-validator-output'
    os.mkdir( json_output_dir )

    command_line = 'export PYTHONPATH={tools_dir}/python_cc_reader${{PYTHONPATH+:$PYTHONPATH}} && cd {rosetta_dir}/source && python ../../tools/clang_ast_transform/run_on_all_ccfiles_w_fork.py -e '\
    '"python  {rosetta_dir}/../tools/clang_ast_transform/run_serialization_validator_on_file.py --executable_path {executable} --json_output_path {json_output_dir} --filename" -n {jobs}'.format(**vars())

    res, output = execute('Running...', command_line, return_='tuple')
    output = 'Running: ' + command_line + '\n' + output + '\n\n'

    json_files = [f for f in os.listdir(json_output_dir) if f.endswith('.json')]
    if json_files:
        state = _S_failed_ if res else _S_passed_
        tests = {}
        results = {_TestsKey_:tests}

        for f in json_files:
            jr = json.load( file(json_output_dir+'/'+f) )
            for k in jr:
                key = k[ len('src/') : ].replace('/', '_')
                if _StateKey_ not in jr[k]  or  jr[k][_StateKey_] != _S_passed_:
                    tests[key] = jr[k]
                    state = _S_failed_
                else:
                    output += 'Test {} is passed!\n'.format(key)

    else:
        state, results = _S_script_failed_, {}
        output += 'ERROR: serialization_validator tool run but i could not find json files... Terminating with script failure!'

    return {_StateKey_ : state,  _ResultsKey_ : results,  _LogKey_ : output }



def run(test, rosetta_dir, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    if test == 'serialization': return run_serialization_test(rosetta_dir, working_dir, platform, config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    else: raise BenchmarkError('Build script does not support TestSuite-like run!')
