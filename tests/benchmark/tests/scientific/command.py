#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   command.py
## @brief  scientific/command.py
## Python script for running multi-step scientific tests
## @author Sergey Lyskov

import os, json

import imp
imp.load_source(__name__, '/'.join(__file__.split('/')[:-1]) +  '/../__init__.py')  # A bit of Python magic here, what we trying to say is this: from ../__init__ import *, but init path is calculated relatively to this location


_api_version_ = '1.0'  # api version


def symlink(source, dest):
    ''' Similar to os.symlink but if dest is alread exisist and if type of *source and *dest does not match or if link points to different location : remove dest first
    '''
    relative_source = os.path.relpath( os.path.abspath(source), os.path.dirname(dest) )

    if os.path.islink(dest):
        if os.readlink(dest) == relative_source: return
        else: os.unlink(dest)
    elif os.path.isdir(dest): shutil.rmtree(dest)
    elif os.path.isfile(dest): os.remove(dest)

    os.symlink(relative_source, dest)


def symlink_tree(source, dest):
    ''' Similar to symlink(...) above but recursivly recreate dir tree at dest and symlink all files
    '''
    source = os.path.abspath(source)
    for dir_name, dirs, files in os.walk(source):
        prefix = dir_name[len(source):] + '/'

        # for d in dirs:
        #     dst = dest + prefix + d
        #     if not os.path.isdir(dst): os.makedirs(dst)

        for f in files: symlink(source + prefix + f, dest + prefix + f)


def run_multi_step_test(test, rosetta_dir, working_dir, platform, config, hpc_driver, verbose, debug):
    test_source_dir = f'{rosetta_dir}/tests/scientific/tests/{test}'

    symlink_tree(test_source_dir, working_dir)

    # os.mkdir( f'{working_dir}/benchmark' )
    # symlink(rosetta_dir + '/tests/benchmark/__init__.py', working_dir + '/benchmark/__init__.py')
    # symlink(rosetta_dir + '/tests/benchmark/hpc_drivers', working_dir + '/benchmark/hpc_drivers')

    multi_step_config = dict(config, test=test, rosetta_dir=rosetta_dir, working_dir=working_dir, platform=platform, verbose=verbose, debug=debug)

    #print('multi_step_config:', multi_step_config)
    with open(f'{working_dir}/{_multi_step_config_}', 'w') as f: json.dump(multi_step_config, f, sort_keys=True, indent=2)

    python = sys.executable
    scripts = sorted( f for f in os.listdir(working_dir) if f[0].isdigit() and f.endswith('.py') )
    for script in scripts:
        #print(script)
        res, output = execute(f'Running {script}...', f'cd {working_dir} && {python} {script}', return_=tuple, add_message_and_command_line_to_output=True)

        if res: return { _StateKey_ : _S_script_failed_,  _ResultsKey_ : {},
                         _LogKey_   : f'run_multi_step_test for {test} failed while running {script}...Aborting!\n\n{output}\n' }

        if os.path.isfile(f'{working_dir}/{_multi_step_error_}'):
            with open(f'{working_dir}/{_multi_step_error_}') as f: return json.load(f)

    with open(f'{working_dir}/{_multi_step_result_}') as f: return json.load(f)


def run(test, rosetta_dir, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    if test == "cartesian_relax": return run_multi_step_test(test, rosetta_dir, working_dir, platform, config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    if test == "fast_relax": return run_multi_step_test(test, rosetta_dir, working_dir, platform, config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    if test == "fast_relax_5iter": return run_multi_step_test(test, rosetta_dir, working_dir, platform, config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    if test == "enzyme_design": return run_multi_step_test(test, rosetta_dir, working_dir, platform, config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)

    else: raise BenchmarkError(f'Unknown scripts test: {test}!')
