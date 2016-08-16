#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   tests/PyRosetta4.py
## @brief  PyRosetta binding self tests
## @author Sergey Lyskov

import os, os.path, json, commands, shutil, distutils.dir_util

import imp
imp.load_source(__name__, '/'.join(__file__.split('/')[:-1]) +  '/__init__.py')  # A bit of Python magic here, what we trying to say is this: from __init__ import *, but init is calculated from file location

_api_version_ = '1.0'


def run_build_test(rosetta_dir, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    memory = config['memory'];  jobs = config['cpu_count']

    TR = Tracer(verbose)

    TR('Running PyRosetta build test: at working_dir={working_dir!r} with rosetta_dir={rosetta_dir}, platform={platform}, jobs={jobs}, memory={memory}GB, hpc_driver={hpc_driver}...'.format( **vars() ) )

    res, output, build_command_line, pyrosetta_path = build_pyrosetta(rosetta_dir, platform, jobs, config, mode='MinSizeRel', debug=debug)

    for f in os.listdir(pyrosetta_path + '/source'):
        if os.path.islink(pyrosetta_path + '/source/' + f): os.remove(pyrosetta_path + '/source/' + f)
    distutils.dir_util.copy_tree(pyrosetta_path + '/source', working_dir + '/source', update=False)

    res_code = _S_failed_ if res else _S_passed_
    if not res: output = '...\n'+'\n'.join( output.split('\n')[-32:] )  # truncating log for passed builds.
    output = 'Running: {}\n'.format(build_command_line) + output + '\nPyRosetta build path: ' + pyrosetta_path + '\n'  # Making sure that exact command line used is stored
    r = {_StateKey_ : res_code,  _ResultsKey_ : {},  _LogKey_ : output }
    # makeing sure that results could be serialize in to json, but ommiting logs because they could take too much space
    json.dump({_ResultsKey_:r[_ResultsKey_], _StateKey_:r[_StateKey_]}, file(working_dir+'/output.json', 'w'), sort_keys=True, indent=2)

    return r


def run_unit_tests(rosetta_dir, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    memory = config['memory'];  jobs = config['cpu_count']
    if platform['os'] != 'windows': jobs = jobs if memory/jobs >= PyRosetta_unix_memory_requirement_per_cpu else max(1, int(memory/PyRosetta_unix_memory_requirement_per_cpu) )  # PyRosetta require at least X Gb per memory per thread

    TR = Tracer(verbose)

    TR('Running PyRosetta unit tests: at working_dir={working_dir!r} with rosetta_dir={rosetta_dir}, platform={platform}, jobs={jobs}, memory={memory}GB, hpc_driver={hpc_driver}...'.format( **vars() ) )

    res, output, build_command_line, pyrosetta_path = build_pyrosetta(rosetta_dir, platform, jobs, config, mode='MinSizeRel', debug=debug)

    for f in os.listdir(pyrosetta_path + '/source'):
        if os.path.islink(pyrosetta_path + '/source/' + f): os.remove(pyrosetta_path + '/source/' + f)
    distutils.dir_util.copy_tree(pyrosetta_path + '/source', working_dir + '/source', update=False)

    file(working_dir+'/build-log.txt', 'w').write(output)

    if res:
        res_code = _S_build_failed_
        results = {_StateKey_ : res_code,  _ResultsKey_ : {},  _LogKey_ : output }
        json.dump({_ResultsKey_:results[_ResultsKey_], _StateKey_:results[_StateKey_]}, file(working_dir+'/output.json', 'w'), sort_keys=True, indent=2)

    else:

        distr_file_list = os.listdir(pyrosetta_path+'/build')

        #gui_flag = '--enable-gui' if platform['os'] == 'mac' else ''
        gui_flag = ''
        if not res: res, output = execute('Running PyRosetta tests...', 'cd {pyrosetta_path}/build && {python} self-test.py {gui_flag} -j{jobs}'.format(pyrosetta_path=pyrosetta_path, python=platform['python'], jobs=jobs, gui_flag=gui_flag), return_='tuple')

        json_file = pyrosetta_path + '/build/.test.output/.test.results.json'
        results = json.load( file(json_file) )

        execute('Deleting PyRosetta tests output...', 'cd {pyrosetta_path}/build && {python} self-test.py --delete-tests-output'.format(pyrosetta_path=pyrosetta_path, python=platform['python']), return_='tuple')
        extra_files = [f for f in os.listdir(pyrosetta_path+'/build') if f not in distr_file_list]  # not f.startswith('.test.')  and
        if extra_files:
            results['results']['tests']['self-test'] = dict(state='failed', log='self-test.py scripts failed to delete files: ' + ' '.join(extra_files))
            results[_StateKey_] = 'failed'

        if not res: output = '...\n'+'\n'.join( output.split('\n')[-32:] )  # truncating log for passed builds.
        output = 'Running: {}\n'.format(build_command_line) + output  # Making sure that exact command line used is stored

        #r = {_StateKey_ : res_code,  _ResultsKey_ : {},  _LogKey_ : output }
        results[_LogKey_] = output

        # makeing sure that results could be serialize in to json, but ommiting logs because they could take too much space
        json.dump({_ResultsKey_:results[_ResultsKey_], _StateKey_:results[_StateKey_]}, file(working_dir+'/output.json', 'w'), sort_keys=True, indent=2)

    return results




def run(test, rosetta_dir, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    ''' Run single test.
        Platform is a dict-like object, mandatory fields: {os='Mac', compiler='gcc'}
    '''
    if   test =='build': return run_build_test(rosetta_dir, working_dir, platform, config=config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    elif test =='unit':  return run_unit_tests(rosetta_dir, working_dir, platform, config=config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    else: raise BenchmarkError('Unknow PyRosetta test: {}!'.format(test))
