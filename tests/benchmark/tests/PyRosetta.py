#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   tests/PyRosetta.py
## @brief  PyRosetta binding self tests
## @author Sergey Lyskov

import os, os.path, json, commands

import imp
imp.load_source(__name__, '/'.join(__file__.split('/')[:-1]) +  '/__init__.py')  # A bit of Python magic here, what we trying to say is this: from __init__ import *, but init is calculated from file location

_api_version_ = '1.0'  # api version


def run_test(test, rosetta_dir, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    ''' Run single test.
        Platform is a dict-like object, mandatory fields: {os='Mac', compiler='gcc'}
    '''
    memory = config['memory']
    jobs = config['cpu_count']

    jobs = jobs if memory/jobs >= 1.0 else max(1, int(memory) )  # PyRosetta builds require at least 1Gb per memory per thread

    TR = Tracer(verbose)

    TR('Running test: "{test}" at working_dir={working_dir!r} with rosetta_dir={rosetta_dir}, platform={platform}, jobs={jobs}, memory={memory}GB, hpc_driver={hpc_driver}...'.format( **vars() ) )

    compiler = platform['compiler']
    extras   = ','.join(platform['extras'])

    command_line = 'cd {rosetta_dir}/source && BuildPyRosetta.sh -u --monolith -j{jobs}'.format(rosetta_dir=rosetta_dir, compiler=compiler, jobs=jobs, extras=extras)

    buildings_path_output = execute('Getting buindings build path...', command_line + ' --print-build-path', return_='tuple')
    buildings_path = buildings_path_output[1].split()[-1]
    if not (buildings_path  and  os.path.isdir(buildings_path) ): raise BenchmarkError('Could not retrieve valid PyRosetta bindings binary path!\nCommand line:{}\nResult:{}\n'.format(command_line, buildings_path_output))

    print buildings_path
    TR('Bindings build path is:{}'.format(buildings_path))

    res, output = execute('Compiling...', 'cd {}/source && {}'.format(rosetta_dir, command_line), return_='tuple')

    if res:  res, output = execute('Compiling...', 'cd {}/source && {}'.format(rosetta_dir, command_line.format(compiler=compiler, jobs=1, extras=extras)), return_='tuple')

    file(working_dir+'/build-log.txt', 'w').write(output)

    if not res: res, output = execute('Running PyRosetta tests...', 'cd {buildings_path} && python TestBindings.py'.format(buildings_path=buildings_path), return_='tuple')

    res_code = _S_failed_ if res else _S_finished_
    if not res: output = '...\n'+'\n'.join( output.split('\n')[-32:] )  # truncating log for passed builds.
    output = 'Running: {}\n'.format(command_line) + output  # Making sure that exact command line used is stored
    r = {_StateKey_ : res_code,  _ResultsKey_ : {},  _LogKey_ : output }
    # makeing sure that results could be serialize in to json, but ommiting logs because they could take too much space
    json.dump({_ResultsKey_:r[_ResultsKey_], _StateKey_:r[_StateKey_]}, file(working_dir+'/output.json', 'w'), sort_keys=True, indent=2)

    return r


def run(test, rosetta_dir, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    if not test: return run_test('', rosetta_dir, working_dir, platform, config=config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    else: raise BenchmarkError('Unknow PyRosetta test: {}!'.format(test))
