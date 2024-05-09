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

import sys, os, os.path, json, shutil
import codecs

# if sys.version_info[:2] < (3, 10):
#     from setuptools import distutils.dir_util
# else:
#     import distutils.dir_util

try: from setuptools.distutils import dir_util as dir_util_module
except ModuleNotFoundError: from distutils import dir_util as dir_util_module

# A bit of Python magic here, what we trying to say is this: from __init__ import *, but init is calculated from file location
import importlib.util, sys
importlib.util.spec_from_file_location(__name__, '/'.join(__file__.split('/')[:-1]) +  '/__init__.py').loader.exec_module(sys.modules[__name__])

_api_version_ = '1.1'


def run_build_test(rosetta_dir, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    memory = config['memory'];  jobs = config['cpu_count']

    TR = Tracer(verbose)

    TR('Running PyRosetta build test: at working_dir={working_dir!r} with rosetta_dir={rosetta_dir}, platform={platform}, jobs={jobs}, memory={memory}GB, hpc_driver={hpc_driver}...'.format( **vars() ) )

    result = build_pyrosetta(rosetta_dir, platform, jobs, config, mode='MinSizeRel', skip_compile=debug)

    for f in os.listdir(result.pyrosetta_path + '/source'):
        if os.path.islink(result.pyrosetta_path + '/source/' + f): os.remove(result.pyrosetta_path + '/source/' + f)
    dir_util_module.copy_tree(result.pyrosetta_path + '/source', working_dir + '/source', update=False)

    res_code = _S_failed_ if result.exitcode else _S_passed_
    if not result.exitcode: result.output = '...\n'+'\n'.join( result.output.split('\n')[-32:] )  # truncating log for passed builds.
    result.output = 'Running: {}\n'.format(result.command_line) + result.output + '\nPyRosetta build path: ' + result.pyrosetta_path + '\n'  # Making sure that exact command line used is stored
    r = {_StateKey_ : res_code,  _ResultsKey_ : {},  _LogKey_ : result.output }
    # makeing sure that results could be serialize in to json, but ommiting logs because they could take too much space
    with open(working_dir+'/output.json', 'w') as f: json.dump({_ResultsKey_:r[_ResultsKey_], _StateKey_:r[_StateKey_]}, f, sort_keys=True, indent=2)

    return r


def run_unit_tests(rosetta_dir, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    memory = config['memory'];  jobs = config['cpu_count']
    if platform['os'] != 'windows': jobs = jobs if memory/jobs >= PyRosetta_unix_memory_requirement_per_cpu else max(1, int(memory/PyRosetta_unix_memory_requirement_per_cpu) )  # PyRosetta require at least X Gb per memory per thread

    TR = Tracer(verbose)

    TR('Running PyRosetta unit tests: at working_dir={working_dir!r} with rosetta_dir={rosetta_dir}, platform={platform}, jobs={jobs}, memory={memory}GB, hpc_driver={hpc_driver}...'.format( **vars() ) )

    result = build_pyrosetta(rosetta_dir, platform, jobs, config, mode='MinSizeRel', skip_compile=config.get('skip_compile') )

    for f in os.listdir(result.pyrosetta_path + '/source'):
        if os.path.islink(result.pyrosetta_path + '/source/' + f): os.remove(result.pyrosetta_path + '/source/' + f)
    dir_util_module.copy_tree(result.pyrosetta_path + '/source', working_dir + '/source', update=False)

    codecs.open(working_dir+'/build-log.txt', 'w', encoding='utf-8', errors='backslashreplace').write(result.output)

    if result.exitcode:
        res_code = _S_build_failed_
        results = {_StateKey_ : res_code,  _ResultsKey_ : {},  _LogKey_ : result.output }
        with open(working_dir+'/output.json', 'w') as f: json.dump({_ResultsKey_:results[_ResultsKey_], _StateKey_:results[_StateKey_]}, f, sort_keys=True, indent=2)

    else:

        distr_file_list = os.listdir(result.pyrosetta_path+'/build')

        packages = ' '.join( get_required_pyrosetta_python_packages_for_testing(platform) ).replace('>', '=').replace('<', '=')

        python_virtual_environment = setup_persistent_python_virtual_environment(result.python_environment, packages)

        additional_flags = ' --timeout 4096' if platform['os'].startswith('aarch64') else ''

        #gui_flag = '--enable-gui' if platform['os'] == 'mac' else ''
        gui_flag, res, output = '', result.exitcode, result.output
        command_line = f'{python_virtual_environment.activate} && cd {result.pyrosetta_path}/build && {python_virtual_environment.python} {rosetta_dir}/source/test/timelimit.py 128 {python_virtual_environment.python} self-test.py {gui_flag} -j{jobs}{additional_flags}'
        output += '\nRunning PyRosetta tests: ' + command_line + '\n'

        res, o = execute('Running PyRosetta tests...', command_line, return_='tuple')
        output += o

        if res:
            results = {_StateKey_ : _S_script_failed_,  _ResultsKey_ : {},  _LogKey_ : f'{output}\n\nPyRosetta self-test.py script terminated with non-zero exit code, terminating with script failure!\n' }

        else:
            json_file = result.pyrosetta_path + '/build/.test.output/.test.results.json'
            with open(json_file) as f: results = json.load(f)

            execute('Deleting PyRosetta tests output...', 'cd {pyrosetta_path}/build && unset PYTHONPATH && unset __PYVENV_LAUNCHER__ && {python} self-test.py --delete-tests-output'.format(pyrosetta_path=result.pyrosetta_path, python=result.python), return_='tuple')
            extra_files = [f for f in os.listdir(result.pyrosetta_path+'/build') if f not in distr_file_list]  # not f.startswith('.test.')  and
            if extra_files:
                results['results']['tests']['self-test'] = dict(state='failed', log='self-test.py scripts failed to delete files: ' + ' '.join(extra_files))
                results[_StateKey_] = 'failed'

            if results[_StateKey_] == _S_passed_: output = '...\n'+'\n'.join( output.split('\n')[-32:] )  # truncating log for passed builds.
            output = 'Running: {}\n'.format(result.command_line) + output  # Making sure that exact command line used is stored

            #r = {_StateKey_ : res_code,  _ResultsKey_ : {},  _LogKey_ : output }
            results[_LogKey_] = output

            # makeing sure that results could be serialize in to json, but ommiting logs because they could take too much space
            with open(working_dir+'/output.json', 'w') as f: json.dump({_ResultsKey_:results[_ResultsKey_], _StateKey_:results[_StateKey_]}, f, sort_keys=True, indent=2)

    return results



def run_notebook_tests(rosetta_dir, working_dir, platform, config, hpc_driver, verbose, debug):
    memory = config['memory'];  jobs = config['cpu_count'];  skip_compile = config.get('skip_compile', False)
    TR = Tracer(verbose)

    packages = 'ipython pyrosettacolabsetup nbconvert'  # base packages for test
    packages += ' matplotlib biopython blosc dask distributed jupyter numpy pandas py3Dmol scipy traitlets graphviz seaborn dask-jobqueue attrs billiard GitPython fsspec' # extra packages for various notebooks
    #packages += ' traitlets==5.3.0 ipywidgets==8.0.1 jupyterlab-widgets==3.0.2 prompt-toolkit==3.0.30 widgetsnbextension==4.0.2'

    P = build_and_install_pyrosetta(working_dir, rosetta_dir, platform, jobs, config, mode='MinSizeRel', packages=packages, skip_compile=skip_compile)

    codecs.open(working_dir+'/build-log.txt', 'w', encoding='utf-8', errors='backslashreplace').write(P.output)

    if P.exitcode:
        results = {_StateKey_ : _S_build_failed_,  _ResultsKey_ : {},  _LogKey_ : P.output }

    else:
        notebooks_source_prefix = f'{rosetta_dir}/PyRosetta.notebooks/notebooks'
        notebooks_path = f'{working_dir}/notebooks'
        shutil.copytree(notebooks_source_prefix, notebooks_path)

        notebooks = [ f[:-len('.ipynb')] for f in os.listdir(notebooks_path) if f.endswith('.ipynb') and f not in ('index.ipynb', 'toc.ipynb') ]
        TR(f'notebooks: {notebooks}')

        jobs = {}
        for n in notebooks:
            command_line = f'cd {notebooks_path} && {P.python_virtual_environment.python} -m nbconvert --to python {n}.ipynb && export DEBUG="DEBUG" && {rosetta_dir}/source/test/timelimit.py 12 {P.python_virtual_environment.bin}/ipython --HistoryManager.enabled=False {n}.py'
            jobs[n] = command_line

        notebook_test_results = parallel_execute('notebook_tests', jobs, rosetta_dir, working_dir, config['cpu_count'], time=60)

        sub_tests = {}
        test_state = False
        for n, d in notebook_test_results.items():
            res = d['result']
            output = d['output']

            output = jobs[n] + '\n' + output
            sub_tests[n] = {_StateKey_ : _S_failed_ if res else _S_passed_, _LogKey_ : output }
            test_state |= res
            with open(f'{notebooks_path}/{n}.output', 'w') as f: f.write(output)


        # for n in notebooks:
        #     command_line = f'cd {notebooks_path} && {P.python_virtual_environment.python} -m nbconvert --to script {n}.ipynb && {rosetta_dir}/source/test/timelimit.py 8 {P.python_virtual_environment.python} {n}.py'
        #     res, output = execute(f'Running converting and running notebook {n}...', command_line, return_='tuple', add_message_and_command_line_to_output=True)
        #     sub_tests[n] = {_StateKey_ : _S_failed_ if res else _S_passed_, _LogKey_ : output }
        #     test_state |= res
        #     with open(f'{notebooks_path}/{n}.output', 'w') as f: f.write(output)

        results = {_StateKey_ : _S_failed_ if test_state else _S_passed_, _ResultsKey_: {_TestsKey_ : sub_tests}, _LogKey_ : P.output }

    if not config['emulation'] and os.path.isdir(P.python_virtual_environment.root): shutil.rmtree(P.python_virtual_environment.root)

    with open(working_dir+'/output.json', 'w') as f: json.dump(results, f, sort_keys=True, indent=2)

    return results


def run(test, repository_root, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    ''' Run single test.
        Platform is a dict-like object, mandatory fields: {os='Mac', compiler='gcc'}
    '''
    if   test =='build': return run_build_test(repository_root, working_dir, platform, config=config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    elif test =='unit':  return run_unit_tests(repository_root, working_dir, platform, config=config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    elif test =='notebook':  return run_notebook_tests(repository_root, working_dir, platform, config=config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    else: raise BenchmarkError('Unknow PyRosetta test: {}!'.format(test))
