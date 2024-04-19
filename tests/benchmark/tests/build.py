#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   build.py
## @brief  Rosetta/PyRosetta build tests
## @author Sergey Lyskov

import os, json, shutil
import codecs

# A bit of Python magic here, what we trying to say is this: from __init__ import *, but init is calculated from file location
import importlib.util, sys
importlib.util.spec_from_file_location(__name__, '/'.join(__file__.split('/')[:-1]) +  '/__init__.py').loader.exec_module(sys.modules[__name__])

_api_version_ = '1.1'

# _TestSuite_ = False  # Set to True for TestSuite-like tests (Unit, Integration, Sfxn_fingerprint) and False other wise


tests = dict(
    debug     = NT(command='{python} ./scons.py bin cxx={compiler} extras={extras} -j{jobs}', incremental=True),
    release   = NT(command='{python} ./scons.py bin cxx={compiler} extras={extras} mode=release -j{jobs}', incremental=True),
    static    = NT(command='{python} ./scons.py bin cxx={compiler} extras={extras} mode=release -j{jobs}', incremental=True),

    ninja_debug    = NT(command='./ninja_build.py debug    -v -remake -j{jobs}', incremental=True),
    ninja_release  = NT(command='./ninja_build.py release  -v -remake -j{jobs}', incremental=True),
    ninja_graphics = NT(command='./ninja_build.py graphics -v -remake -j{jobs}', incremental=True),

    #PyRosetta = NT(command='BuildPyRosetta.sh -u --monolith -j{jobs}', incremental=True),

    header    = NT(command='./scons.py unit_test_platform_only ; cd src && python3 ./../../tools/python_cc_reader/test_all_headers_compile_w_fork.py -n {jobs}', incremental=False),
    levels    = NT(command= ' && '.join(PRE_COMPILE_SETUP_SCRIPTS) + ' && {activate} && cd src && python ./../../tools/python_cc_reader/library_levels.py', incremental=False),

    ui  = NT(command='cd src/ui && {python} update_ui_project.py && cd ../../build && mkdir -p ui.{platform_suffix}.debug && cd ui.{platform_suffix}.debug && {qmake} -r ../qt/qt.pro {qt_extras}&& make -j{jobs}', incremental=True),

    xcode = NT(command='cd xcode && {python} make_project.py all && xcodebuild -scheme rosetta_scripts -configuration Debug -jobs {jobs} build', incremental=True),  #  -alltargets

    #py3_debug = NT(command='python3 ./scons.py bin cxx={compiler} extras={extras} -j{jobs}', incremental=True),
)

# def set_up():
#     pass
# def tear_down():
#     pass
# def rollover():
#     pass

# def get_tests():
#     return tests.keys()


def run_test(test, rosetta_dir, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    ''' Run single test.
        Platform is a dict-like object, mandatory fields: {os='Mac', compiler='gcc'}
    '''
    memory = config['memory']
    jobs = config['cpu_count']

    if platform['os'] != 'windows'  and  test.startswith('PyRosetta'): jobs = jobs if memory/jobs >= 1.0 else max(1, int(memory) )  # PyRosetta builds require at least 1Gb per memory per thread

    TR = Tracer(verbose)

    TR('Running test: "{test}" at working_dir={working_dir!r} with rosetta_dir={rosetta_dir}, platform={platform}, jobs={jobs}, memory={memory}GB, hpc_driver={hpc_driver}...'.format( **vars() ) )

    compiler = platform['compiler']
    extras   = ','.join(platform['extras'])
    platform_suffix = platform_to_pretty_string(platform)

    skip_compile = config.get('skip_compile', False) # Don't skip the actual build we're testing, just the compilation/installation of other things

    activate = 'ACTIVATE-WAS-NOT-SET'

    if skip_compile:
        python = sys.executable

    else:
        python_environment = local_python_install(platform, config)
        python = python_environment.python

        python_extra_packages = dict(
            levels='toolz==0.10',
        ).get(test)

        if python_extra_packages:
            python_virtual_environment = setup_persistent_python_virtual_environment(python_environment, python_extra_packages)
            activate = python_virtual_environment.activate

    if 'qmake' in config:
        qmake = config['qmake']
        qt_extras = '-spec linux-clang ' if (compiler == 'clang' and platform['os'] == 'linux') else ''
    else:
        qmake = 'qmake-path-is-not-defined-in-daemon-config'
        gt_extras = 'none-because-qmake-path-was-not-defined-in-daemon-config'

    command_line = tests[test].command.format( **vars() )

    if debug: res, output = 0, 'build.py: debug is enabled, skippig build phase...\n'
    else: res, output = execute('Compiling...', 'cd {}/source && {}'.format(rosetta_dir, command_line), return_='tuple')

    # re-running builds in case we got error - so we can get nice error message
    if res and tests[test].incremental:  res, output = execute('Compiling...', 'cd {}/source && {}'.format(rosetta_dir, tests[test].command.format( **vars() ) ), return_='tuple')

    #codecs.open(working_dir+'/build-log.txt', 'w', encoding='utf-8', errors='replace').write( u'Running: {}\n{}\n'.format(command_line, output) )
    codecs.open(working_dir+'/build-log.txt', 'w', encoding='utf-8', errors='backslashreplace').write( 'Running: {}\n{}\n'.format(command_line, output) )

    res_code = _S_failed_ if res else _S_passed_

    if not res and len(output) > 64*1024: output = '...truncated...\n'+'\n'.join( output.split('\n')[-32:] )  # truncating log for passed builds.

    # Moved to benchmark.py
    # if len(output) > 1024*1024*1:  # truncating logs if they too large (more then 1Mb in size)...
    #     lines = output.split('\n')
    #     output = '\n'.join( lines[:32] + ['...truncated...'] + lines[-32:] )

    output = 'Running: {}\n'.format(command_line) + output  # Making sure that exact command line used is stored

    r = {_StateKey_ : res_code,  _ResultsKey_ : {},  _LogKey_ : output }

    # makeing sure that results could be serialize in to json, but ommiting logs because they could take too much space
    with open(working_dir+'/output.json', 'w') as f: json.dump({_ResultsKey_:r[_ResultsKey_], _StateKey_:r[_StateKey_]}, f, sort_keys=True, indent=2)

    return r



def run_test_on_fresh_clone(test, rosetta_dir, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    ''' create a fresh clone and run specified test on it
    '''
    clean_main = 'clean_main'

    execute('Cloning main...', 'cd {working_dir} && git clone {rosetta_dir} {clean_main}'.format(**vars()) )


    convert_submodule_urls_from_ssh_to_https(f'{working_dir}/clean_main')

    res = run_test(test, working_dir + '/' + clean_main, working_dir, platform, config, hpc_driver, verbose, debug)

    shutil.rmtree( working_dir + '/' + clean_main )  # removing local clone to avoid storing main/ files as test-result-files

    return res


def run_test_suite(rosetta_dir, working_dir, platform, jobs=1, hpc_driver=None, verbose=False, debug=False):
    raise BenchmarkError('Build script does not support TestSuite-like run!')


def run(test, repository_root, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    if test == "cppcheck":
        # Backwards compatibility shim
        from .code_quality import run as code_quality_run
        return code_quality_run(test, repository_root, working_dir, platform, config, hpc_driver, verbose, debug)
    elif test and test.startswith('clean.'): return run_test_on_fresh_clone( test[ len('clean.') : ], repository_root, working_dir, platform, config=config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    elif test: return run_test(test, repository_root, working_dir, platform, config=config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    else: return run_test_suite(repository_root, working_dir, platform, jobs=jobs, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
