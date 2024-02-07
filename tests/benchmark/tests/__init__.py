#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   tests/__init__.py
## @brief  Common constats and types for all test types
## @author Sergey Lyskov

import os, time, sys, shutil, codecs, urllib.request, imp, subprocess, json, hashlib  # urllib.error, urllib.parse,
import platform as  platform_module
import types as types_module

# ⚔ do not change wording below, it have to stay in sync with upstream (up to benchmark-model).
# Copied from benchmark-model, standard state code's for tests results.

__all__ = ['execute',
           '_S_Values_', '_S_draft_', '_S_queued_', '_S_running_', '_S_passed_', '_S_failed_', '_S_build_failed_', '_S_script_failed_',
           '_StateKey_', '_ResultsKey_', '_LogKey_', '_DescriptionKey_', '_TestsKey_',
           '_multi_step_config_', '_multi_step_error_', '_multi_step_result_',
           'to_bytes',
]

_S_draft_                 = 'draft'
_S_queued_                = 'queued'
_S_running_               = 'running'
_S_passed_                = 'passed'
_S_failed_                = 'failed'
_S_build_failed_          = 'build failed'
_S_script_failed_         = 'script failed'
_S_queued_for_comparison_ = 'queued for comparison'

_S_Values_ = [_S_draft_, _S_queued_, _S_running_, _S_passed_, _S_failed_, _S_build_failed_, _S_script_failed_, _S_queued_for_comparison_]

_IgnoreKey_      = 'ignore'
_StateKey_       = 'state'
_ResultsKey_     = 'results'
_LogKey_         = 'log'
_DescriptionKey_ = 'description'
_TestsKey_       = 'tests'
_SummaryKey_     = 'summary'
_FailedKey_      = 'failed'
_TotalKey_       = 'total'
_PlotsKey_       = 'plots'
_FailedTestsKey_ = 'failed_tests'
_HtmlKey_        = 'html'

# file names for multi-step test files
_multi_step_config_ = 'config.json'
_multi_step_error_  = 'error.json'
_multi_step_result_ = 'result.json'

PyRosetta_unix_memory_requirement_per_cpu = 3.9  # Memory per sub-process in Gb's
PyRosetta_unix_unit_test_memory_requirement_per_cpu = 3.0  # Memory per sub-process in Gb's for running PyRosetta unit tests

# Commands to run all the scripts needed for setting up Rosetta compiles. (Run from main/source directory)
PRE_COMPILE_SETUP_SCRIPTS = [ "./update_options.sh", "./update_submodules.sh", "./update_ResidueType_enum_files.sh", "python version.py" ]

DEFAULT_PYTHON_VERSION='3.9'

# Standard funtions and classes below ---------------------------------------------------------------------------------

class BenchmarkError(Exception):
    def __init__(self, value): self.value = value
    def __repr__(self): return self.value
    def __str__(self): return self.value


class NT:  # named tuple
    def __init__(self, **entries): self.__dict__.update(entries)
    def __repr__(self):
        r = 'NT: |'
        for i in dir(self):
            print(i)
            if not i.startswith('__') and i != '_as_dict' and not isinstance(getattr(self, i), types_module.MethodType): r += '%s --> %s, ' % (i, getattr(self, i))
        return r[:-2]+'|'

    @property
    def _as_dict(self):
        return { a: getattr(self, a) for a in dir(self) if not a.startswith('__') and a != '_as_dict' and not isinstance(getattr(self, a), types_module.MethodType)}


def Tracer(verbose=False):
    return print if verbose else lambda x: None


def to_unicode(b):
    ''' Conver bytes to string and handle the errors. If argument is already in string - do nothing
    '''
    #return b if type(b) == unicode else unicode(b, 'utf-8', errors='replace')
    return b if type(b) == str else str(b, 'utf-8', errors='backslashreplace')


def to_bytes(u):
    ''' Conver string to bytes and handle the errors. If argument is already of type bytes - do nothing
    '''
    return u if type(u) == bytes else u.encode('utf-8', errors='backslashreplace')


''' Python-2 version
def execute(message, commandline, return_=False, until_successes=False, terminate_on_failure=True, add_message_and_command_line_to_output=False):
    message, commandline = to_unicode(message), to_unicode(commandline)

    TR = Tracer()
    TR(message);  TR(commandline)
    while True:
        (res, output) = commands.getstatusoutput(commandline)
        # Subprocess results will always be a bytes-string.
        # Probably ASCII, but may have some Unicode characters.
        # A UTF-8 decode will probably get decent results 99% of the time
        # and the replace option will gracefully handle the rest.
        output = to_unicode(output)

        TR(output)

        if res and until_successes: pass  # Thats right - redability COUNT!
        else: break

        print( "Error while executing %s: %s\n" % (message, output) )
        print( "Sleeping 60s... then I will retry..." )
        time.sleep(60)

    if add_message_and_command_line_to_output: output = message + '\nCommand line: ' + commandline + '\n' + output

    if return_ == 'tuple': return(res, output)

    if res and terminate_on_failure:
        TR("\nEncounter error while executing: " + commandline)
        if return_==True: return res
        else:
            print("\nEncounter error while executing: " + commandline + '\n' + output)
            raise BenchmarkError("\nEncounter error while executing: " + commandline + '\n' + output)

    if return_ == 'output': return output
    else: return res
'''

def execute_through_subprocess(command_line):
    # exit_code, output = subprocess.getstatusoutput(command_line)

    # p = subprocess.Popen(command_line, bufsize=0, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # output, errors = p.communicate()
    # output = (output + errors).decode(encoding='utf-8', errors='backslashreplace')
    # exit_code = p.returncode

    # previous 'main' version based on subprocess module. Main issue that output of segfaults will not be captured since they generated by shell
    p = subprocess.Popen(command_line, bufsize=0, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output, errors = p.communicate()
    # output = output + errors # ← we redirected stderr into same pipe as stdcout so errors is None, - no need to concatenate
    output = output.decode(encoding='utf-8', errors='backslashreplace')
    exit_code = p.returncode

    return exit_code, output


def execute_through_pexpect(command_line):
    import pexpect

    child = pexpect.spawn('/bin/bash', ['-c', command_line])
    child.expect(pexpect.EOF)
    output = child.before.decode(encoding='utf-8', errors='backslashreplace')
    child.close()
    exit_code = child.signalstatus or child.exitstatus

    return exit_code, output


def execute_through_pty(command_line):
    import pty, select

    if sys.platform == "darwin":

        master, slave = pty.openpty()
        p = subprocess.Popen(command_line, shell=True, stdout=slave, stdin=slave,
                             stderr=subprocess.STDOUT, close_fds=True)

        buffer = []
        while True:
            try:
                if select.select([master], [], [], 0.2)[0]:  # has something to read
                    data = os.read(master, 1 << 22)
                    if data: buffer.append(data)

                elif (p.poll() is not None)  and  (not select.select([master], [], [], 0.2)[0] ): break  # process is finished and output buffer if fully read

            except OSError: break  # OSError will be raised when child process close PTY descriptior

        output = b''.join(buffer).decode(encoding='utf-8', errors='backslashreplace')

        os.close(master)
        os.close(slave)

        p.wait()
        exit_code = p.returncode

        '''
        buffer = []
        while True:
            if select.select([master], [], [], 0.2)[0]:  # has something to read
                data = os.read(master, 1 << 22)
                if data: buffer.append(data)
                # else: break  # # EOF - well, technically process _should_ be finished here...

            # elif time.sleep(1) or (p.poll() is not None): # process is finished (sleep here is intentional to trigger race condition, see solution for this on the next few lines)
            #     assert not select.select([master], [], [], 0.2)[0]  # should be nothing left to read...
            #     break

            elif (p.poll() is not None)  and  (not select.select([master], [], [], 0.2)[0] ): break  # process is finished and output buffer if fully read

        assert not select.select([master], [], [], 0.2)[0]  # should be nothing left to read...

        os.close(slave)
        os.close(master)

        output = b''.join(buffer).decode(encoding='utf-8', errors='backslashreplace')
        exit_code = p.returncode
        '''

    else:

        master, slave = pty.openpty()
        p = subprocess.Popen(command_line, shell=True, stdout=slave, stdin=slave,
                             stderr=subprocess.STDOUT, close_fds=True)

        os.close(slave)

        buffer = []
        while True:
            try:
                data = os.read(master, 1 << 22)
                if data: buffer.append(data)
            except OSError: break  # OSError will be raised when child process close PTY descriptior

        output = b''.join(buffer).decode(encoding='utf-8', errors='backslashreplace')

        os.close(master)

        p.wait()
        exit_code = p.returncode

    return exit_code, output



def execute(message, command_line, return_='status', until_successes=False, terminate_on_failure=True, silent=False, silence_output=False, silence_output_on_errors=False, add_message_and_command_line_to_output=False):
    if not silent: print(message);  print(command_line); sys.stdout.flush();
    while True:

        #exit_code, output = execute_through_subprocess(command_line)
        #exit_code, output = execute_through_pexpect(command_line)
        exit_code, output = execute_through_pty(command_line)

        if (exit_code  and  not silence_output_on_errors) or  not (silent or silence_output): print(output); sys.stdout.flush();

        if exit_code and until_successes: pass  # Thats right - redability COUNT!
        else: break

        print( "Error while executing {}: {}\n".format(message, output) )
        print("Sleeping 60s... then I will retry...")
        sys.stdout.flush();
        time.sleep(60)

    if add_message_and_command_line_to_output: output = message + '\nCommand line: ' + command_line + '\n' + output

    if return_ == 'tuple'  or  return_ == tuple: return(exit_code, output)

    if exit_code and terminate_on_failure:
        print("\nEncounter error while executing: " + command_line)
        if return_==True: return True
        else:
            print('\nEncounter error while executing: ' + command_line + '\n' + output);
            raise BenchmarkError('\nEncounter error while executing: ' + command_line + '\n' + output)

    if return_ == 'output': return output
    else: return exit_code


def parallel_execute(name, jobs, rosetta_dir, working_dir, cpu_count, time=16):
    ''' Execute command line in parallel on local host
        time specifies the upper limit for cpu-usage runtime (in minutes) for any one process in the parallel execution.

        jobs should be dict with following structure:
        {
            'job-string-id-1’: command_line-1,
            'job-string-id-2’: command_line-2,
            ...
        }

        return: dict with jobs-id's as keys and value as dict with 'output' and 'result' keys:
        {
            "job-string-id-1": {
                "output": "stdout + stdderr output of command_line-1",
                "result": <integer exit code for command_line-1>
            },
            "c2": {
                "output": "stdout + stdderr output of command_line-2",
                "result": <integer exit code for command_line-2>
            },
            ...
        }
    '''
    job_file_name = working_dir + '/' + name
    with open(job_file_name + '.json', 'w') as f: json.dump(jobs, f, sort_keys=True, indent=2) # JSON handles unicode internally
    if time is not None:
        allowed_time = int(time*60)
        ulimit_command = f'ulimit -t {allowed_time} && '
    else:
        ulimit_command = ''
    command = f'cd {working_dir} && ' + ulimit_command + f'{rosetta_dir}/tests/benchmark/util/parallel.py -j{cpu_count} {job_file_name}.json'
    execute("Running {} in parallel with {} CPU's...".format(name, cpu_count), command )

    with open(job_file_name+'.results.json') as f: return json.load(f)


def calculate_extension(platform, mode='release'):
    ''' Calculate and return extension for Rosetta executables that will match one generated by scons '''

    if platform['os'] in ['windows']:     os = 'windows'
    elif platform['os'] in ['mac', 'm1']: os = 'macos'
    else:                                 os = 'linux'

    extras = ''.join( sorted(platform['extras']) ) if platform['extras'] else 'default'  # platform and therfor extras should be comming from database and therefor it should be already sorted

    return extras + "." + os + platform['compiler'] + mode


def calculate_unique_prefix_path(platform, config):
    ''' calculate path for prefix location that is unique for this machine and OS
    '''
    hostname = os.uname()[1]
    return config['prefix'] + '/' + hostname + '/' + platform['os']


def platform_to_pretty_string(platform):
    ''' Take platform as json object and return normalized human-readable string '''
    return '{}.{}{}{}'.format(platform['os'], platform['compiler'], ('.'+'.'.join(platform['extras']) if 'extras' in platform  and  platform['extras'] else ''), ('.python'+platform['python'] if 'python' in platform else ''))

def setup_for_compile(rosetta_dir):
    execute('Updating options, ResidueTypes and version info...', f'cd {rosetta_dir}/source && ' + ' && '.join(PRE_COMPILE_SETUP_SCRIPTS) )

def build_rosetta(rosetta_dir, platform, config, mode='release', build_unit=False, verbose=False):
    ''' Compile Rosetta binaries on a given platform return (res, output, build_command_line) '''

    compiler = platform['compiler']
    extras   = ','.join(platform['extras'])
    jobs = config['cpu_count']
    skip_compile = config.get('skip_compile', False)

    if skip_compile:
        python = "<python_executable>"
    else:
        python = local_python_install(platform, config).python # Will install a local python

    # removing all symlinks from bin/ and then building binaries...
    build_command_line = f'find bin -type l ! -name ".*" -exec rm {{}} \\; ; {python} ./scons.py bin mode={mode} cxx={compiler} extras={extras} -j{jobs}'
    if build_unit:
        build_command_line += f' && {python} ./scons.py cxx={compiler} mode={mode} extras={extras} cat=test -j{jobs}'

    if skip_compile:
        build_command_line = "SKIPPED: " + build_command_line
        res, output = 0, 'build_rosetta: skip_compile is enabled, skipping build...\n'
    else:
        res, output = execute('Compiling...', 'cd {}/source && {}'.format(rosetta_dir, build_command_line), return_='tuple')

    return res, output, build_command_line



def get_required_pyrosetta_python_packages_for_testing(platform):
    ''' return list of Python packages that is required to run PyRosetta for given platform

        IMPORTANT: each package should have version specification in it by either using ==, >= or <=
        so we can have reproducible test enviroment

        IMPORTANT: there should be no spaces between package name and version number
    '''
    # not available in standard Conda channels:
    #    blosc==1.8.3         \
    #    py3Dmol>=0.8.0       \

    python_version = tuple( map(int, platform.get('python', DEFAULT_PYTHON_VERSION).split('.') ) )


    if python_version <= (3, 8):
        packages = 'attrs>=19.3.0 billiard>=3.6.3.0 cloudpickle>=1.4.1 dask>=2.16.0 dask-jobqueue>=0.7.0 distributed>=2.16.0 gitpython>=3.1.1 jupyter>=1.0.0 traitlets>=4.3.3 blosc>=1.8.3 numpy>=1.17.3 pandas>=0.25.2 scipy>=1.4.1 python-xz>=0.4.0'
        if python_version == (3, 7): packages = " ".join(map(lambda p: "cloudpickle<=0.7.0" if p.startswith("cloudpickle") else p, packages.split()))
        elif python_version == (3, 6): packages = " ".join(map(lambda p: "cloudpickle<=0.7.0" if p.startswith("cloudpickle") else p.split(">=")[0], packages.split()))

    elif python_version == (3, 9): packages = 'numpy>=1.20.2'
    elif python_version == (3, 10): packages = 'numpy>=1.22.3'
    elif python_version == (3, 11): packages = 'numpy>=1.23.5'
    elif python_version == (3, 12): packages = 'numpy>=1.26.0'
    else: packages = 'numpy>=1.23'

    if platform['os'] == 'mac' and python_version == (3, 7): packages = packages.replace('blosc>=1.8.3', 'blosc>=1.10.6')
    if platform['os'] == 'mac' and python_version == (3, 8): packages = packages.replace('blosc>=1.8.3', 'blosc>=1.10.6')

    packages = packages.split() if 'serialization' in platform['extras'] and platform.get('python', DEFAULT_PYTHON_VERSION)[:2] != '2.' else []

    if python_version >= (3, 7):
        for p in packages: assert '=' in p

    return packages



def get_required_pyrosetta_python_packages_for_release_package(platform, conda):
    ''' return list of Python packages that is required to run PyRosetta for given platform

        IMPORTANT: each package should have version specification in it by either using ==, >= or <=
        so we can have reproducible test enviroment

        IMPORTANT: there should be no spaces between package name and version number

        IMPORTANT: if PyRosetta build without serialization support or on Python-2 platform then we DO NOT REQUIRE any packages
    '''

    python_version = tuple( map(int, platform.get('python', DEFAULT_PYTHON_VERSION).split('.') ) )

    if python_version < (3, 9):
        packages = '\
        blosc>=1.8.3         \
        cloudpickle>=1.4.1   \
        dask>=2.16.0         \
        dask-jobqueue>=0.7.0 \
        distributed>=2.16.0  \
        jupyter>=1.0.0       \
        numpy>=1.17.3        \
        pandas>=0.25.2       \
        python-xz>=0.4.0     \
        scipy>=1.4.1         \
        traitlets>=4.3.3     \
        '
        if python_version == (3, 7): packages = " ".join(map(lambda p: "cloudpickle<=0.7.0" if p.startswith("cloudpickle") else p, packages.split()))
        elif python_version == (3, 6): packages = " ".join(map(lambda p: "cloudpickle<=0.7.0" if p.startswith("cloudpickle") else p.split(">=")[0], packages.split()))

    else:
        #packages = 'numpy>=1.19.2' if platform['os'] == 'mac' else 'numpy>=1.19.2'
        packages = 'numpy>=1.20.2'

    if conda:
        packages = packages.replace('blosc', 'python-blosc')
        packages = " ".join(map(lambda p: "xz" if p.startswith("python-xz") else p, packages.split()))

    packages = packages.split() if 'serialization' in platform['extras'] and platform.get('python', DEFAULT_PYTHON_VERSION)[:2] != '2.' else []

    if not conda:
        for p in packages: assert '=' in p

    return packages


def build_pyrosetta(rosetta_dir, platform, jobs, config, mode='MinSizeRel', options='', conda=None, verbose=False, skip_compile=False, version=None):
    ''' Compile Rosetta binaries on a given platform return NT(exitcode, output, build_command_line, pyrosetta_path, python) '''

    #binder = install_llvm_tool('binder', source_location='{}/source/src/python/PyRosetta/binder'.format(rosetta_dir), config=config)

    py_env = conda if conda else local_python_install(platform, config)

    #print(sysconfig.get_config_vars())
    #CONFINCLUDEPY

    extra = ' --python-include-dir={py_env.python_include_dir} --python-lib={py_env.python_lib_dir}'.format(**vars())
    # if platform['os'] == 'mac'  and  platform['python'].startswith('python3'):
    #     python_prefix = execute('Getting {} prefix path...'.format(platform['python']), '{}-config --prefix'.format(platform['python']), return_='output')
    #     extra += ' --python-include-dir={0}/include/python3.5m --python-lib={0}/lib/libpython3.5.dylib'.format(python_prefix)

    if 'cxx11thread'   in platform['extras']: extra += ' --multi-threaded'
    if 'serialization' in platform['extras']: extra += ' --serialization'

    if version: extra += " --version '{version}'".format(**vars())

    command_line = f'cd {rosetta_dir}/source/src/python/PyRosetta && {py_env.python} build.py -j{jobs} --compiler {platform["compiler"]} --type {mode}{extra} {options}'

    pyrosetta_path = execute('Getting PyRosetta build path...', command_line + ' --print-build-root', return_='output').split()[-1]

    if skip_compile:
        res, output = 0, '__init__.py:build_pyrosetta: skip_compile is enabled, skipping build...\n'
    else:
        res, output = execute('Building PyRosetta {}...'.format(mode), command_line, return_='tuple', add_message_and_command_line_to_output=True)

    return NT(exitcode=res, output=output, command_line=command_line, pyrosetta_path=pyrosetta_path, python=py_env.python, root=py_env.root, python_environment=py_env)



def build_and_install_pyrosetta(working_dir, rosetta_dir, platform, jobs, config, mode='MinSizeRel', options='', packages='', conda=None, verbose=False, skip_compile=False, version=None):
    ''' Compile PyRosetta on configured platform and install it into current virtual environment and return (res, output, build_command_line) '''

    python_environment = local_python_install(platform, config)
    python_virtual_environment_path = working_dir+'/.python_virtual_environment'
    python_virtual_environment = setup_python_virtual_environment(python_virtual_environment_path, python_environment, packages=packages)

    pyrosetta_package_path = python_environment.root + '/.pyrosetta_package'
    if os.path.isdir(pyrosetta_package_path): shutil.rmtree(pyrosetta_package_path)

    options = f'--create-package {pyrosetta_package_path}' + ( ' -sd' if skip_compile else '' )
    r = build_pyrosetta(rosetta_dir=rosetta_dir, platform=platform, jobs=config['cpu_count'], config=config, mode=mode, options=options, skip_compile=False)

    if not r.exitcode: execute('Installing PyRosetta...', f'cd {pyrosetta_package_path}/setup && {python_virtual_environment.python} setup.py install')

    return NT(exitcode = r.exitcode, output = r.output, command_line = r.command_line, python_virtual_environment = python_virtual_environment)



def install_llvm_tool_deprecated(name, source_location, platform, config, clean=True):
    ''' Install and update (if needed) custom LLVM tool at given prefix (from config).
        Return absolute path to executable on success and raise BenchmarkError exception on failure (do not catch this! if you really need 'normal' exit from this function on failure - refactor it instead)
    '''

    prefix = calculate_unique_prefix_path(platform, config)

    jobs = config['cpu_count']

    release = 'release_37'
    git_checkout = '( git checkout {0} && git reset --hard {0} )'.format(release) if clean else 'git checkout {}'.format(release)

    if not os.path.isdir(prefix): os.makedirs(prefix)

    if not os.path.isdir(prefix+'/llvm'): execute('Clonning llvm...', 'cd {} && git clone http://llvm.org/git/llvm.git llvm'.format(prefix) )
    execute('Checking out LLVM revision: {}...'.format(release), 'cd {prefix}/llvm && ( {git_checkout} || ( git fetch && {git_checkout} ) )'.format(prefix=prefix, git_checkout=git_checkout) )

    if not os.path.isdir(prefix+'/llvm/tools/clang'): execute('Clonning clang...', 'cd {}/llvm/tools && git clone http://llvm.org/git/clang.git clang'.format(prefix) )
    execute('Checking out Clang revision: {}...'.format(release), 'cd {prefix}/llvm/tools/clang && ( {git_checkout} || ( git fetch && {git_checkout} ) )'.format(prefix=prefix, git_checkout=git_checkout) )

    if not os.path.isdir(prefix+'/llvm/tools/clang/tools/extra'): execute('Clonning clang...', 'cd {}/llvm/tools/clang/tools && git clone http://llvm.org/git/clang-tools-extra.git extra'.format(prefix) )
    execute('Checking out Clang-tools revision: {}...'.format(release), 'cd {prefix}/llvm/tools/clang/tools/extra && ( {git_checkout} || ( git fetch && {git_checkout} ) )'.format(prefix=prefix, git_checkout=git_checkout) )

    tool_link_path = '{prefix}/llvm/tools/clang/tools/extra/{name}'.format(prefix=prefix, name=name)
    if os.path.islink(tool_link_path): os.unlink(tool_link_path)
    os.symlink(source_location, tool_link_path)

    cmake_lists = prefix + '/llvm/tools/clang/tools/extra/CMakeLists.txt'
    tool_build_line = 'add_subdirectory({})\n'.format(name)

    for line in codecs.open(cmake_lists, encoding='utf-8', errors='replace'):
        if line == tool_build_line: break
    else:
        with codecs.open(cmake_lists, 'w', encoding='utf-8', errors='replace') as f: f.write( codecs.open(cmake_lists, encoding='utf-8', errors='replace').read() + tool_build_line)

    build_dir = prefix+'/llvm/build-' + release
    if not os.path.isdir(build_dir): os.makedirs(build_dir)
    execute('Building tool: {}...'.format(name), 'cd {build_dir} && cmake -DCMAKE_BUILD_TYPE=Release .. && make -j{jobs}'.format(build_dir=build_dir, jobs=jobs) )

    # build_dir = prefix+'/llvm/build-ninja-' + release
    # if not os.path.isdir(build_dir): os.makedirs(build_dir)
    # execute('Building tool: {}...'.format(name), 'cd {build_dir} && cmake -DCMAKE_BUILD_TYPE=Release .. -G Ninja && ninja -j{jobs}'.format(build_dir=build_dir, jobs=jobs) )

    executable = build_dir + '/bin/' + name
    if not os.path.isfile(executable): raise BenchmarkError("\nEncounter error while running install_llvm_tool: Build is complete but executable {} is not there!!!".format(executable) )

    return executable


def install_llvm_tool(name, source_location, platform, config, clean=True):
    ''' Install and update (if needed) custom LLVM tool at given prefix (from config).
        Return absolute path to executable on success and terminate with error on failure
    '''
    debug = False

    prefix_root = calculate_unique_prefix_path(platform, config)

    if not os.path.isdir(prefix_root): os.makedirs(prefix_root)

    llvm_version='3.7.1'
    prefix = prefix_root + '/llvm-' + llvm_version

    build_dir = prefix+'/llvm-' + llvm_version + '.' + platform_to_pretty_string(platform) + ('.debug' if debug else '.release')

    res, output = execute('Getting source HEAD commit SHA1...', 'cd {} && git rev-parse HEAD'.format(source_location), return_='tuple', silent=True)
    if res: binder_head = 'unknown'
    else: binder_head = output.split('\n')[0]

    signature = dict(config = 'LLVM install by install_llvm_tool version: 1.1, HTTPS', binder = binder_head, llvm_version=llvm_version, compiler=platform['compiler'], gcc_install_prefix=None)
    signature_file_name = build_dir + '/.signature.json'

    disk_signature = dict(config = 'unknown', binder = 'unknown')
    if os.path.isfile(signature_file_name):
        with open(signature_file_name) as f: disk_signature = json.load(f)

    if signature == disk_signature:
        print('LLVM:{} + Binder install is detected at {}, skipping LLVM installation and Binder building procedures...\n'.format(llvm_version, build_dir))

    else:
        print('LLVM build detected, but config/binder version has changed, perfoming a clean rebuild...')
        if os.path.isdir(build_dir): shutil.rmtree(build_dir)

        def get_cmake_compiler_options():
            ''' Get cmake compiler flags from Options.compiler '''
            if platform['compiler'] == 'clang': return ' -DCMAKE_C_COMPILER=`which clang` -DCMAKE_CXX_COMPILER=`which clang++`'
            elif platform['compiler'] == 'gcc': return ' -DCMAKE_C_COMPILER=`which gcc` -DCMAKE_CXX_COMPILER=`which g++`'
            else: return ''

        clang_path = "{prefix}/tools/clang".format(**locals())

        if not os.path.isfile(prefix + '/CMakeLists.txt'): execute('Download LLVM source...', 'cd {prefix_root} && curl https://releases.llvm.org/{llvm_version}/llvm-{llvm_version}.src.tar.xz | tar -Jxom && mv llvm-{llvm_version}.src {prefix}'.format(**locals()) )

        if not os.path.isdir(clang_path): execute('Download Clang source...', 'cd {prefix_root} && curl https://releases.llvm.org/{llvm_version}/cfe-{llvm_version}.src.tar.xz | tar -Jxom && mv cfe-{llvm_version}.src {clang_path}'.format(**locals()) )

        if not os.path.isdir(prefix+'/tools/clang/tools/extra'): os.makedirs(prefix+'/tools/clang/tools/extra')


        # if signature['config'] != disk_signature['config']:
        #     print( 'LLVM build detected, but config version mismatched: was:"{}" current:"{}", perfoming a clean rebuild...'.format(disk_signature['config'], signature['config']) )
        #     if os.path.isdir(build_dir): shutil.rmtree(build_dir)
        # else: print( 'Binder build detected, but source version mismatched: was:{} current:{}, rebuilding...'.format(disk_signature['binder'], signature['binder']) )
        '''
        git_checkout = '( git checkout {0} && git reset --hard {0} )'.format(release) if clean else 'git checkout {}'.format(release)

        #if os.path.isdir(prefix) and  (not os.path.isdir(prefix+'/.git')): shutil.rmtree(prefix)  # removing old style checkoiut

        if not os.path.isdir(prefix):
            print( 'No LLVM:{} + Binder install is detected! Going to check out LLVM and install Binder. This procedure will require ~1Gb of free disk space and will only be needed to be done once...\n'.format(release) )
            os.makedirs(prefix)

        if not os.path.isdir(prefix+'/.git'): execute('Clonning llvm...', 'cd {} && git clone git@github.com:llvm-mirror/llvm.git .'.format(prefix) )
        execute('Checking out LLVM revision: {}...'.format(release), 'cd {prefix} && ( {git_checkout} || ( git fetch && {git_checkout} ) )'.format(prefix=prefix, git_checkout=git_checkout) )

        if not os.path.isdir(prefix+'/tools/clang'): execute('Clonning clang...', 'cd {}/tools && git clone git@github.com:llvm-mirror/clang.git clang'.format(prefix) )
        execute('Checking out Clang revision: {}...'.format(release), 'cd {prefix}/tools/clang && ( {git_checkout} || ( git fetch && {git_checkout} ) )'.format(prefix=prefix, git_checkout=git_checkout) )

        if not os.path.isdir(prefix+'/tools/clang/tools/extra'): os.makedirs(prefix+'/tools/clang/tools/extra')
        '''

        tool_link_path = '{prefix}/tools/clang/tools/extra/{name}'.format(prefix=prefix, name=name)
        if os.path.islink(tool_link_path): os.unlink(tool_link_path)
        os.symlink(source_location, tool_link_path)

        cmake_lists = prefix + '/tools/clang/tools/extra/CMakeLists.txt'
        tool_build_line = 'add_subdirectory({})'.format(name)

        if not os.path.isfile(cmake_lists):
            with open(cmake_lists, 'w') as f: f.write(tool_build_line + '\n')

        build_config = '-DCMAKE_BUILD_TYPE={build_type}'.format(build_type='Debug' if debug else 'Release')
        build_config += get_cmake_compiler_options()

        if not os.path.isdir(build_dir): os.makedirs(build_dir)
        execute(
            'Building tool: {}...'.format(name), # -DLLVM_TEMPORARILY_ALLOW_OLD_TOOLCHAIN=1
            'cd {build_dir} && cmake -G Ninja {build_config} -DLLVM_ENABLE_EH=1 -DLLVM_ENABLE_RTTI=ON {gcc_install_prefix} .. && ninja {name} tools/clang/lib/Headers/clang-headers {jobs}'.format( # was 'binder clang', we need to build Clang so lib/clang/<version>/include is also built
                name=name,
                build_dir=build_dir, build_config=build_config,
                jobs=f"-j{config['cpu_count']}",
                gcc_install_prefix = ''), # ='-DGCC_INSTALL_PREFIX='+Options.gcc_install_prefix if Options.gcc_install_prefix else ''),
            silence_output=True)

        print()
        # build_dir = prefix+'/build-ninja-' + release
        # if not os.path.isdir(build_dir): os.makedirs(build_dir)
        # execute('Building tool: {}...'.format(name), 'cd {build_dir} && cmake -DCMAKE_BUILD_TYPE={build_type} .. -G Ninja && ninja -j{jobs}'.format(build_dir=build_dir, jobs=Options.jobs, build_type='Debug' if debug else 'Release')) )

        with open(signature_file_name, 'w') as f: json.dump(signature, f, sort_keys=True, indent=2)

    executable = build_dir + '/bin/' + name
    if not os.path.isfile(executable): print("\nEncounter error while running install_llvm_tool: Build is complete but executable {} is not there!!!".format(executable) ); sys.exit(1)

    return executable



def get_python_include_and_lib(python):
    ''' calculate python include dir and lib dir from given python executable path
    '''
    #python = os.path.realpath(python)
    python_bin_dir = python.rpartition('/')[0]
    python_config = f'{python} {python}-config' if python.endswith('2.7') else f'{python}-config'

    #if not os.path.isfile(python_config): python_config = python_bin_dir + '/python-config'

    info = execute('Getting python configuration info...', f'unset __PYVENV_LAUNCHER__ && cd {python_bin_dir} && PATH=.:$PATH && {python_config} --prefix --includes', return_='output').replace('\r', '').split('\n')  # Python-3 only: --abiflags
    python_prefix = info[0]
    python_include_dir = info[1].split()[0][len('-I'):]
    python_lib_dir = python_prefix + '/lib'
    #python_abi_suffix = info[2]
    #print(python_include_dir, python_lib_dir)

    return NT(python_include_dir=python_include_dir, python_lib_dir=python_lib_dir)


def local_open_ssl_install(prefix, build_prefix, jobs):
    ''' install OpenSSL at given prefix, return url of source archive
    '''
    #with tempfile.TemporaryDirectory('open_ssl_build', dir=prefix) as build_prefix:

    url = 'https://www.openssl.org/source/openssl-1.1.1b.tar.gz'
    #url = 'https://www.openssl.org/source/openssl-3.0.0.tar.gz'


    archive = build_prefix + '/' + url.split('/')[-1]
    build_dir = archive.rpartition('.tar.gz')[0]
    if os.path.isdir(build_dir): shutil.rmtree(build_dir)

    with open(archive, 'wb') as f:
        response = urllib.request.urlopen(url)
        f.write( response.read() )

    execute('Unpacking {}'.format(archive), 'cd {build_prefix} && tar -xvzf {archive}'.format(**vars()) )

    execute('Configuring...', f'cd {build_dir} && ./config --prefix={prefix}')
    execute('Building...',    f'cd {build_dir} && make -j{jobs}')
    execute('Installing...',  f'cd {build_dir} && make -j{jobs} install')

    return url


def remove_pip_and_easy_install(prefix_root_path):
    ''' remove `pip` and `easy_install` executable from given Python / virtual-environments install
    '''
    for f in os.listdir(prefix_root_path + '/bin'):  # removing all pip's and easy_install's to make sure that environment is immutable
        for p in ['pip', 'easy_install']:
            if f.startswith(p): os.remove(prefix_root_path + '/bin/' + f)



def local_python_install(platform, config):
    ''' Perform local install of given Python version and return path-to-python-interpreter, python_include_dir, python_lib_dir
        If previous install is detected skip installiation.
        Provided Python install will _persistent_ and _immutable_
    '''
    jobs = config['cpu_count']
    compiler, cpp_compiler = ('clang', 'clang++') if platform['os'] == 'mac' else ('gcc', 'g++')  # disregarding platform compiler setting and instead use default compiler for platform

    python_version = platform.get('python', DEFAULT_PYTHON_VERSION)

    if python_version.endswith('.s'):
        assert python_version == f'{sys.version_info.major}.{sys.version_info.minor}.s'
        #root = executable.rpartition('/bin/python')[0]
        h = hashlib.md5(); h.update( (sys.executable + sys.version).encode('utf-8', errors='backslashreplace') ); hash = h.hexdigest()
        return NT(
            python = sys.executable,
            root = None,
            python_include_dir = None,
            python_lib_dir = None,
            version = python_version,
            url = None,
            platform = platform,
            config = config,
            hash = hash,
        )

    # deprecated, no longer needed
    # python_version = {'python2'   : '2.7',
    #                   'python2.7' : '2.7',
    #                   'python3'   : '3.5',
    # }.get(python_version, python_version)

    # for security reasons we only allow installs for version listed here with hand-coded URL's
    python_sources = {
        '2.7' : 'https://www.python.org/ftp/python/2.7.18/Python-2.7.18.tgz',

        '3.5'  : 'https://www.python.org/ftp/python/3.5.9/Python-3.5.9.tgz',
        '3.6'  : 'https://www.python.org/ftp/python/3.6.15/Python-3.6.15.tgz',
        '3.7'  : 'https://www.python.org/ftp/python/3.7.14/Python-3.7.14.tgz',
        '3.8'  : 'https://www.python.org/ftp/python/3.8.14/Python-3.8.14.tgz',
        '3.9'  : 'https://www.python.org/ftp/python/3.9.14/Python-3.9.14.tgz',
        '3.10' : 'https://www.python.org/ftp/python/3.10.10/Python-3.10.10.tgz',
        '3.11' : 'https://www.python.org/ftp/python/3.11.2/Python-3.11.2.tgz',
        '3.12' : 'https://www.python.org/ftp/python/3.12.0/Python-3.12.0.tgz',
    }

    # map of env -> ('shell-code-before ./configure', 'extra-arguments-for-configure')
    extras = {
        ('m1',) :           ('unset CPATH LDFLAGS CPPFLAGS && MACOSX_DEPLOYMENT_TARGET={}'.format(platform_module.mac_ver()[0]), ''),
        ('mac',) :          ('MACOSX_DEPLOYMENT_TARGET={}'.format(platform_module.mac_ver()[0]), ''),
        ('linux',  '2.7') : ('', '--enable-unicode=ucs4'),
        ('ubuntu', '2.7') : ('', '--enable-unicode=ucs4'),
    }

    #packages = '' if (python_version[0] == '2' or  python_version == '3.5' ) and  platform['os'] == 'mac' else 'pip setuptools wheel' # 2.7 is now deprecated on Mac so some packages could not be installed
    packages = 'setuptools'

    url = python_sources[python_version]

    extra = extras.get( (platform['os'],)  , ('', '') )
    extra = extras.get( (platform['os'], python_version) , extra)

    extra = ('unset __PYVENV_LAUNCHER__ && ' + extra[0], extra[1])

    options = '--with-ensurepip' #'--without-ensurepip'
    signature = f'v1.5.1 url: {url}\noptions: {options}\ncompiler: {compiler}\nextra: {extra}\npackages: {packages}\n'

    h = hashlib.md5(); h.update( signature.encode('utf-8', errors='backslashreplace') ); hash = h.hexdigest()

    root = calculate_unique_prefix_path(platform, config) + '/python-' + python_version + '.' +  compiler + '/' + hash

    signature_file_name = root + '/.signature'

    #activate   = root + '/bin/activate'
    executable = root + '/bin/python' + python_version

    # if os.path.isfile(executable)  and  (not execute('Getting python configuration info...', '{executable}-config --prefix --includes'.format(**vars()), terminate_on_failure=False) ):
    #     print('found executable!')
    #     _, executable_version = execute('Checking Python interpreter version...', '{executable} --version'.format(**vars()), return_='tuple')
    #     executable_version = executable_version.split()[-1]
    # else: executable_version = ''
    # print('executable_version: {}'.format(executable_version))
    #if executable_version != url.rpartition('Python-')[2][:-len('.tgz')]:

    if os.path.isfile(signature_file_name) and open(signature_file_name).read() == signature:
        #print('Install for Python-{} is detected, skipping installation procedure...'.format(python_version))
        pass

    else:
        print( 'Installing Python-{python_version}, using {url} with extra:{extra}...'.format( **vars() ) )

        if os.path.isdir(root): shutil.rmtree(root)

        build_prefix = os.path.abspath(root + '/../build-python-{}'.format(python_version) )

        if not os.path.isdir(root): os.makedirs(root)
        if not os.path.isdir(build_prefix): os.makedirs(build_prefix)

        platform_is_mac = True if platform['os'] in ['mac', 'm1'] else False
        platform_is_linux = not platform_is_mac

        #if False and platform['os'] == 'mac' and platform_module.machine() == 'arm64' and tuple( map(int, python_version.split('.') ) ) >= (3, 9):
        if ( platform['os'] == 'mac' and python_version == '3.6' ) \
           or ( platform_is_linux and python_version in ['3.10', '3.11', '3.12'] ):
            open_ssl_url = local_open_ssl_install(root, build_prefix, jobs)
            options += f' --with-openssl={root} --with-openssl-rpath=auto'
            #signature += 'OpenSSL install: ' + open_ssl_url + '\n'

        archive = build_prefix + '/' + url.split('/')[-1]
        build_dir = archive.rpartition('.tgz')[0]
        if os.path.isdir(build_dir): shutil.rmtree(build_dir)

        with open(archive, 'wb') as f:
            #response = urllib2.urlopen(url)
            response = urllib.request.urlopen(url)
            f.write( response.read() )

        #execute('Execution environment:', 'env'.format(**vars()) )

        execute('Unpacking {}'.format(archive), 'cd {build_prefix} && tar -xvzf {archive}'.format(**vars()) )

        #execute('Building and installing...', 'cd {} && CC={compiler} CXX={cpp_compiler} {extra[0]} ./configure {extra[1]} --prefix={root} && {extra[0]} make -j{jobs} && {extra[0]} make install'.format(build_dir, **locals()) )
        execute('Configuring...', 'cd {} && CC={compiler} CXX={cpp_compiler} {extra[0]} ./configure {options} {extra[1]} --prefix={root}'.format(build_dir, **locals()) )
        execute('Building...', 'cd {} && {extra[0]} make -j{jobs}'.format(build_dir, **locals()) )
        execute('Installing...', 'cd {} && {extra[0]} make -j{jobs} install'.format(build_dir, **locals()) )

        shutil.rmtree(build_prefix)

        #execute('Updating setuptools...', f'cd {root} && {root}/bin/pip{python_version} install --upgrade setuptools wheel' )

        # if 'certifi' not in packages:
        #     packages += ' certifi'

        if packages: execute( f'Installing packages {packages}...', f'cd {root} && unset __PYVENV_LAUNCHER__ && {root}/bin/pip{python_version} install --upgrade {packages}' )
        #if packages: execute( f'Installing packages {packages}...', f'cd {root} && unset __PYVENV_LAUNCHER__ && {executable} -m pip install --upgrade {packages}' )

        remove_pip_and_easy_install(root)  # removing all pip's and easy_install's to make sure that environment is immutable

        with open(signature_file_name, 'w') as f: f.write(signature)

        print( 'Installing Python-{python_version}, using {url} with extra:{extra}... Done.'.format( **vars() ) )

    il = get_python_include_and_lib(executable)

    return NT(
        python = executable,
        root = root,
        python_include_dir = il.python_include_dir,
        python_lib_dir = il.python_lib_dir,
        version = python_version,
        url = url,
        platform = platform,
        config = config,
        hash = hash,
    )



def setup_python_virtual_environment(working_dir, python_environment, packages=''):
    ''' Deploy Python virtual environment at working_dir
    '''

    python = python_environment.python

    execute('Setting up Python virtual environment...', 'unset __PYVENV_LAUNCHER__ && {python} -m venv --clear {working_dir}'.format(**vars()) )

    activate = f'unset __PYVENV_LAUNCHER__ && . {working_dir}/bin/activate'

    bin=working_dir+'/bin'

    if packages: execute('Installing packages: {}...'.format(packages), 'unset __PYVENV_LAUNCHER__ && {bin}/python {bin}/pip install --upgrade pip setuptools && {bin}/python {bin}/pip install --progress-bar off {packages}'.format(**vars()) )
    #if packages: execute('Installing packages: {}...'.format(packages), '{bin}/pip{python_environment.version} install {packages}'.format(**vars()) )

    return NT(activate = activate, python = bin + '/python', root = working_dir, bin = bin)



def setup_persistent_python_virtual_environment(python_environment, packages):
    ''' Setup _persistent_ and _immutable_ Python virtual environment which will be saved between test runs
    '''

    if python_environment.version.startswith('2.'):
        assert not packages, f'ERROR: setup_persistent_python_virtual_environment does not support Python-2.* with non-empty package list!'
        return NT(activate = ':', python = python_environment.python, root = python_environment.root, bin = python_environment.root + '/bin')

    else:
        #if 'certifi' not in packages: packages += ' certifi'

        h = hashlib.md5()
        h.update(f'v1.0.0 platform: {python_environment.platform} python_source_url: {python_environment.url} python-hash: {python_environment.hash} packages: {packages}'.encode('utf-8', errors='backslashreplace') )
        hash = h.hexdigest()

        prefix = calculate_unique_prefix_path(python_environment.platform, python_environment.config)

        root = os.path.abspath( prefix + '/python_virtual_environments/' + '/python-' + python_environment.version + '/' + hash )
        signature_file_name = root + '/.signature'
        signature = f'setup_persistent_python_virtual_environment v1.0.0\npython: {python_environment.hash}\npackages: {packages}\n'

        activate = f'unset __PYVENV_LAUNCHER__ && . {root}/bin/activate'
        bin = f'{root}/bin'

        if os.path.isfile(signature_file_name) and open(signature_file_name).read() == signature: pass
        else:
            if os.path.isdir(root): shutil.rmtree(root)
            setup_python_virtual_environment(root, python_environment, packages=packages)
            remove_pip_and_easy_install(root)  # removing all pip's and easy_install's to make sure that environment is immutable
            with open(signature_file_name, 'w') as f: f.write(signature)

        return NT(activate = activate, python = bin + '/python', root = root, bin = bin, hash = hash)



def generate_version_information(rosetta_dir, **kwargs):
    ''' Generate standard Rosetta version structure and save it to JSON file if file_name is provided.

    This is a light wrapper around the generate_version_information() function in Rosetta/main/source/version.py -- see there for the interface definition.
    '''

    version = imp.load_source('version', rosetta_dir + '/source/version.py')
    return version.generate_version_information(rosetta_dir=rosetta_dir, **kwargs)



def _get_path_to_conda_root(platform, config):
    ''' Perform local (prefix) install of miniconda and return NT(activate, conda_root_dir, conda)
        this function is for inner use only, - to setup custom conda environment inside your test use `setup_conda_virtual_environment` defined below
    '''
    miniconda_sources = {
        'mac'    : 'https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh',
        'linux'  : 'https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh',
        'aarch64': 'https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh',
        'ubuntu' : 'https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh',
        'm1'     : 'https://repo.anaconda.com/miniconda/Miniconda3-py38_4.10.1-MacOSX-arm64.sh',
    }

    conda_sources = {
        'mac'    : 'https://repo.continuum.io/archive/Anaconda3-2018.12-MacOSX-x86_64.sh',
        'linux'  : 'https://repo.continuum.io/archive/Anaconda3-2018.12-Linux-x86_64.sh',
        'ubuntu' : 'https://repo.continuum.io/archive/Anaconda3-2018.12-Linux-x86_64.sh',
    }

    #platform_os = 'm1' if platform_module.machine() == 'arm64' else platform['os']
    #url = miniconda_sources[ platform_os ]

    platform_os = platform['os']
    for o in 'alpine centos ubuntu'.split():
        if platform_os.startswith(o): platform_os = 'linux'

    url = miniconda_sources[platform_os]

    version = '1'
    channels = ''  # conda-forge

    #packages = ['conda-build gcc libgcc', 'libgcc=5.2.0'] # libgcc installs is workaround for "Anaconda libstdc++.so.6: version `GLIBCXX_3.4.20' not found", see: https://stackoverflow.com/questions/48453497/anaconda-libstdc-so-6-version-glibcxx-3-4-20-not-found
    #packages = ['conda-build gcc'] # libgcc installs is workaround for "Anaconda libstdc++.so.6: version `GLIBCXX_3.4.20' not found", see: https://stackoverflow.com/questions/48453497/anaconda-libstdc-so-6-version-glibcxx-3-4-20-not-found
    packages = ['conda-build anaconda-client conda-verify',]

    signature = f'url: {url}\nversion: {version}\channels: {channels}\npackages: {packages}\n'

    root = calculate_unique_prefix_path(platform, config) + '/conda'

    signature_file_name = root + '/.signature'

    # presense of __PYVENV_LAUNCHER__,PYTHONHOME, PYTHONPATH sometimes confuse Python so we have to unset them
    unset = 'unset __PYVENV_LAUNCHER__ && unset PYTHONHOME && unset PYTHONPATH'
    activate = unset + ' && . ' + root + '/bin/activate'

    executable = root + '/bin/conda'


    if os.path.isfile(signature_file_name) and open(signature_file_name).read() == signature:
        print( f'Install for MiniConda is detected, skipping installation procedure...' )

    else:
        print( f'Installing MiniConda, using {url}...' )

        if os.path.isdir(root): shutil.rmtree(root)

        build_prefix = os.path.abspath(root + f'/../build-conda' )

        #if not os.path.isdir(root): os.makedirs(root)
        if not os.path.isdir(build_prefix): os.makedirs(build_prefix)

        archive = build_prefix + '/' + url.split('/')[-1]

        with open(archive, 'wb') as f:
            response = urllib.request.urlopen(url)
            f.write( response.read() )

        execute('Installing conda...', f'cd {build_prefix} && {unset} && bash {archive} -b -p {root}' )

        # conda update --yes --quiet -n base -c defaults conda

        if channels: execute(f'Adding extra channles {channels}...', f'cd {build_prefix} && {activate} && conda config --add channels {channels}' )

        for p in packages: execute(f'Installing conda packages: {p}...', f'cd {build_prefix} && {activate} && conda install --quiet --yes {p}' )

        shutil.rmtree(build_prefix)

        with open(signature_file_name, 'w') as f: f.write(signature)

        print( f'Installing MiniConda, using {url}... Done.' )

    execute(f'Updating conda base...', f'{activate} && conda update --all --yes' )
    return NT(conda=executable, root=root, activate=activate, url=url)



def setup_conda_virtual_environment(working_dir, platform, config, packages=''):
    ''' Deploy Conda virtual environment at working_dir
    '''
    conda_root_env = _get_path_to_conda_root(platform, config)
    activate = conda_root_env.activate

    python_version = platform.get('python', DEFAULT_PYTHON_VERSION)

    prefix = os.path.abspath( working_dir + '/.conda-python-' + python_version )

    command_line = f'conda create --quiet --yes --prefix {prefix} python={python_version}'

    execute( f'Setting up Conda for Python-{python_version} virtual environment...', f'cd {working_dir} && {activate} && ( {command_line} || ( conda clean --yes && {command_line} ) )' )

    activate = f'{activate} && conda activate {prefix}'

    if packages: execute( f'Setting up extra packages {packages}...', f'cd {working_dir} && {activate} && conda install --quiet --yes {packages}' )

    python = prefix + '/bin/python' + python_version

    il = get_python_include_and_lib(python)

    return NT(
        activate = activate,
        root = prefix,
        python = python,
        python_include_dir = il.python_include_dir,
        python_lib_dir = il.python_lib_dir,
        version = python_version,
        activate_base = conda_root_env.activate,
        url = prefix, # conda_root_env.url,
        platform=platform,
        config=config,
    )



class FileLock():
    ''' Implementation of file-lock object that could be use with Python `with` statement
    '''

    def __init__(self, file_name):
        self.locked = False
        self.file_name = file_name


    def __enter__(self):
        if not self.locked: self.acquire()
        return self


    def __exit__(self, exc_type, exc_value, traceback):
        if self.locked: self.release()


    def __del__(self):
        self.release()


    def acquire(self):
        while True:
            try:
                os.close( os.open(self.file_name, os.O_CREAT | os.O_EXCL, mode=0o600) )
                self.locked = True
                break

            except FileExistsError as e:
                time.sleep(60)


    def release(self):
        if self.locked:
            os.remove(self.file_name)
            self.locked = False



def convert_submodule_urls_from_ssh_to_https(repository_root):
    ''' switching submodules URL to HTTPS so we can clone without SSH key
    '''
    with open(f'{repository_root}/.gitmodules') as f: m = f.read()
    with open(f'{repository_root}/.gitmodules', 'w') as f:
        f.write(
            m
            .replace('url = git@github.com:',       'url = https://github.com/')
            .replace('url = ../../RosettaCommons/', 'url = https://github.com/RosettaCommons/')
            .replace('url = ../../../',             'url = https://github.com/RosettaCommons/')
            .replace('url = ../../',                'url = https://github.com/RosettaCommons/')
            .replace('url = ../',                   'url = https://github.com/RosettaCommons/')
        )
