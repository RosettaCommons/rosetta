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

import os, time, sys, shutil, codecs, urllib.request, imp, subprocess  # urllib.error, urllib.parse,
import platform as  platform_module
import types as types_module

# ⚔ do not change wording below, it have to stay in sync with upstream (up to benchmark-model).
# Copied from benchmark-model, standard state code's for tests results.

__all__ = ['execute',
           '_S_Values_', '_S_draft_', '_S_queued_', '_S_running_', '_S_passed_', '_S_failed_', '_S_build_failed_', '_S_script_failed_',
           '_StateKey_', '_ResultsKey_', '_LogKey_',
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

PyRosetta_unix_memory_requirement_per_cpu = 2.5  # Memory per sub-process in Gb's
PyRosetta_unix_unit_test_memory_requirement_per_cpu = 3.0  # Memory per sub-process in Gb's for running PyRosetta unit tests


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
        time is specify upper time limit in minutes after which jobs will be automatically terminated

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
    allowed_time = int(time*60)
    job_file_name = working_dir + '/' + name
    with file(job_file_name + '.json', 'w') as f: json.dump(jobs, f, sort_keys=True, indent=2) # JSON handles unicode internally
    execute("Running {} in parallel with {} CPU's...".format(name, cpu_count), 'cd {working_dir} && ulimit -t {allowed_time} && {rosetta_dir}/tests/benchmark/util/parallel.py -j{cpu_count} {job_file_name}.json'.format(**vars()))

    return json.load( file(job_file_name+'.results.json') )


def calculate_extension(platform, mode='release'):
    ''' Calculate and return extension for Rosetta executables that will match one generated by scons '''

    if platform['os'] in 'linux ubuntu': os = 'linux'
    elif platform['os'] == 'mac':        os = 'macos'
    else:                                os = platform['os']

    extras = ''.join(platform['extras']) if platform['extras'] else 'default'  # platform and therfor extras should be comming from database and therefor it should be already sorted

    return extras + "." + os + platform['compiler'] + mode


def platform_to_pretty_string(platform):
    ''' Take platform as json object and return normalized human-readable string '''
    return '{}.{}{}{}'.format(platform['os'], platform['compiler'], ('.'+'.'.join(platform['extras']) if 'extras' in platform  and  platform['extras'] else ''), ('.python'+platform['python'] if 'python' in platform else ''))


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

    pyrosetta_path = execute('Getting PyRosetta build path...', command_line + ' --print-build-root', return_='output').split()[0]

    if skip_compile:
        res, output = 0, '__init__.py:build_pyrosetta: skip_compile is enabled, skipping build...\n'
    else:
        res, output = execute('Building PyRosetta {}...'.format(mode), command_line, return_='tuple', add_message_and_command_line_to_output=True)

    return NT(exitcode=res, output=output, command_line=command_line, pyrosetta_path=pyrosetta_path, python=py_env.python, root=py_env.root)



def install_llvm_tool(name, source_location, config, clean=True):
    ''' Install and update (if needed) custom LLVM tool at given prefix (from config).
        Return absolute path to executable on success and raise BenchmarkError exception on failure (do not catch this! if you really need 'normal' exit from this function on failure - refactor it instead)
    '''
    prefix = config['prefix']
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


def get_python_include_and_lib(python):
    ''' calculate python include dir and lib dir from given python executable executable
    '''
    python_bin_dir = python.rpartition('/')[0]
    info = execute('Getting python configuration info...', 'unset __PYVENV_LAUNCHER__ && cd {python_bin_dir} && PATH=.:$PATH && {python}-config --prefix --includes'.format(**vars()), return_='output').replace('\r', '').split('\n')  # Python-3 only: --abiflags
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



def local_python_install(platform, config):
    ''' Perform local install of given Python version and return path-to-python-interpreter, python_include_dir, python_lib_dir
        If previous install is detected skip installiation.
    '''
    prefix = config['prefix']
    jobs = config['cpu_count']
    compiler, cpp_compiler = ('clang', 'clang++') if platform['os'] == 'mac' else ('gcc', 'g++')  # disregarding platform compiler setting and instead use default compiler for platform

    python_version = platform.get('python', '3.6')

    # deprecated, no longer needed
    # python_version = {'python2'   : '2.7',
    #                   'python2.7' : '2.7',
    #                   'python3'   : '3.5',
    # }.get(python_version, python_version)

    # for security reasons we only allow installs for version listed here with hand-coded URL's
    python_sources = {
        '2.7' : 'https://www.python.org/ftp/python/2.7.14/Python-2.7.14.tgz',

        '3.5' : 'https://www.python.org/ftp/python/3.5.5/Python-3.5.5.tgz',
        '3.6' : 'https://www.python.org/ftp/python/3.6.8/Python-3.6.8.tgz',
        '3.7' : 'https://www.python.org/ftp/python/3.7.3/Python-3.7.3.tgz',
    }

    # map of env -> ('shell-code-before ./configure', 'extra-arguments-for-configure')
    extras = {
        #('mac',) :          ('__PYVENV_LAUNCHER__="" MACOSX_DEPLOYMENT_TARGET={}'.format(platform_module.mac_ver()[0]), ''),  # __PYVENV_LAUNCHER__ now used by-default for all platform installs
        ('mac',) :          ('MACOSX_DEPLOYMENT_TARGET={}'.format(platform_module.mac_ver()[0]), ''),
        ('linux',  '2.7') : ('', '--enable-unicode=ucs4'),
        ('ubuntu', '2.7') : ('', '--enable-unicode=ucs4'),
    }

    packages = None #if python_version.startswith('2.') else 'pydoc'  # Python-2 does not install pip as default

    url = python_sources[python_version]

    extra = extras.get( (platform['os'],)  , ('', '') )
    extra = extras.get( (platform['os'], python_version) , extra)

    extra = ('unset __PYVENV_LAUNCHER__ && ' + extra[0], extra[1])

    signature = 'v1.1 url: {url}\ncompiler: {compiler}\nextra: {extra}\npackages: {packages}\n'.format( **vars() )

    machine_name = os.uname()[1]
    suffix = platform['os'] + '.' + machine_name

    root = os.path.abspath(prefix + '/' + suffix + '/python-' + python_version + '.' +  compiler)

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

        if platform['os'] == 'mac' and python_version == '3.6':
            open_ssl_url = local_open_ssl_install(root, build_prefix, jobs)
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
        execute('Configuring...', 'cd {} && CC={compiler} CXX={cpp_compiler} {extra[0]} ./configure {extra[1]} --prefix={root}'.format(build_dir, **locals()) )
        execute('Building...', 'cd {} && {extra[0]} make -j{jobs}'.format(build_dir, **locals()) )
        execute('Installing...', 'cd {} && {extra[0]} make -j{jobs} install'.format(build_dir, **locals()) )

        shutil.rmtree(build_prefix)

        #if packages: execute('Installing packages {}...'.format(packages), 'cd {root} && {root}/bin/pip{python_version} install {packages}'.format(**vars()) )

        with open(signature_file_name, 'w') as f: f.write(signature)

        print( 'Installing Python-{python_version}, using {url} with extra:{extra}... Done.'.format( **vars() ) )

    il = get_python_include_and_lib(executable)

    return NT(python=executable, root=root, python_include_dir=il.python_include_dir, python_lib_dir=il.python_lib_dir, version=python_version)



def setup_python_virtual_environment(working_dir, python_environment, packages=''):
    ''' Deploy Python virtual environment at working_dir
    '''

    python = python_environment.python

    execute('Setting up Python virtual environment...', '{python} -m venv {working_dir}'.format(**vars()) )

    activate = working_dir + '/bin/activate'

    bin=working_dir+'/bin'

    if packages: execute('Installing packages: {}...'.format(packages), '{bin}/python {bin}/pip install --upgrade pip && {bin}/python {bin}/pip install {packages}'.format(**vars()) )
    #if packages: execute('Installing packages: {}...'.format(packages), '{bin}/pip{python_environment.version} install {packages}'.format(**vars()) )

    return NT(activate = activate, python = bin + '/python', root = working_dir, bin = bin)


def generate_version_information(rosetta_dir, url=None, branch=None, package=None, revision=None, date=None, file_name=None):
    ''' Generate standard Rosetta version structure and save it to JSON file if file_name is provided
    '''

    version = imp.load_source('version', rosetta_dir + '/source/version.py')
    return version.generate_version_information(rosetta_dir=rosetta_dir, url=url, branch=branch, package=package, revision=revision, date=date, file_name=file_name)



def _get_path_to_conda_root(platform, config):
    ''' Perform local (prefix) install of miniconda and return NT(activate, conda_root_dir, conda)
        this function is for inner use only, - to setup custom conda environment inside your test use `setup_conda_virtual_environment` below
    '''
    prefix = config['prefix']

    miniconda_sources = {
        'mac'    : 'https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh',
        'linux'  : 'https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh',
        'ubuntu' : 'https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh',
    }

    conda_sources = {
        'mac'    : 'https://repo.continuum.io/archive/Anaconda3-2018.12-MacOSX-x86_64.sh',
        'linux'  : 'https://repo.continuum.io/archive/Anaconda3-2018.12-Linux-x86_64.sh',
        'ubuntu' : 'https://repo.continuum.io/archive/Anaconda3-2018.12-Linux-x86_64.sh',
    }


    #url = conda_sources[ platform['os'] ]
    url = miniconda_sources[ platform['os'] ]

    version = '1'
    channels = ''  #'conda-forge'

    #packages = ['conda-build gcc libgcc', 'libgcc=5.2.0'] # libgcc installs is workaround for "Anaconda libstdc++.so.6: version `GLIBCXX_3.4.20' not found", see: https://stackoverflow.com/questions/48453497/anaconda-libstdc-so-6-version-glibcxx-3-4-20-not-found
    #packages = ['conda-build gcc'] # libgcc installs is workaround for "Anaconda libstdc++.so.6: version `GLIBCXX_3.4.20' not found", see: https://stackoverflow.com/questions/48453497/anaconda-libstdc-so-6-version-glibcxx-3-4-20-not-found
    packages = ['conda-build anaconda-client',]

    signature = f'url: {url}\nversion: {version}\channels: {channels}\npackages: {packages}\n'

    machine_name = os.uname()[1]
    suffix = platform['os'] + '.' + machine_name

    root = os.path.abspath( prefix + '/' + suffix + '/conda' )

    signature_file_name = root + '/.signature'

    # presense of __PYVENV_LAUNCHER__,PYTHONHOME, PYTHONPATH sometimes confuse Python so we have to unset them
    unset = 'unset __PYVENV_LAUNCHER__ && unset PYTHONHOME && unset PYTHONPATH'
    activate = unset + ' && source ' + root + '/bin/activate'

    executable = root + '/bin/conda'


    if os.path.isfile(signature_file_name) and open(signature_file_name).read() == signature:
        print( f'Install for MiniConda is detected, skipping installation procedure...' )
        pass

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

    return NT(conda=executable, root=root, activate=activate)



def setup_conda_virtual_environment(working_dir, platform, config, packages=''):
    ''' Deploy Conda virtual environment at working_dir
    '''
    conda_root_env = _get_path_to_conda_root(platform, config)
    activate = conda_root_env.activate

    python_version = platform.get('python', '3.6')

    prefix = os.path.abspath( working_dir + '/.conda-python-' + python_version )

    command_line = f'conda create --quiet --yes --prefix {prefix} python={python_version}'

    execute( f'Setting up Conda for Python-{python_version} virtual environment...', f'cd {working_dir} && {activate} && ( {command_line} || ( conda -y clean && {command_line} ) )' )

    activate = f'{activate} && conda activate {prefix}'

    if packages: execute( f'Setting up extra packages {packages}...', f'cd {working_dir} && {activate} && conda install --quiet --yes {packages}' )

    python = prefix + '/bin/python' + python_version

    il = get_python_include_and_lib(python)

    return NT( activate = activate, root = prefix, python = python, python_include_dir=il.python_include_dir, python_lib_dir=il.python_lib_dir, version=python_version, activate_base = conda_root_env.activate)
