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


import commands
import codecs

# ⚔ do not change wording below, it have to stay in sync with upstream (up to benchmark-model).
# Copied from benchmark-model, standard state code's for tests results.

__all__ = ['_S_Values_', '_S_draft_', '_S_queued_', '_S_running_', '_S_passed_', '_S_failed_', '_S_build_failed_', '_S_script_failed_',
           '_StateKey_', '_ResultsKey_', '_LogKey_'
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

PyRosetta_unix_memory_requirement_per_cpu = 2.5  # Memory per sub-process in Gb's
PyRosetta_unix_unit_test_memory_requirement_per_cpu = 3.0  # Memory per sub-process in Gb's for running PyRosetta unit tests


# Standard funtions and classes below ---------------------------------------------------------------------------------

class BenchmarkError(Exception):
    def __init__(self, value): self.value = value
    def __str__(self): return self.value


class NT:  # named tuple
    def __init__(self, **entries): self.__dict__.update(entries)
    def __repr__(self):
        r = 'NT: |'
        for i in dir(self):
            if not i.startswith('__') and not isinstance(getattr(self, i), types.MethodType): r += '%s --> %s, ' % (i, getattr(self, i))
        return r[:-2]+'|'


def Tracer(verbose=False):
    def print_(x): print x
    return print_ if verbose else lambda x: None


def execute(message, commandline, return_=False, until_successes=False, terminate_on_failure=True):
    TR = Tracer()
    TR(message);  TR(commandline)
    while True:
        (res, output) = commands.getstatusoutput(commandline)
        # Subprocess results will always be a bytes-string.
        # Probably ASCII, but may have some Unicode characters.
        # A UTF-8 decode will probably get decent results 99% of the time
        # and the replace option will gracefully handle the rest.
        output = output.decode(encoding="utf-8", errors="replace")

        TR(output)

        if res and until_successes: pass  # Thats right - redability COUNT!
        else: break

        print "Error while executing %s: %s\n" % (message, output)
        print "Sleeping 60s... then I will retry..."
        time.sleep(60)

    if return_ == 'tuple': return(res, output)

    if res and terminate_on_failure:
        TR("\nEncounter error while executing: " + commandline)
        if return_==True: return True
        else: raise BenchmarkError("\nEncounter error while executing: " + commandline + '\n' + output)

    if return_ == 'output': return output
    else: return False


def parallel_execute(name, jobs, rosetta_dir, working_dir, cpu_count, time=16):
    ''' Execute command line in parallel on local host
        time is specify upper time limit in minutes after which jobs will be automatically terminated

        jobs should be dict with following structure:
        {
            ‘job-string-id-1’: command_line-1,
            ‘job-string-id-2’: command_line-2,
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



def build_rosetta(rosetta_dir, platform, jobs, mode='release', verbose=False, debug=False):
    ''' Compile Rosetta binaries on a given platform return (res, output, build_command_line) '''

    compiler = platform['compiler']
    extras   = ','.join(platform['extras'])

    # removing all symlinks from bin/ and then building binaries...
    build_command_line = 'find bin -type l ! -name ".*" -exec rm {{}} \\; && ./scons.py bin mode={mode} cxx={compiler} extras={extras} -j{jobs}'.format(jobs=jobs, mode=mode, compiler=compiler, extras=extras)

    if debug:
        res, output = 0, '__init__.py:build_rosetta: debug is enabled, skipping build...\n'
    else:
        res, output = execute('Compiling...', 'cd {}/source && {}'.format(rosetta_dir, build_command_line), return_='tuple')

    return res, output, build_command_line



def build_pyrosetta(rosetta_dir, platform, jobs, config, mode='MinSizeRel', verbose=False, debug=False):
    ''' Compile Rosetta binaries on a given platform return (res, output, build_command_line, pyrosetta_path) '''

    #binder = install_llvm_tool('binder', source_location='{}/source/src/python/PyRosetta/binder'.format(rosetta_dir), config=config)

    command_line = 'cd {rosetta_dir}/source/src/python/PyRosetta && {python} build.py -j{jobs} --compiler {compiler} --type {mode}'.format(rosetta_dir=rosetta_dir, python=platform['python'], jobs=jobs, compiler=platform['compiler'], mode=mode)

    pyrosetta_path = execute('Getting PyRosetta build path...', command_line + ' --print-build-root', return_='output')

    if debug:
        res, output = 0, '__init__.py:build_pyrosetta: debug is enabled, skipping build...\n'
    else:
        res, output = execute('Building PyRosetta {}...'.format(mode), command_line, return_='tuple')

    return res, output, command_line, pyrosetta_path



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
    tool_build_line = 'add_subdirectory({})'.format(name)

    for line in codecs.open(cmake_lists, encoding='utf-8', errors='replace'):
        if line == tool_build_line: break
    else:
        with codecs.open(cmake_lists, 'w', encoding='utf-8', errors='replace') as f: f.write( codecs.open(cmake_lists, encoding='utf-8', errors='replace').read() + tool_build_line + '\n' )

    build_dir = prefix+'/llvm/build-' + release
    if not os.path.isdir(build_dir): os.makedirs(build_dir)
    execute('Building tool: {}...'.format(name), 'cd {build_dir} && cmake -DCMAKE_BUILD_TYPE=Release .. && make -j{jobs}'.format(build_dir=build_dir, jobs=jobs) )

    # build_dir = prefix+'/llvm/build-ninja-' + release
    # if not os.path.isdir(build_dir): os.makedirs(build_dir)
    # execute('Building tool: {}...'.format(name), 'cd {build_dir} && cmake -DCMAKE_BUILD_TYPE=Release .. -G Ninja && ninja -j{jobs}'.format(build_dir=build_dir, jobs=jobs) )

    executable = build_dir + '/bin/' + name
    if not os.path.isfile(executable): raise BenchmarkError("\nEncounter error while running install_llvm_tool: Build is complete but executable {} is not there!!!".format(executable) )

    return executable
