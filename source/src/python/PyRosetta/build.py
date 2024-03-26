#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   build.py
## @brief  Build PyRosetta
## @author Sergey Lyskov

# https://cmake.org/files/v3.5/cmake-3.5.2.tar.gz
# /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain

from __future__ import print_function

import os, sys, argparse, platform, subprocess, shutil, codecs, json, hashlib

try:
    from setuptools.distutils import dir_util as dir_util_module
except ModuleNotFoundError:
    from distutils import dir_util as dir_util_module

from collections import OrderedDict

if sys.platform.startswith("linux"): Platform = "linux" # can be linux1, linux2, etc
elif sys.platform == "darwin" : Platform = "macos"
elif sys.platform == "cygwin" : Platform = "cygwin"
elif sys.platform == "win32" : Platform = "windows"
else: Platform = "unknown"
PlatformBits = platform.architecture()[0][:2]

_script_version_ = '1.1'  # used as seed for dependency tracking

if os.path.isfile('.hostname'):
    with open('.hostname') as f: _machine_name_ = f.read().split()[0]
else: _machine_name_ = os.getenv('HOSTNAME') or platform.node()

_python_version_ = '{}.{}'.format(sys.version_info.major, sys.version_info.minor)  # should be formatted: 2.7 or 3.5
#_python_version_ = '{}.{}.{}'.format(sys.version_info.major, sys.version_info.minor, sys.version_info.micro)  # should be formatted: 2.7.6 or 3.5.0

#_pybind11_version_ = 'fa6a4241326a361fc23915f6a82c1e48667de668'

_banned_dirs_ = 'external/bcl src/utility/pointer'.split()  # src/utility/keys src/utility/options src/basic/options src/protocols/jd3
_banned_headers_ = 'utility/py/PyHelper.hh utility/keys/KeyCount.hh utility/keys/KeyLookup.functors.hh'
_banned_headers_ +=' core/scoring/fiber_diffraction/FiberDiffractionKernelGpu.hh' # GPU code
_banned_headers_ +=' basic/database/DatabaseSessionLoader.hh' # TEMP deprecated code from the old ResourceManager ...
_banned_headers_ +=' basic/database/DatabaseSessionLoaderCreator.hh' # TEMP ... that could still be useful if slightly edited...
_banned_headers_ +=' basic/database/DatabaseSessionOptions.hh' # TEMP ... and turned into a protocols::parser::DataLoader

_banned_headers_ +=' protocols/jd3/JobOutputWritter.hh protocols/jd3/standard/PDBPoseOutputSpecification.hh protocols/jd3/standard/StandardResultOutputter.hh'  # protocols/jd3/job_distributors/MPIWorkPartitionJobDistributor.hh

#_banned_headers_ +=' json.hpp'

# Setting output to be printed in unicode regardless of terminal settings
if sys.version_info[0] == 2: sys.stdout = codecs.getwriter('utf8')(sys.stdout)
else: sys.stdout = codecs.getwriter('utf-8')(sys.stdout.detach());

# moved to config file _banned_namespaces_ = 'utility::options'.split()

def is_dir_banned(directory):
    for b in _banned_dirs_:
        if b in directory: return True
    return False


def get_rosetta_system_include_directories():
    ''' return list of include directories for compilation '''
    r = 'external external/include external/boost_submod external/dbio external/dbio/sqlite3 external/libxml2/include external/rdkit external/bcl/include'.split()
    return r

def get_rosetta_include_directories():
    ''' return list of include directories for compilation '''
    r = 'src'.split()
    r.append('src/platform/'+Platform)
    return r


def get_defines():
    ''' return list of #defines '''
    defines = 'PYROSETTA NOCRASHREPORT BOOST_ERROR_CODE_HEADER_ONLY BOOST_SYSTEM_NO_DEPRECATED BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS PTR_STD' # USEBCL
    if Platform == 'macos': defines += ' UNUSUAL_ALLOCATOR_DECLARATION'
    if Options.type in 'Release MinSizeRel': defines += ' NDEBUG BCL_NO_OS_SIGNAL_HANDLING'
    if Options.serialization: defines += ' SERIALIZATION'
    if Options.multi_threaded: defines += ' MULTI_THREADED'
    if Options.torch: defines += ' USE_PYTORCH'
    if Options.tensorflow: defines += ' USE_TENSORFLOW USE_TENSORFLOW_CPU'
    #if Options.hdf5: defines += ' USEHDF5'
    return defines.split()


def execute(message, command_line, return_='status', until_successes=False, terminate_on_failure=True, silent=False, silence_output=False):
    if not silent: print(message);  print(command_line); sys.stdout.flush();
    while True:

        p = subprocess.Popen(command_line, bufsize=0, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, errors = p.communicate()

        output = output + errors

        output = output.decode(encoding="utf-8", errors="replace")

        exit_code = p.returncode

        if exit_code  or  not (silent or silence_output): print(output); sys.stdout.flush();

        if exit_code and until_successes: pass  # Thats right - redability COUNT!
        else: break

        print( "Error while executing {}: {}\n".format(message, output) )
        print("Sleeping 60s... then I will retry...")
        sys.stdout.flush();
        time.sleep(60)

    if return_ == 'tuple': return(exit_code, output)

    if exit_code and terminate_on_failure:
        print("\nEncounter error while executing: " + command_line)
        if return_==True: return True
        else: print("\nEncounter error while executing: " + command_line + '\n' + output); sys.exit(1)

    if return_ == 'output': return output
    else: return False


def update_source_file(file_name, data):
    ''' write data to a file but only if file does not exist or it content different from data
        return True if file was written and False otherwise
    '''
    if( not os.path.isfile(file_name)  or  open(file_name).read() != data ):
        print('Writing '+ file_name)
        with open(file_name, 'w') as f: f.write(data)
        return True
    return False


def get_compiler_family():
    ''' Try to guess compiler family using Options.compiler '''
    if 'clang' in Options.compiler: return 'clang'
    if 'gcc'   in Options.compiler: return 'gcc'
    if 'cl'    in Options.compiler: return 'cl'
    return 'unknown'

def get_cmake_compiler_options():
    ''' Get cmake compiler flags from Options.compiler '''
    if Platform == "linux" and Options.compiler == 'clang': return ' -DCMAKE_C_COMPILER=`which clang` -DCMAKE_CXX_COMPILER=`which clang++`'
    if Platform == "linux" and Options.compiler == 'gcc': return ' -DCMAKE_C_COMPILER=`which gcc` -DCMAKE_CXX_COMPILER=`which g++`'

    return ''

def install_llvm_tool(name, source_location, prefix_root, debug, compiler, jobs, gcc_install_prefix, clean=True):
    ''' Install and update (if needed) custom LLVM tool at given prefix (from config).
        Return absolute path to executable on success and terminate with error on failure
    '''
    if not os.path.isdir(prefix_root): os.makedirs(prefix_root)

    # llvm_version='9.0.0'  # v8 and v9 can not be build with Clang-3.4, we if need upgrade to v > 7 then we should probably dynamicly change LLVM version based on complier versions
    # llvm_version='7.1.0'  # compiling v7.* on clang-3.4 lead to lockup while compiling tools/clang/lib/Sema/SemaChecking.cpp
    llvm_version, headers = ('13.0.0', 'tools/clang/lib/Headers/clang-resource-headers clang') if Platform == 'macos' and platform.machine() == 'arm64' else ('6.0.1', 'tools/clang/lib/Headers/clang-headers')
    #llvm_version, headers = ('13.0.0', 'tools/clang/lib/Headers/clang-resource-headers clang') if Platform == 'macos' else ('6.0.1', 'tools/clang/lib/Headers/clang-headers')
    #llvm_version, headers = ('13.0.0', 'tools/clang/lib/Headers/clang-resource-headers clang')

    if Platform == 'macos'  and  tuple( platform.mac_ver()[0].split('.') ) >= ('13','3') : llvm_version, headers = ('13.0.0', 'tools/clang/lib/Headers/clang-resource-headers clang')
    #print(f'{Platform=} {llvm_version=}')

    prefix = prefix_root + '/llvm-' + llvm_version

    build_dir = prefix+'/llvm-' + llvm_version + '.' + platform.platform() + ('.debug' if debug else '.release') # + '.' + _machine_name_

    res, output = execute('Getting binder HEAD commit SHA1...', 'cd {} && git rev-parse HEAD'.format(source_location), return_='tuple', silent=True)
    if res: binder_head = 'unknown'
    else: binder_head = output.split('\n')[0]

    signature = dict(config = 'LLVM install by install_llvm_tool version: 1.5.1, HTTPS', binder = binder_head, llvm_version=llvm_version, compiler=compiler, gcc_install_prefix=gcc_install_prefix)
    signature_file_name = build_dir + '/.signature.json'

    disk_signature = dict(config = 'unknown', binder = 'unknown')
    if os.path.isfile(signature_file_name):
        with open(signature_file_name) as f: disk_signature = json.load(f)

    if signature == disk_signature:
        print('LLVM:{} + Binder install is detected at {}, skipping LLVM installation and Binder building procedures...\n'.format(llvm_version, build_dir))

    else:
        print('LLVM build detected, but config/binder version has changed, perfoming a clean rebuild...')
        if os.path.isdir(build_dir): shutil.rmtree(build_dir)

        clang_path = "{prefix}/tools/clang".format(**locals())

        llvm_url, clang_url = {
            '6.0.1'  : ('https://releases.llvm.org/6.0.1/llvm-6.0.1.src.tar.xz', 'https://releases.llvm.org/6.0.1/cfe-6.0.1.src.tar.xz'),
            '13.0.0' : ('https://github.com/llvm/llvm-project/releases/download/llvmorg-13.0.0/llvm-13.0.0.src.tar.xz', 'https://github.com/llvm/llvm-project/releases/download/llvmorg-13.0.0/clang-13.0.0.src.tar.xz'),
            '14.0.6' : ('https://github.com/llvm/llvm-project/releases/download/llvmorg-14.0.6/llvm-14.0.6.src.tar.xz', 'https://github.com/llvm/llvm-project/releases/download/llvmorg-14.0.6/clang-14.0.6.src.tar.xz'),
            '15.0.7' : ('https://github.com/llvm/llvm-project/releases/download/llvmorg-15.0.7/llvm-15.0.7.src.tar.xz', 'https://github.com/llvm/llvm-project/releases/download/llvmorg-15.0.7/clang-15.0.7.src.tar.xz'),
            '16.0.0' : ('https://github.com/llvm/llvm-project/releases/download/llvmorg-16.0.0/llvm-16.0.0.src.tar.xz', 'https://github.com/llvm/llvm-project/releases/download/llvmorg-16.0.0/clang-16.0.0.src.tar.xz'),
        }[llvm_version]

        if not os.path.isfile(prefix + '/CMakeLists.txt'):
            #execute('Download LLVM source...', 'cd {prefix_root} && curl https://releases.llvm.org/{llvm_version}/llvm-{llvm_version}.src.tar.xz | tar -Jxom && mv llvm-{llvm_version}.src {prefix}'.format(**locals()) )
            execute('Download LLVM source...', 'cd {prefix_root} && mkdir llvm-{llvm_version}.src && curl -LJ {llvm_url} | tar --strip-components=1 -Jxom -C llvm-{llvm_version}.src && mv llvm-{llvm_version}.src {prefix}'.format(**locals()) )

        if not os.path.isdir(clang_path):
            #execute('Download Clang source...', 'cd {prefix_root} && curl https://releases.llvm.org/{llvm_version}/cfe-{llvm_version}.src.tar.xz | tar -Jxom && mv cfe-{llvm_version}.src {clang_path}'.format(**locals()) )
            execute('Download Clang source...', 'cd {prefix_root} && mkdir clang-{llvm_version}.src && curl -LJ {clang_url} | tar --strip-components=1 -Jxom -C clang-{llvm_version}.src && mv clang-{llvm_version}.src {clang_path}'.format(**locals()) )

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

        config = '-DCMAKE_BUILD_TYPE={build_type}'.format(build_type='Debug' if debug else 'Release')
        config += get_cmake_compiler_options()

        if not os.path.isdir(build_dir): os.makedirs(build_dir)
        execute(
            'Building tool: {}...'.format(name), # -DLLVM_TEMPORARILY_ALLOW_OLD_TOOLCHAIN=1
            'cd {build_dir} && cmake -G Ninja {config} -DLLVM_ENABLE_EH=1 -DLLVM_ENABLE_RTTI=ON {gcc_install_prefix} .. && ninja binder {headers} {jobs}'.format( # was 'binder clang', we need to build Clang so lib/clang/<version>/include is also built
                build_dir=build_dir, config=config,
                jobs=f'-j{jobs}' if jobs else '',
                gcc_install_prefix='-DGCC_INSTALL_PREFIX='+gcc_install_prefix if gcc_install_prefix else '',
                headers=headers,
            ),
            silence_output=True)
        print()
        # build_dir = prefix+'/build-ninja-' + release
        # if not os.path.isdir(build_dir): os.makedirs(build_dir)
        # execute('Building tool: {}...'.format(name), 'cd {build_dir} && cmake -DCMAKE_BUILD_TYPE={build_type} .. -G Ninja && ninja -j{jobs}'.format(build_dir=build_dir, jobs=Options.jobs, build_type='Debug' if debug else 'Release')) )

        with open(signature_file_name, 'w') as f: json.dump(signature, f, sort_keys=True, indent=2)

    executable = build_dir + '/bin/' + name
    if not os.path.isfile(executable): print("\nEncounter error while running install_llvm_tool: Build is complete but executable {} is not there!!!".format(executable) ); sys.exit(1)

    return executable


# def install_pybind11(prefix, clean=True):
#     ''' Download and install PyBind11 library at given prefix. Install version specified by _pybind11_version_ sha1
#     '''
#     #git_checkout = '( git fetch && git checkout {0} && git reset --hard {0} && git pull )'.format(_pybind11_version_) if clean else 'git checkout {}'.format(_pybind11_version_)
#     git_checkout = '( git fetch && git reset --hard {0} )'.format(_pybind11_version_) if clean else 'git checkout {}'.format(_pybind11_version_)

#     if not os.path.isdir(prefix): os.makedirs(prefix)
#     package_dir = prefix + '/pybind11'

#     if not os.path.isdir(package_dir): execute('Clonning pybind11...', 'cd {} && git clone git@github.com:RosettaCommons/pybind11.git'.format(prefix) )
#     execute('Checking out PyBind11 revision: {}...'.format(_pybind11_version_), 'cd {package_dir} && ( {git_checkout} )'.format(package_dir=package_dir, git_checkout=git_checkout), silent=True)
#     print()

#     include = package_dir + '/include/pybind11/pybind11.h'
#     if not os.path.isfile(include): print("\nEncounter error while running install_pybind11: Install is complete but include file {} is not there!!!".format(include) ); sys.exit(1)

#     return package_dir + '/include'


def get_compiler_version():
    compiler = get_compiler_family()

    if compiler == 'gcc':
        _, version = execute('Getting gcc version...'.format(**locals()), 'gcc -dumpfullversion -dumpversion', return_='tuple', silent=True)

    elif compiler == 'clang':
        _, version_output = execute('Getting clang version...'.format(**locals()), 'clang --version', return_='tuple', silent=True)

        if version_output.startswith('Apple clang version') or version_output.startswith('Apple LLVM version'): version = version_output.split()[3]
        elif version_output.startswith('Ubuntu clang version'): version = version_output.split()[3].partition('-')[0]
        else: version = version_output.split()[2].split('-')[0]

    else:
        version = 'unknown'

    return version.split()[0]  # removing new lines if any


def get_binding_build_root(rosetta_source_path, source=False, build=False, documentation=False):
    ''' Calculate bindings build path using current platform and compilation settings and create dir if it is not there yet '''

    p = os.path.join(rosetta_source_path, 'build/PyRosetta')

    #p =  os.path.join(p, Platform + ( '.' + Options.build_suffix if  Options.build_suffix else '') + '/' + get_compiler_family() + '-' + get_compiler_version() + '/python-' + _python_version_)
    p =  os.path.join(p, platform.platform() + '/' + get_compiler_family() + '-' + get_compiler_version() + '/python-' + _python_version_)

    p = os.path.join(p, Options.type.lower()
                     + ('.serialization' if Options.serialization else '')
                     + ('.thread' if Options.multi_threaded else '')
                     + ('.torch' if Options.torch else '')
                     + ('.tensorflow' if Options.tensorflow else '')
                     + ('.annotate' if Options.annotate_includes else '')
                     + ('.trace' if Options.trace else '')
                     )

    source_p = p + '/source'
    build_p  = p + '/build'

    if not os.path.isdir(source_p): os.makedirs(source_p)
    if not os.path.isdir(build_p): os.makedirs(build_p)

    if source: return source_p
    if build: return build_p

    if documentation: # optional, avoid creation if not requested. Also: no incremental builds for documentation so we always return clean dir
        documentation_p  = p + '/documentation'
        if os.path.isdir(documentation_p): shutil.rmtree(documentation_p)
        os.makedirs(documentation_p)
        return documentation_p

    return p


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

        for d in dirs:
            dst = dest + prefix + d
            if not os.path.isdir(dst): os.makedirs(dst)

        for f in files: symlink(source + prefix + f, dest + prefix + f)



def link_supplemental_files(rosetta_source_path):
    prefix = get_binding_build_root(rosetta_source_path, build=True)
    source = rosetta_source_path + '/src/python/PyRosetta/src'

    # remove broken symlinks if any (possibly from previous build attempts)
    for dir_name, dirs, files in os.walk(prefix):
        for p in dirs+files:
            path = dir_name + '/' + p
            if os.path.islink(path) and not os.path.exists(path): os.unlink(path)

    #distutils.dir_util.copy_tree(source, prefix, update=False)
    symlink_tree(source, prefix)

    # database_dest = prefix + '/pyrosetta/database'
    # if os.path.islink(database_dest): os.unlink(database_dest)
    # elif os.path.isdir(database_dest): shutil.rmtree(database_dest)
    # os.symlink('../../../../../../../../../database', database_dest)
    symlink(rosetta_source_path + '/../database', prefix + '/pyrosetta/database')

    symlink(rosetta_source_path + '/../LICENSE.PyRosetta.md', prefix + '/pyrosetta/LICENSE.PyRosetta.md')

    #if not os.path.islink(prefix + '/apps'): os.symlink('../../../../../../../scripts/PyRosetta/public', prefix + '/apps')  # creating link to PyRosetta apps dir
    symlink(rosetta_source_path + '/scripts/PyRosetta/public', prefix + '/apps')



def setup_source_directory_links(rosetta_source_path):
    prefix = get_binding_build_root(rosetta_source_path, source=True)

    for d in ['src', 'external']:
        source_path = os.path.relpath(os.path.join(rosetta_source_path, d), prefix)
        s = os.path.join(prefix, d)
        if os.path.islink(s): os.unlink(s)
        os.symlink(source_path, s)

    if False and Options.external_link:
        target_lib_path = os.path.relpath( os.path.join(rosetta_source_path, "cmake", "build_PyRosetta.{}".format(Options.external_link) ), prefix)

        s = os.path.join(prefix, "lib")
        if os.path.islink(s): os.unlink(s)
        os.symlink(target_lib_path, s)

    # if Platform not in ('windows', 'cygwin'):
    #     database_path = os.path.relpath(os.path.join(rosetta_source_path, "../database"), prefix)
    #     s = os.path.join(prefix, "database")
    #     if os.path.islink(s): os.unlink(s)
    #     os.symlink(database_path, s)


def generate_rosetta_external_cmake_files(rosetta_source_path, prefix):
    libs = OrderedDict([ ('cppdb', ['dbio/cppdb/atomic_counter', 'dbio/cppdb/conn_manager', 'dbio/cppdb/driver_manager', 'dbio/cppdb/frontend', 'dbio/cppdb/backend',
                                    'dbio/cppdb/mutex', 'dbio/cppdb/pool', 'dbio/cppdb/shared_object', 'dbio/cppdb/sqlite3_backend', 'dbio/cppdb/utils'] ),
                         ('sqlite3', ['dbio/sqlite3/sqlite3'] ), ])

    defines = dict(cppdb = 'CPPDB_EXPORTS CPPDB_DISABLE_SHARED_OBJECT_LOADING CPPDB_DISABLE_THREAD_SAFETY CPPDB_WITH_SQLITE3' +
                           ' CPPDB_LIBRARY_PREFIX=\\"lib\\" CPPDB_LIBRARY_SUFFIX=\\".dylib\\" CPPDB_SOVERSION=\\"0\\" ' +
                           ' CPPDB_MAJOR=0 CPPDB_MINOR=3 CPPDB_PATCH=0 CPPDB_VERSION=\\"0.3.0\\"',
                   sqlite3 = 'SQLITE_DISABLE_LFS SQLITE_OMIT_LOAD_EXTENSION SQLITE_THREADSAFE=0')

    scons_file_extension = '.external.settings'
    external_scons_files = [f for f in os.listdir(rosetta_source_path+'/external') if f.endswith(scons_file_extension)  and (Options.zmq  or  not f.startswith('libzmq.') ) and (not f.startswith('zlib.') ) ]
    for scons_file in external_scons_files:
        G = {}
        sources = []
        with open(rosetta_source_path+'/external/'+scons_file) as f:
            code = compile(f.read(), rosetta_source_path+'/external/'+scons_file, 'exec')
            exec(code, G)

        for dir_ in G['sources']:
            for f in G['sources'][dir_]:
                if '.' in f: sources.append( dir_ + '/' + f )
                else: sources.append( dir_ + '/' + f + '.cc')

            for h in os.listdir(rosetta_source_path+'/external/' + dir_):
                if h.endswith('.hh') or h.endswith('.h'): sources.append(dir_ + '/' + h)

            lib_name = scons_file[:-len(scons_file_extension)]
            libs[lib_name] = sources
            defines[lib_name] = ' '.join( G.get('defines',[]) )

    modified = False

    source_extensions = dict(sqlite3='.c', cppdb='.cpp')

    for l in libs:
        t  = 'add_library({} OBJECT\n{})\n\n'.format(l, '\n'.join( [ rosetta_source_path + '/external/' + s +  source_extensions.get(l, '') for s in libs[l]] ))
        t += 'set_property(TARGET {} PROPERTY POSITION_INDEPENDENT_CODE ON)\n'.format(l)
        if defines[l]: t += 'target_compile_options({} PRIVATE {})\n'.format(l, ' '.join( ['-D'+d for d in defines[l].split()] ) )   #  target_compile_definitions
        #t += 'target_compile_options({} PUBLIC -fPIC {})\n'.format(l, ' '.join([ '-D'+d for d in defines[l].split() ] ) )   #  target_compile_definitions

        #t += 'set_target_properties(cifparse PROPERTIES COMPILE_FLAGS "-Wno-implicit-function-declaration")\n'
        t += 'target_compile_options({l} PRIVATE "-Wno-implicit-function-declaration")\n'.format(l=l)

        modified |= update_source_file(prefix + l + '.cmake', t)

    return list( libs.keys() ), modified


def generate_rosetta_cmake_files(rosetta_source_path, prefix):
    scons_file_prefix = './../../../src/'
    lib_suffix = '.src.settings'
    libs = []

    all_libs = [ f[:-len(lib_suffix)] for f in os.listdir(scons_file_prefix)
                 if f.endswith(lib_suffix) and f not in ['apps.src.settings', 'pilot_apps.src.settings', 'devel.src.settings',]
    ]

    def key(k):
        if k.startswith('ObjexxFCL'): i = '0'
        if k.startswith('utility'):   i = '1'
        if k.startswith('numeric'):   i = '2'
        if k.startswith('basic'):     i = '3'
        if k.startswith('core'):      i = '4'
        if k.startswith('protocols'): i = '5' + k.split('.')[1]
        return i+k

    all_libs.sort(key=key, reverse=True)

    modified = False
    for lib in all_libs:
        #if lib not in 'ObjexxFCL utility numeric basic': continue
        #if lib.startswith('protocols'): continue

        G = {}
        sources = []

        with open(scons_file_prefix + lib + lib_suffix) as f:
               code = compile(f.read(), scons_file_prefix + lib + lib_suffix, 'exec')
               exec(code, G)

        for dir_ in G['sources']:
            for f in G['sources'][dir_]:
                #if not f.endswith('.cu'): sources.append( dir_ + '/' + f + '.cc')
                if f.endswith('.cu'): sources.append( dir_ + '/' + f )
                else: sources.append( dir_ + '/' + f + '.cc')

            for h in os.listdir(scons_file_prefix + dir_):
                if h.endswith('.hh') or h.endswith('.h'): sources.append(dir_ + '/' + h)


        sources.sort()



        t  = 'add_library({} OBJECT\n{})\n\n'.format(lib, '\n'.join( [ rosetta_source_path + '/src/' + s for s in sources] ))
        #t += 'set_property(TARGET {} PROPERTY POSITION_INDEPENDENT_CODE ON)\n'.format(lib)
        t += 'set_target_properties({} PROPERTIES POSITION_INDEPENDENT_CODE ON LINKER_LANGUAGE CXX)\n'.format(lib)
        #t += '\ntarget_compile_options({} PUBLIC -fPIC)\n'.format(lib)  # Enable Position Independent Code generation for libraries

        if Options.torch or Options.tensorflow: t += f'target_compile_options({lib} PRIVATE -std=c++14)\n'

        modified |= update_source_file(prefix + lib + '.cmake', t)

        libs.append(lib)

    return libs, modified



def generate_cmake_file(rosetta_source_path, extra_sources):
    ''' generate cmake file in bindings source location '''

    prefix = get_binding_build_root(rosetta_source_path, source=True) + '/'

    l1, m1 = generate_rosetta_cmake_files(rosetta_source_path, prefix)
    l2, m2 = generate_rosetta_external_cmake_files(rosetta_source_path, prefix)
    libs, modified = l1 + l2, m1 or m2

    with open('rosetta.cmake') as f: cmake = f.read()

    build_config = '# PyRosetta build config: {})'.format(
        dict(compiler = execute('Getting compiler path...', 'which gcc && which g++' if Options.compiler == 'gcc' else 'which clang && which clang++', return_ = 'output', silent = True),
             gcc_install_prefix = Options.gcc_install_prefix,
             pybind11 = Options.pybind11,
             python_include = Options.python_include_dir,
             python_lib = Options.python_lib)
    )

    if False and Options.external_link:
        rosetta_cmake = """
            include_directories(SYSTEM {system_include})
            include_directories({rosetta_include})
            add_definitions({defs})
            link_directories(lib)

            set(PYROSETTA_EXTERNAL_LINK ON)
            """.format(
            system_include  = ' '.join( get_rosetta_system_include_directories() + [Options.pybind11] ),
            rosetta_include = ' '.join( get_rosetta_include_directories() ),
            defs = ' '.join([ '-D'+d for d in get_defines()])
        )
        cmake = cmake.replace('#%__Rosetta_cmake_instructions__%#', rosetta_cmake)
        cmake = cmake.replace('#%__PyRosetta_sources__%#', '\n'.join(extra_sources))
        cmake = cmake.replace('#%__Rosetta_libraries__%#', ' '.join(libs + ["${LINK_EXTERNAL_LIBS}"]))
        cmake = cmake.replace('#%__PyRosetta_build_config__%#', build_config)

        modified |= update_source_file(prefix + 'CMakeLists.txt', cmake)

    else:
        rosetta_cmake =  ''.join( ['include({}.cmake)\n'.format(l) for l in libs] )
        rosetta_cmake += '\ninclude_directories(SYSTEM {})\n\n'.format( ' '.join(get_rosetta_system_include_directories() + [Options.pybind11] ) )
        rosetta_cmake += '\ninclude_directories({})\n\n'.format( ' '.join( get_rosetta_include_directories() ) )
        rosetta_cmake += 'add_definitions({})\n'.format(' '.join([ '-D'+d for d in get_defines()] ) )

        cmake = cmake.replace('#%__Rosetta_cmake_instructions__%#', rosetta_cmake)
        cmake = cmake.replace('#%__PyRosetta_sources__%#', '\n'.join(extra_sources + ['$<TARGET_OBJECTS:{}>'.format(l) for l in libs] ) )  # cmake = cmake.replace('#%__PyRosetta_sources__%#', '\n'.join([ os.path.abspath(prefix + f) for f in extra_sources]))
        cmake = cmake.replace('#%__Rosetta_libraries__%#', '')  # cmake = cmake.replace('#%__Rosetta_libraries__%#', ' '.join(libs))
        cmake = cmake.replace('#%__PyRosetta_build_config__%#', build_config)
        cmake = cmake.replace('#%__PyRosetta_compile_options__%#', '-std=c++14' if Options.torch or Options.tensorflow else '')
        cmake = cmake.replace('#%__PyRosetta_target_link_libraries__%#',
                              ( 'c10 torch torch_cpu torch_global_deps' if Options.torch else '' )
                              + ( ' tensorflow' if Options.tensorflow else '' )
                              )

        modified |= update_source_file(prefix + 'CMakeLists.txt', cmake)

    return modified


def cmake_needs_to_be_run(rosetta_source_path):
    ''' Perform some basic sanity checks on build directory to see if cmake was run before. This should help for cases when build script was aborted during rebuilds.
    '''
    cmake_list = get_binding_build_root(rosetta_source_path, source=True) + '/CMakeLists.txt'
    build_ninja = get_binding_build_root(rosetta_source_path, build=True) + '/build.ninja'

    if ( not os.path.isfile(build_ninja) ) or \
       ( os.path.getmtime(cmake_list) >=  os.path.getmtime(build_ninja) ): return True
    else: return False


def run_cmake(rosetta_source_path):
    ''' Build generate bindings '''
    prefix = get_binding_build_root(rosetta_source_path, build=True) + '/'

    config = '-DCMAKE_BUILD_TYPE={}'.format(Options.type)
    config += ' -DPYROSETTA_STRIP_MODULE={build_type}'.format(build_type="TRUE" if Options.strip_module else "FALSE")
    config += get_cmake_compiler_options()

    # if we ever go back to use static libs for intermediates: fix this on Mac: -DCMAKE_RANLIB=/usr/bin/ranlib -DCMAKE_AR=/usr/bin/ar

    execute('Running CMake...', 'cd {prefix} && cmake -G Ninja {} -DPYROSETTA_PYTHON_VERSION={python_version}{py_lib}{py_include}{gcc_install_prefix} ../source'.format(config, prefix=prefix, python_version=_python_version_,
                                                                                                                                                                        py_lib=' -DPYTHON_LIBRARY='+Options.python_lib if Options.python_lib else '',
                                                                                                                                                                        py_include=' -DPYTHON_INCLUDE_DIR='+Options.python_include_dir if Options.python_include_dir else '',
                                                                                                                                                                        gcc_install_prefix=' -DGCC_INSTALL_PREFIX='+Options.gcc_install_prefix if Options.gcc_install_prefix else ''))
    sys.stdout.flush()


def generate_bindings(rosetta_source_path):
    ''' Generate bindings using binder tools and return list of source files '''
    setup_source_directory_links(rosetta_source_path)
    maybe_version = '--version {} '.format(Options.version) if Options.version else ''
    submodule_extras = "PyRosetta bcl"
    if Options.zmq:
        submodule_extras = submodule_extras + " zeromq"
    execute('Updating compilation submodules, version, options and residue-type-enum files...', 'cd {rosetta_source_path} && ./update_submodules.sh {submodule_extras} && ./version.py {maybe_version} && ./update_options.sh && ./update_ResidueType_enum_files.sh'.format(**vars()) )

    prefix = get_binding_build_root(rosetta_source_path, source=True) + '/'

    # generate include file that contains all headers
    skip_extensions = (".fwd.hh", ".impl.hh", ".py.hh")

    signature = hashlib.md5()
    def signature_update(s): signature.update( s.encode('ascii', errors='replace') )
    signature_update(_script_version_)

    signature_file_name = prefix + '.signature'
    if os.path.isfile(signature_file_name):
        with open(signature_file_name) as f: disk_signature = f.read()
    else: disk_signature = None


    all_includes, serialization_instantiation = [], []
    for path in 'ObjexxFCL utility numeric basic core protocols'.split()[:-Options.skip_namespaces]:
        for dir_name, _, files in os.walk(rosetta_source_path + '/src/' + path):
            for f in sorted(files):
                if not is_dir_banned(dir_name):
                    if f.endswith('.hh')  and not any(f.endswith(ex) for ex in skip_extensions):
                        header = dir_name[len(rosetta_source_path+'/src/'):] + '/' + f
                        if header not in _banned_headers_  and  not header.startswith('basic/options/keys/OptionKeys.cc.gen'):
                            #print(header)
                            all_includes.append(header)

                    elif f.endswith('.cc')  and  Options.serialization:
                        #print('{} {}'.format(dir_name, f))
                        with codecs.open(dir_name + '/' + f, encoding='utf-8', errors='replace') as fcc:
                            for line in fcc.read().split('\n'):
                                source = dir_name[len(rosetta_source_path+'/src/'):] + '/' + f
                                #if line.split('(')[0] in 'SAVE_AND_LOAD_SERIALIZABLE EXTERNAL_SAVE_AND_LOAD_SERIALIZABLE SAVE_AND_LOAD_AND_CONSTRUCT_SERIALIZABLE INSTANTIATE_FOR_OUTPUT_ARCHIVES INSTANTIATE_FOR_INPUT_ARCHIVES'.split(): serialization_instantiation.append(line)
                                if 'SAVE_AND_LOAD_SERIALIZABLE' in line  and  '_SAVE_AND_LOAD_SERIALIZABLE' not in line :
                                    serialization_instantiation.append( line.replace('SAVE_AND_LOAD_SERIALIZABLE', 'EXTERN_SAVE_AND_LOAD_SERIALIZABLE') + '  // file:{}\n'.format(source) )

    all_includes.sort()
    serialization_instantiation.sort()

    include = prefix + 'all_rosetta_includes.hh'
    with open(include, 'w') as fh:
        for i in all_includes: fh.write( '#include <{}>\n'.format(i) )
        for s in serialization_instantiation: fh.write(s)

    for i in all_includes: signature_update(i); signature_update( str( os.path.getmtime(rosetta_source_path+'/src/'+i) ) )
    for s in serialization_instantiation: signature_update(s)

    binder_config = list(Options.binder_config)
    config = ''
    for config_file in binder_config:
        config += '\n# ' + config_file + '\n'
        with open(config_file) as f: config += f.read()

    if 'clang' not in Options.compiler: config += open('rosetta.gcc.config').read()
    with open(prefix + 'rosetta.config', 'w') as f: f.write(config)
    signature_update(config)

    includes = ''.join( [' -isystem '+i for i in get_rosetta_system_include_directories()] ) + ''.join( [' -I'+i for i in get_rosetta_include_directories()] )
    defines  = ''.join( [' -D'+d for d in get_defines()] ) + ' -DPYROSETTA_BINDER'

    if Platform == 'macos':
        if False and tuple( map(int, platform.mac_ver()[0].split('.') ) ) < (11, 4):
            includes += ' -isystem /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1'
            includes += ' -isystem /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/usr/include'
        else:
            # includes += ' -isystem /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk'
            # includes += ' -isystem /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include'
            # includes += ' -isystem /Library/Developer/CommandLineTools/usr/lib/clang/12.0.5/include'
            # includes += ' -isystem /Library/Developer/CommandLineTools/usr/include/c++/v1'
            # includes += ' -isystem `xcrun --show-sdk-path`/usr/include/c++/v1'
            # includes += ' -isystem `xcrun --show-sdk-path`/usr/include'
            includes += ' -isysroot `xcrun --show-sdk-path`'

    if Options.binder_options:      include  += ' ' + Options.binder_options
    if Options.binder_llvm_options: includes += ' ' + Options.binder_llvm_options

    cpp_standard = 'c++14' if Options.torch else 'c++11'
    config='./rosetta.config'
    annotate=' --annotate-includes' if Options.annotate_includes else ''
    trace=' --trace' if Options.trace else ''

    binder_command_line_options = f'--config {config} --root-module rosetta --prefix {prefix}{annotate}{trace} {include} -- -std={cpp_standard} {includes} {defines}'
    signature_update(binder_command_line_options)

    signature = signature.hexdigest()
    if signature != disk_signature:
        execute('Generating bindings...', 'cd {prefix} && {binder} {binder_command_line_options}'.format(prefix=prefix, binder=Options.binder, binder_command_line_options=binder_command_line_options) )
        with open(signature_file_name, 'w') as f: f.write(signature)
    else: print('No changes in include files detected, skipping Binder run...')

    with open(prefix+'rosetta.sources') as f: sources = f.read().split()
    modified = generate_cmake_file(rosetta_source_path, sources)

    if modified or (signature != disk_signature) or cmake_needs_to_be_run(rosetta_source_path): run_cmake(rosetta_source_path)
    else: print('No changes in source files detected, skipping CMake run...')


def build_generated_bindings(rosetta_source_path):
    ''' Build generate bindings '''
    prefix = get_binding_build_root(rosetta_source_path, build=True) + '/'
    execute('Building...', 'cd {prefix} && ninja {jobs}'.format(prefix=prefix, jobs="-j{}".format(Options.jobs) if Options.jobs else ""))


def generate_version_file(rosetta_source_path, file_name):
    if Options.version: shutil.copy(Options.version, file_name)
    else: shutil.copy(rosetta_source_path + '/.version.json', file_name)

    # release_json = rosetta_source_path + './../.release.json'
    # moved to main/version.json: elif os.path.isfile(release_json): shutil.copy(release_json, file_name)
        # with open(rosetta_source_path + '/.version.json') as f: info = json.load(f)
        # with open(file_name, 'w') as f:
        #     json.dump(dict(version = info['ver'],
        #                    date    = info['commit_date'],
        #                    source = dict(main = info['ver']),
        #                    release = 'unknown', # all release info fields set to unknown
        #                    year    = 'unknown',
        #                    week    = 'unknown',
        #     ), f, sort_keys=True, indent=2)


_documentation_version_template_ = '''

.. _version-json-file:

Version
-------

.. code-block:: javascript

'''

def generate_documentation(rosetta_source_path, path, version):
    path = os.path.abspath(path)
    print('Creating PyRosetta documentation at: {}...'.format(path))

    if not os.path.isdir(path): os.makedirs(path)

    source_prefix = get_binding_build_root(rosetta_source_path, source=True)
    build_prefix = get_binding_build_root(rosetta_source_path, build=True)
    documentation_prefix = get_binding_build_root(rosetta_source_path, documentation=True)

    dir_util_module.copy_tree(rosetta_source_path + '/src/python/PyRosetta/documentation', documentation_prefix, update=False)

    version_file = path + '/version.json'
    generate_version_file(rosetta_source_path, version_file)

    with open(version_file) as f: info = json.load(f)

    index = documentation_prefix + '/source/index.rst'
    with open(index) as f: text = f.read()
    with open(index, 'w') as f:
        s = 'Version: {package} (for full version info see :ref:`version-json-file`)'.format(**info) if info['package'] else 'No explicit version info was provided during build phase, see :ref:`version-json-file` for version information of individual packages'.format(**info)
        f.write(s + '\n\n' + text)
        f.write(_documentation_version_template_)
        for line in open(version_file): f.write('    '+line)


    documentation_source_prefix = documentation_prefix + '/source'
    generator = rosetta_source_path + '/src/python/PyRosetta/binder/sphinx-doc-generator.py'

    python = sys.executable
    execute('Generating Sphinx files for PyRosetta Python code...', '{python} {generator} -o {documentation_source_prefix} --javascript-path {documentation_source_prefix}/_static --javascript-web-path _static/ {build_prefix}/pyrosetta'.format(**vars()))
    execute('Generating Sphinx files for PyRosetta Rosetta code...', '{python} {generator} -o {documentation_source_prefix} --javascript-path {documentation_source_prefix}/_static --javascript-web-path _static/ --root pyrosetta.rosetta {source_prefix}/rosetta.modules'.format(**vars()))  # --depth 1

    execute('Building documentation...', 'cd {documentation_prefix} && make clean && PATH={python_path}:$PATH PYTHONPATH={build_prefix}:$PYTHONPATH make html'.format(python_path=os.path.dirname(python), **vars()))

    css = documentation_prefix + '/build/html/_static/basic.css'
    if os.path.isfile(css):
        with open(css) as f: css_data = f.read()
        with open(css, 'w') as f: f.write( css_data.replace('white-space: nowrap;', '') )

    dir_util_module.copy_tree(documentation_prefix + '/build/html', path, update=False)

    print('Creating PyRosetta documentation at: {}... Done!'.format(path))


def create_package(rosetta_source_path, path):
    print('Creating Python package at: {}...'.format(path));  sys.stdout.flush()

    package_prefix = path + '/setup'
    if os.path.isdir(package_prefix): shutil.rmtree(package_prefix)
    os.makedirs(package_prefix)

    for f in 'self-test.py PyMOL-RosettaServer.py PyMOL-RosettaServer.python3.py PyMOL-Rosetta-relay-client.python2.py PyMOL-Rosetta-relay-client.python3.py'.split(): shutil.copy(rosetta_source_path + '/src/python/PyRosetta/src/' + f, path)

    for d in 'demo test'.split(): dir_util_module.copy_tree(rosetta_source_path + '/src/python/PyRosetta/src/' + d, path + '/' + d, update=False)

    dir_util_module.copy_tree(rosetta_source_path + '/scripts/PyRosetta/public', path + '/apps', update=False)

    build_prefix = get_binding_build_root(rosetta_source_path, build=True)

    for f in 'setup.py setup.cfg ez_setup.py'.split(): shutil.copy(build_prefix + '/' + f, package_prefix)
    #shutil.copy(rosetta_source_path + '/../LICENSE.PyRosetta.md', package_prefix)

    for d in ['pyrosetta', 'rosetta']:
        dir_util_module.copy_tree(build_prefix + '/' + d, package_prefix + '/' + d, update=False)

        #symlink = path + '/' + d
        #if os.path.islink(symlink): os.unlink(symlink)
        os.symlink('./setup/' + d, path + '/' + d)

        for dir_name, dirs, files in os.walk(package_prefix + '/' + d):
            #if '__pycache__' in dirs: shutil.rmtree(dir_name + '/__pycache__')
            for d in '__pycache__ .git'.split():
                if d in dirs:
                    shutil.rmtree(dir_name + '/' + d)

            for f in files:
                if f.endswith('.pyc') or f in ['.git']: os.remove(dir_name + '/' + f)

    for d in 'pyrosetta/database/additional_protocol_data'.split(): shutil.rmtree(package_prefix + '/' + d)

    generate_version_file(rosetta_source_path, path + '/version.json')


def create_wheel(rosetta_source_path, wheel_path):
    print("Creating Python wheel at: {}...".format(wheel_path));  sys.stdout.flush()

    build_prefix = get_binding_build_root(rosetta_source_path, build=True)

    execute('Building wheel...', 'cd {build_prefix} && {python} setup.py bdist_wheel -d {output_path}'.format(python=sys.executable, output_path=os.path.abspath(wheel_path), **vars()))


def main(args):
    ''' PyRosetta building script '''

    parser = argparse.ArgumentParser()
    parser.add_argument('-j', '--jobs', default=1, const=0, nargs="?", type=int, help="Number of processors to use on when building, use '-j' with no arguments to launch job-per-core. (default: 1) ")
    parser.add_argument('-s', '--skip-generation-phase', action="store_true", help="Assume that bindings code is already generaded and skipp the Binder call's")
    parser.add_argument('-d', '--skip-building-phase', action="store_true", help="Assume that bindings code is already generaded and skipp the Binder call's")
    parser.add_argument("--type", default='Release', choices=['Release', 'Debug', 'MinSizeRel', 'RelWithDebInfo'], help="Specify build type")

    parser.add_argument("--strip-module", dest="strip_module", action="store_true", help="Strip symbols from compiled modules to produce minimum file size in release mode. (default on)")
    parser.add_argument("--no-strip-module", dest="strip_module", action="store_false", help="Do not strip symbols from compiled modules to support tracebacks in release mode.")
    parser.set_defaults(strip_module = True)

    parser.add_argument('--compiler', default='clang', help='Compiler to use, defualt is clang')
    parser.add_argument('--binder', default='', help='Path to Binder tool. If none is given then download, build and install binder into main/source/build/prefix. Use "--binder-debug" to control which mode of binder (debug/release) is used.')
    parser.add_argument("--binder-debug", action="store_true", help="Run binder tool in debug mode (only relevant if no '--binder' option was specified)")
    parser.add_argument("--binder-config", action="append", default=["rosetta.config"], help="Binder config file. [Default='rosetta.config']")
    parser.add_argument("--print-build-root", action="store_true", help="Print path to where PyRosetta binaries will be located with given build options and exit. Use this option to automate package creation.")
    parser.add_argument('--cross-compile', action="store_true", help='Specify for cross-compile build')
    parser.add_argument('--pybind11', default='', help='Path to pybind11 source tree')
    parser.add_argument('--annotate-includes', action="store_true", help='Annotate includes in generated PyRosetta source files')
    parser.add_argument('--trace', action="store_true", help='Binder will add trace output to to generated PyRosetta source files')

    parser.add_argument('-p', '--create-package', default='', help='Create PyRosetta Python package at specified path (default is to skip creating package)')
    parser.add_argument('--create-wheel', default='', help='Create python wheel in the specified directory. (default is to skip creating wheel)')

    parser.add_argument('--python-include-dir', default=None, help='Path to python C headers. Use this if CMake fails to autodetect it')
    parser.add_argument('--python-lib', default=None, help='Path to python library. Use this if CMake fails to autodetect it')

    parser.add_argument('--gcc-install-prefix', default=None, help='Path to GCC install prefix which will be used to determent location of libstdc++ for Binder build. Default is: auto-detected. Use this option if you would like to build Binder with compiler that was side-installed and which LLVM build system failed to identify. To see what path Binder uses for libstdc++ run `binder -- -xc++ -E -v`.')

    parser.add_argument('--serialization', action="store_true", help="Build PyRosetta with serialization enabled (off by default)")
    parser.add_argument('--multi-threaded', action="store_true", help="Build PyRosetta with multi_threaded enabled (off by default)")
    parser.add_argument('--torch', action="store_true", help="Build PyRosetta with lib torch support enabled (off by default)")
    parser.add_argument('--tensorflow', action="store_true", help="Build PyRosetta with lib tensorflow support enabled (off by default)")
    #parser.add_argument('--hdf5', action="store_true", help="Build PyRosetta with HDF5 enabled (off by default)")

    parser.add_argument('--binder-options', default=None, help='Specify Binder extra (non LLVM) options. Use this to specify options specific to Binder.')
    parser.add_argument('--binder-llvm-options', default=None, help='Specify Binder extra (LLVM) options. Use this to point Binder to additional include dir or add/change LLVM flags. For example on Mac with no Xcode install set this to "-isystem /Library/Developer/CommandLineTools/usr/include/c++/v1" to point Binder to correct location of system includes.')

    parser.add_argument('--documentation', help='Generate PyRosetta documentation at specified path (default is to skip documentation creation). Note that Sphinx package need to be installed for this to work correctly.')

    parser.set_defaults(zmq = True)
    parser.add_argument('--no-zmq', dest='zmq', action="store_false", help='Disable building and linking of ZeroMQ library')

    parser.add_argument('--version', help='Supply JSON version file to be used for during package creation and documentation building. File must be in the same format as standard Rosetta .release.json used to mark release versions. If no file is supplied script will fallback to use main/.version.json.')

    parser.add_argument('-n', '--skip-namespaces', default=-16777216, type=int, help="EXPERIMENTAL: Specify number of high-level Rosetta namespaces to skip during generation phase. This allow one to bypass bindings generations for higher level libraries (like protocols, core etc )Default is 0, - do not skip any namespaces.")


    #parser.add_argument('--build-suffix', default=None, help='Specify build suffix that will be be used when creating build directories. Default is None, - use either $HOSTNAME or value provided in local .hostname file.')


    global Options
    Options = parser.parse_args()

    #Options.build_suffix =  _machine_name_ if Options.build_suffix is None else Options.build_suffix

    rosetta_source_path = os.path.abspath('./../../../')

    binding_build_root = get_binding_build_root(rosetta_source_path)

    if Options.print_build_root: print(binding_build_root, end=('\n' if sys.stdout.isatty() else '') ); sys.exit(0)

    print('Creating PyRosetta in "{}" mode in: {}'.format(Options.type, binding_build_root))

    if len(Options.binder_config) > 1 : print('Binder config: ', Options.binder_config)

    link_supplemental_files(rosetta_source_path)

    if Options.skip_generation_phase:
        print('Option --skip-generation-phase is supplied, skipping generation phase...')

    else:
        if (not Options.pybind11) or (not Options.binder):
            res, output = execute('Getting main repository root...', 'cd {rosetta_source_path} && git rev-parse --show-toplevel'.format(**vars()), return_='tuple')
            main_repository_root = output.split()[0] # removing \n at the end

            if res == 0:
                if Options.pybind11: print('NOTE: options `--pybind11` was specified, skipping Pybind11 submodule update...')
                else: execute('Updating Pybind11 Git submodule...', 'cd {}/.. && git submodule update --init -- source/external/pybind11'.format(rosetta_source_path) )  # --recursive

                if Options.binder: print('NOTE: options `--binder` was specified, skipping Binder submodule update...')
                else: execute('Updating Binder Git submodule...', 'cd {}/.. && git submodule update --init -- source/src/python/PyRosetta/binder'.format(rosetta_source_path) )

                if main_repository_root == os.path.abspath(rosetta_source_path + '/../'):
                    output = execute('Checking if Binder submodule present...',  'cd {}/.. && git submodule status'.format(rosetta_source_path), return_='output', silent=True)
                    if 'source/src/python/PyRosetta/binder' not in output: print('ERROR: Binder submodule is not found... terminating...'); sys.exit(1)

        if not Options.pybind11: Options.pybind11 = os.path.abspath(rosetta_source_path + '/external/pybind11/include')  # install_pybind11(rosetta_source_path + '/build/prefix')
        if not Options.binder: Options.binder = install_llvm_tool('binder', rosetta_source_path+'/src/python/PyRosetta/binder/source', rosetta_source_path + '/build/prefix', debug=Options.binder_debug, compiler=Options.compiler, jobs=Options.jobs, gcc_install_prefix=Options.gcc_install_prefix)

        generate_bindings(rosetta_source_path)


    if Options.skip_building_phase: print('Option --skip-building-phase is supplied, skipping building phase...')
    else: build_generated_bindings(rosetta_source_path)

    if Options.documentation: generate_documentation(rosetta_source_path, Options.documentation, Options.version)
    if Options.create_package: create_package(rosetta_source_path, Options.create_package)
    if Options.create_wheel: create_wheel(rosetta_source_path, Options.create_wheel)


if __name__ == "__main__":

    # Check if current working dir is where this file is located...
    if os.path.dirname(os.path.realpath(__file__)) != os.getcwd():
        print('This script must be run from within source/src/python/PyRosetta/ directory! Exiting...')
        sys.exit(1)

    main(sys.argv)
