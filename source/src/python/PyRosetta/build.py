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

import os, sys, argparse, platform, subprocess, imp, shutil, codecs, distutils.dir_util

from collections import OrderedDict

if sys.platform.startswith("linux"): Platform = "linux" # can be linux1, linux2, etc
elif sys.platform == "darwin" : Platform = "macos"
elif sys.platform == "cygwin" : Platform = "cygwin"
elif sys.platform == "win32" : Platform = "windows"
else: Platform = "unknown"
PlatformBits = platform.architecture()[0][:2]

_machine_name_ = os.uname()[1]

_python_version_ = '{}.{}'.format(sys.version_info.major, sys.version_info.minor)  # should be formatted: 2.7 or 3.5
#_python_version_ = '{}.{}.{}'.format(sys.version_info.major, sys.version_info.minor, sys.version_info.micro)  # should be formatted: 2.7.6 or 3.5.0

_pybind11_version_ = '1ee4128cfe8efcab618a980649cad9d830d8b32b'

_banned_dirs_ = 'src/utility/pointer src/protocols/jd3'.split()  # src/utility/keys src/utility/options src/basic/options
_banned_headers_ = 'utility/py/PyHelper.hh utility/keys/KeyCount.hh utility/keys/KeyLookup.functors.hh'
_banned_headers_ +=' core/scoring/fiber_diffraction/FiberDiffractionKernelGpu.hh' # GPU code

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
    r = 'external external/include external/boost_1_55_0 external/dbio external/dbio/sqlite3 external/libxml2/include'.split()
    return r

def get_rosetta_include_directories():
    ''' return list of include directories for compilation '''
    r = 'src'.split()
    r.append('src/platform/'+Platform)
    return r


def get_defines():
    ''' return list of #defines '''
    defines = 'PYROSETTA BOOST_ERROR_CODE_HEADER_ONLY BOOST_SYSTEM_NO_DEPRECATED BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS PTR_STD'
    if Platform == 'macos': defines += ' UNUSUAL_ALLOCATOR_DECLARATION'
    if Options.type in 'Release MinSizeRel': defines += ' NDEBUG'
    return defines.split()


def execute(message, command_line, return_='status', until_successes=False, terminate_on_failure=True, silent=False):
    print(message);  print(command_line); sys.stdout.flush();
    while True:

        p = subprocess.Popen(command_line, bufsize=0, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, errors = p.communicate()

        output = output + errors

        output = output.decode(encoding="utf-8", errors="replace")

        exit_code = p.returncode

        if exit_code  or  not silent: print(output); sys.stdout.flush();

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
    ''' write data to a file but only if file does not exist or it content different from data '''
    if( not os.path.isfile(file_name)  or  open(file_name).read() != data ):
        print('Writing '+ file_name)
        with open(file_name, 'w') as f: f.write(data)


def get_compiler_family():
    ''' Try to guess compiler family using Options.compiler '''
    if 'clang' in Options.compiler: return 'clang'
    if 'gcc'   in Options.compiler: return 'gcc'
    if 'cl'    in Options.compiler: return 'cl'
    return 'unknown'


def install_llvm_tool(name, source_location, prefix, debug, clean=True):
    ''' Install and update (if needed) custom LLVM tool at given prefix (from config).
        Return absolute path to executable on success and terminate with error on failure
    '''
    release = 'release_38'
    prefix += '/llvm-3.8'

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

    for line in open(cmake_lists):
        if line == tool_build_line: break
    else:
        with open(cmake_lists, 'w') as f: f.write( open(cmake_lists).read() + tool_build_line + '\n' )

    build_dir = prefix+'/llvm/build_' + release + '.' + Platform + '.' +_machine_name_ + ('.debug' if debug else '.release')
    if not os.path.isdir(build_dir): os.makedirs(build_dir)
    execute('Building tool: {}...'.format(name), 'cd {build_dir} && cmake -G Ninja -DCMAKE_BUILD_TYPE={build_type} -DLLVM_ENABLE_EH=1 -DLLVM_ENABLE_RTTI=ON {gcc_install_prefix} .. && ninja -j{jobs}'.format(build_dir=build_dir, jobs=Options.jobs, build_type='Debug' if debug else 'Release',
                                                                                                                                                                                                              gcc_install_prefix='-DGCC_INSTALL_PREFIX='+Options.gcc_install_prefix if Options.gcc_install_prefix else ''),
            silent=True)
    print()
    # build_dir = prefix+'/llvm/build-ninja-' + release
    # if not os.path.isdir(build_dir): os.makedirs(build_dir)
    # execute('Building tool: {}...'.format(name), 'cd {build_dir} && cmake -DCMAKE_BUILD_TYPE={build_type} .. -G Ninja && ninja -j{jobs}'.format(build_dir=build_dir, jobs=Options.jobs, build_type='Debug' if debug else 'Release')) )

    executable = build_dir + '/bin/' + name
    if not os.path.isfile(executable): print("\nEncounter error while running install_llvm_tool: Build is complete but executable {} is not there!!!".format(executable) ); sys.exit(1)

    return executable


def install_pybind11(prefix, clean=True):
    ''' Download and install PyBind11 library at given prefix. Install version specified by _pybind11_version_ sha1
    '''
    #git_checkout = '( git fetch && git checkout {0} && git reset --hard {0} && git pull )'.format(_pybind11_version_) if clean else 'git checkout {}'.format(_pybind11_version_)
    git_checkout = '( git fetch && git reset --hard {0} )'.format(_pybind11_version_) if clean else 'git checkout {}'.format(_pybind11_version_)

    if not os.path.isdir(prefix): os.makedirs(prefix)
    package_dir = prefix + '/pybind11'

    if not os.path.isdir(package_dir): execute('Clonning pybind11...', 'cd {} && git clone https://github.com/RosettaCommons/pybind11.git'.format(prefix) )
    execute('Checking out PyBind11 revision: {}...'.format(_pybind11_version_), 'cd {package_dir} && ( {git_checkout} )'.format(package_dir=package_dir, git_checkout=git_checkout), silent=True)
    print()

    include = package_dir + '/include/pybind11/pybind11.h'
    if not os.path.isfile(include): print("\nEncounter error while running install_pybind11: Install is complete but include file {} is not there!!!".format(include) ); sys.exit(1)

    return package_dir + '/include'


def get_binding_build_root(rosetta_source_path, source=False, build=False):
    ''' Calculate bindings build path using current platform and compilation settings and create dir if it is not there yet '''

    p = os.path.join(rosetta_source_path, 'build/PyRosetta')

    #p = os.path.join(p, 'cross_compile' if Options.cross_compile else (Platform+ '/' + options.compiler) )
    p =  os.path.join(p, Platform+ '/' + get_compiler_family() + '/python-' + _python_version_)

    #p = os.path.join(p, 'monolith' if True else 'namespace' )

    p = os.path.join(p, Options.type.lower())

    source_p = p + '/source'
    build_p  = p + '/build'

    if not os.path.isdir(source_p): os.makedirs(source_p)
    if not os.path.isdir(build_p) : os.makedirs(build_p)

    if source: return source_p
    if build:  return build_p

    return p



def copy_supplemental_files(rosetta_source_path):
    prefix = get_binding_build_root(rosetta_source_path, build=True)
    source = rosetta_source_path + '/src/python/PyRosetta/src'

    distutils.dir_util.copy_tree(source, prefix, update=False)

    if Platform not in ['windows', 'cygwin'] and (not os.path.islink(prefix + '/database')): os.symlink('../../../../../../../../database', prefix + '/database')  # creating link to Rosetta database dir


def generate_rosetta_external_cmake_files(rosetta_source_path, prefix):
    libs = OrderedDict([ ('cppdb', ['dbio/cppdb/atomic_counter', 'dbio/cppdb/conn_manager', 'dbio/cppdb/driver_manager', 'dbio/cppdb/frontend', 'dbio/cppdb/backend',
                                    'dbio/cppdb/mutex', 'dbio/cppdb/pool', 'dbio/cppdb/shared_object', 'dbio/cppdb/sqlite3_backend', 'dbio/cppdb/utils'] ),
                         ('sqlite3', ['dbio/sqlite3/sqlite3'] ), ])

    defines = dict(cppdb = 'CPPDB_EXPORTS CPPDB_DISABLE_SHARED_OBJECT_LOADING CPPDB_DISABLE_THREAD_SAFETY CPPDB_WITH_SQLITE3' +
                           ' CPPDB_LIBRARY_PREFIX=\\"lib\\" CPPDB_LIBRARY_SUFFIX=\\".dylib\\" CPPDB_SOVERSION=\\"0\\" ' +
                           ' CPPDB_MAJOR=0 CPPDB_MINOR=3 CPPDB_PATCH=0 CPPDB_VERSION=\\"0.3.0\\"',
                   sqlite3 = 'SQLITE_DISABLE_LFS SQLITE_OMIT_LOAD_EXTENSION SQLITE_THREADSAFE=0')

    scons_file_extension = '.external.settings'
    external_scons_files = [f for f in os.listdir(rosetta_source_path+'/external') if f.endswith(scons_file_extension)]
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
            defines[lib_name] = ' '.join( G['defines'] )

    for l in libs:
        t  = 'add_library({} OBJECT\n{})\n\n'.format(l, '\n'.join( [ rosetta_source_path + '/external/' + s for s in libs[l]] ))
        t += 'set_property(TARGET {} PROPERTY POSITION_INDEPENDENT_CODE ON)\n'.format(l)
        if defines[l]: t += 'target_compile_options({} PRIVATE {})\n'.format(l, ' '.join( ['-D'+d for d in defines[l].split()] ) )   #  target_compile_definitions
        #t += 'target_compile_options({} PUBLIC -fPIC {})\n'.format(l, ' '.join([ '-D'+d for d in defines[l].split() ] ) )   #  target_compile_definitions
        update_source_file(prefix + l + '.cmake', t)

    return list( libs.keys() )


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
        t += 'set_property(TARGET {} PROPERTY POSITION_INDEPENDENT_CODE ON)\n'.format(lib)
        #t += '\ntarget_compile_options({} PUBLIC -fPIC)\n'.format(lib)  # Enable Position Independent Code generation for libraries
        update_source_file(prefix + lib + '.cmake', t)

        libs.append(lib)

    return libs



def generate_cmake_file(rosetta_source_path, extra_sources):
    ''' generate cmake file in bindings source location '''

    prefix = get_binding_build_root(rosetta_source_path, source=True) + '/'

    libs = generate_rosetta_cmake_files(rosetta_source_path, prefix) + generate_rosetta_external_cmake_files(rosetta_source_path, prefix)

    rosetta_cmake =  ''.join( ['include({}.cmake)\n'.format(l) for l in libs] )
    rosetta_cmake += '\ninclude_directories(SYSTEM {})\n\n'.format( ' '.join(get_rosetta_system_include_directories() ) )
    rosetta_cmake += '\ninclude_directories({})\n\n'.format( ' '.join( get_rosetta_include_directories() + [Options.pybind11] ) )
    rosetta_cmake += 'add_definitions({})\n'.format(' '.join([ '-D'+d for d in get_defines()] ) )

    cmake = open('cmake.template').read()

    cmake = cmake.replace('#%__Rosetta_cmake_instructions__%#', rosetta_cmake)
    cmake = cmake.replace( '#%__PyRosetta_sources__%#', '\n'.join(extra_sources + ['$<TARGET_OBJECTS:{}>'.format(l) for l in libs] ) )  # cmake = cmake.replace('#%__PyRosetta_sources__%#', '\n'.join([ os.path.abspath(prefix + f) for f in extra_sources]))
    cmake = cmake.replace('#%__Rosetta_libraries__%#', '')  # cmake = cmake.replace('#%__Rosetta_libraries__%#', ' '.join(libs))

    update_source_file(prefix + 'CMakeLists.txt', cmake)


def generate_bindings(rosetta_source_path):
    ''' Generate bindings using binder tools and return list of source files '''
    copy_supplemental_files(rosetta_source_path)
    execute('Updating version, options and residue-type-enum files...', 'cd {} && ./version.py && ./update_options.sh && ./update_ResidueType_enum_files.sh'.format(rosetta_source_path) )

    prefix = get_binding_build_root(rosetta_source_path, source=True) + '/'

    for d in ['src', 'external']:
        s = prefix + d
        if os.path.islink(s): os.unlink(s)
        os.symlink(rosetta_source_path + '/' + d, s)



    # generate include file that contains all headers
    all_includes = []

    #for path in 'ObjexxFCL utility numeric basic core protocols'.split():
    for path in 'ObjexxFCL utility numeric basic core protocols'.split():
        for dir_name, _, files in os.walk(rosetta_source_path + '/src/' + path):
            for f in sorted(files):
                if f.endswith('.hh')  and  (not f.endswith('.fwd.hh')):
                    #if f == 'exit.hh':
                    #if dir_name.endswith('utility'):
                    if not is_dir_banned(dir_name):
                    #if 'utility/' in dir_name  and  'src/utility/pointer' not in dir_name:
                        header = dir_name[len(rosetta_source_path+'/src/'):] + '/' + f
                        if header not in _banned_headers_  and  not header.startswith('basic/options/keys/OptionKeys.cc.gen'):
                            #print(header)
                            all_includes.append(header)

    all_includes.sort()
    include = prefix + 'all_rosetta_includes.hh'
    with open(include, 'w') as fh:
        for i in all_includes: fh.write( '#include <{}>\n'.format(i) )

    config = open('rosetta.config').read()
    if 'clang' not in Options.compiler: config += open('rosetta.gcc.config').read()
    with open(prefix + 'rosetta.config', 'w') as f: f.write(config)

    includes = ''.join( [' -isystem '+i for i in get_rosetta_system_include_directories()] ) + ''.join( [' -I'+i for i in get_rosetta_include_directories()] )
    defines  = ''.join( [' -D'+d for d in get_defines()] )

    if Platform == 'macos': includes = '-isystem/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/../include/c++/v1' + includes

    execute('Generating bindings...', 'cd {prefix} && {} --config {config} --root-module rosetta --prefix {prefix}{annotate}{trace} {} -- -std=c++11 {} {}'.format(Options.binder, include, includes, defines,
                                                                                                                                                            prefix=prefix, config='./rosetta.config',
                                                                                                                                                            annotate=' --annotate-includes' if Options.annotate_includes else '',
                                                                                                                                                            trace=' --trace' if Options.trace else '',) ) # -stdlib=libc++ -x c++

    sources = open(prefix+'rosetta.sources').read().split()

    generate_cmake_file(rosetta_source_path, sources)



def  build_generated_bindings(rosetta_source_path):
    ''' Build generate bindings '''
    prefix = get_binding_build_root(rosetta_source_path, build=True) + '/'

    config = '-DCMAKE_BUILD_TYPE={}'.format(Options.type)

    if Platform == "linux": config += ' -DCMAKE_C_COMPILER=`which clang` -DCMAKE_CXX_COMPILER=`which clang++`'if Options.compiler == 'clang' else ' -DCMAKE_C_COMPILER=`which gcc` -DCMAKE_CXX_COMPILER=`which g++`'

    # if we ever go back to use static libs for intermediates: fix this on Mac: -DCMAKE_RANLIB=/usr/bin/ranlib -DCMAKE_AR=/usr/bin/ar

    execute('Running CMake...', 'cd {prefix} && cmake -G Ninja {} -DPYROSETTA_PYTHON_VERSION={python_version} {py_lib} {py_include} ../source'.format(config, prefix=prefix, python_version=_python_version_,
                                                                                                                                                      py_lib='-DPYTHON_LIBRARY='+Options.python_lib if Options.python_lib else '',
                                                                                                                                                      py_include='-DPYTHON_INCLUDE_DIR='+Options.python_include_dir if Options.python_include_dir else ''))

    execute('Building...', 'cd {prefix} && ninja -j{jobs}'.format(prefix=prefix, jobs=Options.jobs))


def generate_documentation(rosetta_source_path, path):
    path = os.path.abspath(path)
    print('Creating PyRosetta documentation at: {}...'.format(path))

    source_prefix = get_binding_build_root(rosetta_source_path, source=True)
    build_prefix = get_binding_build_root(rosetta_source_path, build=True)

    pyrosetta_modules = ['pyrosetta']
    for dir_name, _, files in os.walk(build_prefix+'/pyrosetta'):
        for f in sorted(files):
            if f.endswith('.py')  and  (not f.startswith('__')):
                pyrosetta_modules.append( (dir_name[len(build_prefix)+1:] + '/' + f).replace('/', '.')[:-len('.py')] )

    with open(source_prefix + '/rosetta.modules') as fh: modules = ' '.join( pyrosetta_modules + ['rosetta'] + [ 'rosetta.'+ f for f in fh.read().split()] )

    if not os.path.isdir(path): os.makedirs(path)
    execute('Generating PyRosetta documentation...', 'cd {build_prefix} && {pydoc} -w {modules} && mv *.html {path}'.format(pydoc=Options.pydoc, **vars()), silent=True)

    with open(path+'/index.html', 'w') as f: f.write('<!DOCTYPE html><html><head><title>PyRosetta-4 Documentation</title></head> <body><a href="rosetta.html">[Rosetta Module Documentation]</a> <a href="pyrosetta.html">[PyRosetta Module Documentation]</a> </body> </html>')


def create_package(rosetta_source_path, path):
    print('Creating Python package at: {}...'.format(path))

    package_prefix = path + '/setup'
    if not os.path.isdir(package_prefix): os.makedirs(package_prefix)

    for f in 'self-test.py PyMol-RosettaServer.py'.split(): shutil.copy(rosetta_source_path + '/src/python/PyRosetta/src/' + f, path)

    for d in 'demo test'.split(): distutils.dir_util.copy_tree(rosetta_source_path + '/src/python/PyRosetta/src/' + d, path + '/' + d, update=False)

    distutils.dir_util.copy_tree(rosetta_source_path + '/../database', package_prefix + '/database', update=False)
    distutils.dir_util.copy_tree(rosetta_source_path + '/src/python/PyRosetta/package', package_prefix, update=False)

    build_prefix = get_binding_build_root(rosetta_source_path, build=True)
    shutil.copy(build_prefix + '/rosetta.so', package_prefix)
    distutils.dir_util.copy_tree(build_prefix + '/pyrosetta', package_prefix + '/pyrosetta', update=False)

    #generate_documentation(rosetta_source_path, path+'/documentation')



def main(args):
    ''' PyRosetta building script '''

    parser = argparse.ArgumentParser()

    parser.add_argument('-j', '--jobs', default=1, type=int, help="Number of processors to use on when building. (default: use all avaliable memory)")
    parser.add_argument('-s', '--skip-generation-phase', action="store_true", help="Assume that bindings code is already generaded and skipp the Binder call's")
    parser.add_argument('-d', '--skip-building-phase', action="store_true", help="Assume that bindings code is already generaded and skipp the Binder call's")
    parser.add_argument("--type", default='Release', choices=['Release', 'Debug', 'MinSizeRel', 'RelWithDebInfo'], help="Specify build type")
    parser.add_argument('--compiler', default='clang', help='Compiler to use, defualt is clang')
    parser.add_argument('--binder', default='', help='Path to Binder tool. If none is given then download, build and install binder into main/source/build/prefix. Use "--binder-debug" to control which mode of binder (debug/release) is used.')
    parser.add_argument("--binder-debug", action="store_true", help="Run binder tool in debug mode (only relevant if no '--binder' option was specified)")
    parser.add_argument("--print-build-root", action="store_true", help="Print path to where PyRosetta binaries will be located with given options and exit. Use this option to automate package creation.")
    parser.add_argument('--cross-compile', action="store_true", help='Specify for cross-compile build')
    parser.add_argument('--pybind11', default='', help='Path to pybind11 source tree')
    parser.add_argument('--annotate-includes', action="store_true", help='Annotate includes in generated PyRosetta source files')
    parser.add_argument('--trace', action="store_true", help='Binder will add trace output to to generated PyRosetta source files')

    parser.add_argument("--pydoc", default='pydoc', help="Specify pydoc executable to use (default is 'pydoc')")
    parser.add_argument('--documentation', default='', help='Generate PyRosetta documentation at specified path (default is to skip documentation creation)')

    parser.add_argument('--create-package', default='', help='Create PyRosetta Python package at specified path (default is to skip creating package)')

    parser.add_argument('--python-include-dir', default=None, help='Path to python C headers. Use this if CMake fails to autodetect it')
    parser.add_argument('--python-lib', default=None, help='Path to python library. Use this if CMake fails to autodetect it')

    parser.add_argument('--gcc-install-prefix', default=None, help='Path to GCC install prefix which will be used to determent location of libstdc++ for Binder build. Default is: auto-detected. Use this option if you would like to build Binder with compiler that was side-installed and which LLVM build system failed to identify. To see what path Binder uses for libstdc++ run `binder -- -xc++ -E -v`.')

    global Options
    Options = parser.parse_args()

    rosetta_source_path = os.path.abspath('./../../../')

    binding_build_root = get_binding_build_root(rosetta_source_path)

    if Options.print_build_root: print(binding_build_root, end=''); sys.exit(0)

    print('Creating PyRosetta in "{}" mode in: {}'.format(Options.type, binding_build_root))

    if not Options.binder: Options.binder = install_llvm_tool('binder', rosetta_source_path+'/src/python/PyRosetta/binder', rosetta_source_path + '/build/prefix', Options.binder_debug)

    if not Options.pybind11: Options.pybind11 = install_pybind11(rosetta_source_path + '/build/prefix')

    if Options.skip_generation_phase: print('Option --skip-generation-phase is supplied, skipping generation phase...')
    else: generate_bindings(rosetta_source_path)

    if Options.skip_building_phase: print('Option --skip-building-phase is supplied, skipping building phase...')
    else: build_generated_bindings(rosetta_source_path)

    if Options.documentation: generate_documentation(rosetta_source_path, Options.documentation)
    if Options.create_package: create_package(rosetta_source_path, Options.create_package)


if __name__ == "__main__":

    # Check if current working dir is where this file is located...
    if os.path.dirname(os.path.realpath(__file__)) != os.getcwd():
        print('This script must be run from within setup directory! Exiting...')
        sys.exit(1)

    main(sys.argv)
