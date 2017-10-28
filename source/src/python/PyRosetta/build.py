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

import os, sys, argparse, platform, subprocess, imp, shutil, codecs, distutils.dir_util, json

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

_pybind11_version_ = 'fa6a4241326a361fc23915f6a82c1e48667de668'

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
    if Options.serialization: defines += ' SERIALIZATION'
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
    release = 'release_40'
    prefix += '/llvm-4.0'

    build_dir = prefix+'/build_' + release + '.' + Platform + '.' +_machine_name_ + ('.debug' if debug else '.release')

    binder_head = execute('Getting binder HEAD commit SHA1...', 'cd {} && git rev-parse HEAD'.format(source_location), return_='output', silent=True).split('\n')[0]

    signature = dict(config = 'LLVM install by install_llvm_tool version: 1.0', binder = binder_head)
    signature_file_name = build_dir + '/.signature.json'

    disk_signature = dict(config = 'unknown', binder = 'unknown')
    if os.path.isfile(signature_file_name):
        with open(signature_file_name) as f: disk_signature = json.load(f)

    if signature == disk_signature:
        print('LLVM:{} + Binder install is detected at {}, skipping LLVM installation Binder building procedures...\n'.format(release, build_dir))

    else:
        print('LLVM build detected, but config/binder version has changed, perfoming a clean rebuild...')
        if os.path.isdir(build_dir): shutil.rmtree(build_dir)

        # if signature['config'] != disk_signature['config']:
        #     print( 'LLVM build detected, but config version mismatched: was:"{}" current:"{}", perfoming a clean rebuild...'.format(disk_signature['config'], signature['config']) )
        #     if os.path.isdir(build_dir): shutil.rmtree(build_dir)
        # else: print( 'Binder build detected, but source version mismatched: was:{} current:{}, rebuilding...'.format(disk_signature['binder'], signature['binder']) )

        if not os.path.isdir(build_dir): os.makedirs(build_dir)

        git_checkout = '( git checkout {0} && git reset --hard {0} )'.format(release) if clean else 'git checkout {}'.format(release)

        if os.path.isdir(prefix) and  (not os.path.isdir(prefix+'/.git')): shutil.rmtree(prefix)  # removing old style checkoiut

        if not os.path.isdir(prefix):
            print( 'No LLVM:{} + Binder install is detected! Going to check out LLVM and install Binder. This procedure will require ~3Gb of free disk space and will only be needed to be done once...\n'.format(release) )
            os.makedirs(prefix)

        if not os.path.isdir(prefix+'/.git'): execute('Clonning llvm...', 'cd {} && git clone https://github.com/llvm-mirror/llvm.git .'.format(prefix) )
        execute('Checking out LLVM revision: {}...'.format(release), 'cd {prefix} && ( {git_checkout} || ( git fetch && {git_checkout} ) )'.format(prefix=prefix, git_checkout=git_checkout) )

        if not os.path.isdir(prefix+'/tools/clang'): execute('Clonning clang...', 'cd {}/tools && git clone https://github.com/llvm-mirror/clang.git clang'.format(prefix) )
        execute('Checking out Clang revision: {}...'.format(release), 'cd {prefix}/tools/clang && ( {git_checkout} || ( git fetch && {git_checkout} ) )'.format(prefix=prefix, git_checkout=git_checkout) )

        if not os.path.isdir(prefix+'/tools/clang/tools/extra'): os.makedirs(prefix+'/tools/clang/tools/extra')

        tool_link_path = '{prefix}/tools/clang/tools/extra/{name}'.format(prefix=prefix, name=name)
        if os.path.islink(tool_link_path): os.unlink(tool_link_path)
        os.symlink(source_location, tool_link_path)

        cmake_lists = prefix + '/tools/clang/tools/extra/CMakeLists.txt'
        tool_build_line = 'add_subdirectory({})'.format(name)

        if not os.path.isfile(cmake_lists):
            with open(cmake_lists, 'w') as f: f.write(tool_build_line + '\n')


        config = '-DCMAKE_BUILD_TYPE={build_type}'.format(build_type='Debug' if debug else 'Release')
        if Platform == "linux": config += ' -DCMAKE_C_COMPILER=`which clang` -DCMAKE_CXX_COMPILER=`which clang++`'if Options.compiler == 'clang' else ' -DCMAKE_C_COMPILER=`which gcc` -DCMAKE_CXX_COMPILER=`which g++`'

        execute(
            'Building tool: {}...'.format(name),
            'cd {build_dir} && cmake -G Ninja {config} -DLLVM_ENABLE_EH=1 -DLLVM_ENABLE_RTTI=ON {gcc_install_prefix} .. && ninja {jobs}'.format(
                build_dir=build_dir, config=config,
                jobs="-j{}".format(Options.jobs) if Options.jobs else "",
                gcc_install_prefix='-DGCC_INSTALL_PREFIX='+Options.gcc_install_prefix if Options.gcc_install_prefix else ''),
            silence_output=True)
        print()
        # build_dir = prefix+'/build-ninja-' + release
        # if not os.path.isdir(build_dir): os.makedirs(build_dir)
        # execute('Building tool: {}...'.format(name), 'cd {build_dir} && cmake -DCMAKE_BUILD_TYPE={build_type} .. -G Ninja && ninja -j{jobs}'.format(build_dir=build_dir, jobs=Options.jobs, build_type='Debug' if debug else 'Release')) )

        with open(signature_file_name, 'w') as f: json.dump(signature, f, sort_keys=True, indent=2)

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


def get_binding_build_root(rosetta_source_path, source=False, build=False, documentation=False):
    ''' Calculate bindings build path using current platform and compilation settings and create dir if it is not there yet '''

    p = os.path.join(rosetta_source_path, 'build/PyRosetta')

    p =  os.path.join(p, Platform + '/' + get_compiler_family() + '/python-' + _python_version_)

    p = os.path.join(p, Options.type.lower() + ('.serialization' if Options.serialization else '') )

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


def copy_supplemental_files(rosetta_source_path):
    prefix = get_binding_build_root(rosetta_source_path, build=True)
    source = rosetta_source_path + '/src/python/PyRosetta/src'

    distutils.dir_util.copy_tree(source, prefix, update=False)

    database_dest = prefix + '/pyrosetta/database'
    if os.path.islink(database_dest): os.unlink(database_dest)
    elif os.path.isdir(database_dest): shutil.rmtree(database_dest)
    os.symlink('../../../../../../../../../database', database_dest)

    if not os.path.islink(prefix + '/apps'): os.symlink('../../../../../../../scripts/PyRosetta/public', prefix + '/apps')  # creating link to PyRosetta apps dir


def setup_source_directory_links(rosetta_source_path):
    prefix = get_binding_build_root(rosetta_source_path, source=True)

    for d in ['src', 'external']:
        source_path = os.path.relpath(os.path.join(rosetta_source_path, d), prefix)
        s = os.path.join(prefix, d)
        if os.path.islink(s): os.unlink(s)
        os.symlink(source_path, s)

    if Options.external_link:
        target_lib_path = os.path.relpath(
            os.path.join(rosetta_source_path, "cmake", "build_{}".format(Options.external_link)), prefix)

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
        #t += 'set_property(TARGET {} PROPERTY POSITION_INDEPENDENT_CODE ON)\n'.format(lib)
        t += 'set_target_properties({} PROPERTIES POSITION_INDEPENDENT_CODE ON LINKER_LANGUAGE CXX)\n'.format(lib)
        #t += '\ntarget_compile_options({} PUBLIC -fPIC)\n'.format(lib)  # Enable Position Independent Code generation for libraries
        update_source_file(prefix + lib + '.cmake', t)

        libs.append(lib)

    return libs



def generate_cmake_file(rosetta_source_path, extra_sources):
    ''' generate cmake file in bindings source location '''

    prefix = get_binding_build_root(rosetta_source_path, source=True) + '/'

    libs = generate_rosetta_cmake_files(rosetta_source_path, prefix) + generate_rosetta_external_cmake_files(rosetta_source_path, prefix)

    if Options.external_link:
        rosetta_cmake = """
            include_directories(SYSTEM {system_include})
            include_directories({rosetta_include})
            add_definitions({defs})
            link_directories(lib)

            set(PYROSETTA_EXTERNAL_LINK ON)
            """.format(
            system_include = ' '.join(get_rosetta_system_include_directories()),
            rosetta_include = ' '.join( get_rosetta_include_directories() + [Options.pybind11] ),
            defs = ' '.join([ '-D'+d for d in get_defines()])
        )
        cmake = open('cmake.template').read()

        cmake = cmake.replace('#%__Rosetta_cmake_instructions__%#', rosetta_cmake)
        cmake = cmake.replace( '#%__PyRosetta_sources__%#', '\n'.join(extra_sources))
        cmake = cmake.replace('#%__Rosetta_libraries__%#', ' '.join(libs + ["${LINK_EXTERNAL_LIBS}"]))

        update_source_file(prefix + 'CMakeLists.txt', cmake)

    else:
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
    setup_source_directory_links(rosetta_source_path)
    maybe_version = '--version {} '.format(Options.version) if Options.version else ''
    execute('Updating version, options and residue-type-enum files...', 'cd {rosetta_source_path} && ./version.py {maybe_version}&& ./update_options.sh && ./update_ResidueType_enum_files.sh'.format(**vars()) )

    prefix = get_binding_build_root(rosetta_source_path, source=True) + '/'

    # generate include file that contains all headers
    skip_extensions = (".fwd.hh", ".impl.hh", ".py.hh")

    all_includes, serialization_instantiation = [], []
    for path in 'ObjexxFCL utility numeric basic core protocols'.split():
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

    config = open(Options.binder_config).read()
    if 'clang' not in Options.compiler: config += open('rosetta.gcc.config').read()
    with open(prefix + 'rosetta.config', 'w') as f: f.write(config)

    includes = ''.join( [' -isystem '+i for i in get_rosetta_system_include_directories()] ) + ''.join( [' -I'+i for i in get_rosetta_include_directories()] )
    defines  = ''.join( [' -D'+d for d in get_defines()] )

    if Platform == 'macos': includes += ' -isystem /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/../include/c++/v1'

    if Options.binder_extra_options: includes += ' ' + Options.binder_extra_options

    execute(
    'Generating bindings...',
    'cd {prefix} && {} --config {config} --root-module rosetta --prefix {prefix}{annotate}{trace} {} -- -std=c++11 {} {}'.format(
        Options.binder, include, includes, defines,
        prefix=prefix,
        config='./rosetta.config',
        annotate=' --annotate-includes' if Options.annotate_includes else '',
        trace=' --trace' if Options.trace else '',
    ))

    sources = open(prefix+'rosetta.sources').read().split()

    generate_cmake_file(rosetta_source_path, sources)

def build_generated_bindings(rosetta_source_path):
    ''' Build generate bindings '''
    prefix = get_binding_build_root(rosetta_source_path, build=True) + '/'

    config = '-DCMAKE_BUILD_TYPE={}'.format(Options.type)

    if Platform == "linux": config += ' -DCMAKE_C_COMPILER=`which clang` -DCMAKE_CXX_COMPILER=`which clang++`'if Options.compiler == 'clang' else ' -DCMAKE_C_COMPILER=`which gcc` -DCMAKE_CXX_COMPILER=`which g++`'

    # if we ever go back to use static libs for intermediates: fix this on Mac: -DCMAKE_RANLIB=/usr/bin/ranlib -DCMAKE_AR=/usr/bin/ar

    execute('Running CMake...', 'cd {prefix} && cmake -G Ninja {} -DPYROSETTA_PYTHON_VERSION={python_version} {py_lib} {py_include} {gcc_install_prefix} ../source'.format(config, prefix=prefix, python_version=_python_version_,
                                                                                                                                                                           py_lib='-DPYTHON_LIBRARY='+Options.python_lib if Options.python_lib else '',
                                                                                                                                                                           py_include='-DPYTHON_INCLUDE_DIR='+Options.python_include_dir if Options.python_include_dir else '',
                                                                                                                                                                           gcc_install_prefix='-DGCC_INSTALL_PREFIX='+Options.gcc_install_prefix if Options.gcc_install_prefix else ''))

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

    distutils.dir_util.copy_tree(rosetta_source_path + '/src/python/PyRosetta/documentation', documentation_prefix, update=False)

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

    execute('Building documentation...', 'cd {documentation_prefix} && PATH={python_path}:$PATH && PYTHONPATH={build_prefix}:PYTHONPATH && make clean && make html'.format(python_path=os.path.dirname(python), **vars()))

    distutils.dir_util.copy_tree(documentation_prefix + '/build/html', path, update=False)

    print('Creating PyRosetta documentation at: {}... Done!'.format(path))


def create_package(rosetta_source_path, path):
    print('Creating Python package at: {}...'.format(path))

    package_prefix = path + '/setup'
    if not os.path.isdir(package_prefix): os.makedirs(package_prefix)

    for f in 'self-test.py PyMOL-RosettaServer.py'.split(): shutil.copy(rosetta_source_path + '/src/python/PyRosetta/src/' + f, path)

    for d in 'demo test'.split(): distutils.dir_util.copy_tree(rosetta_source_path + '/src/python/PyRosetta/src/' + d, path + '/' + d, update=False)

    distutils.dir_util.copy_tree(rosetta_source_path + '/scripts/PyRosetta/public', path + '/apps', update=False)

    build_prefix = get_binding_build_root(rosetta_source_path, build=True)

    for f in 'setup.py setup.cfg ez_setup.py'.split(): shutil.copy(build_prefix + '/' + f, package_prefix)

    for d in ['pyrosetta', 'rosetta']:
        distutils.dir_util.copy_tree(build_prefix + '/' + d, package_prefix + '/' + d, update=False)

        #symlink = path + '/' + d
        #if os.path.islink(symlink): os.unlink(symlink)
        os.symlink('./setup/' + d, path + '/' + d)

        for dir_name, dirs, files in os.walk(package_prefix + '/' + d):
            if '__pycache__' in dirs: shutil.rmtree(dir_name + '/__pycache__')
            for f in files:
                if f.endswith('.pyc'): os.remove(dir_name + '/' + f)

    generate_version_file(rosetta_source_path, path + '/version.json')


def main(args):
    ''' PyRosetta building script '''

    parser = argparse.ArgumentParser()
    parser.add_argument('-j', '--jobs', default=1, const=0, nargs="?", type=int, help="Number of processors to use on when building, use '-j' with no arguments to launch job-per-core. (default: 1) ")
    parser.add_argument('-s', '--skip-generation-phase', action="store_true", help="Assume that bindings code is already generaded and skipp the Binder call's")
    parser.add_argument('-d', '--skip-building-phase', action="store_true", help="Assume that bindings code is already generaded and skipp the Binder call's")
    parser.add_argument("--type", default='Release', choices=['Release', 'Debug', 'MinSizeRel', 'RelWithDebInfo'], help="Specify build type")
    parser.add_argument('--compiler', default='clang', help='Compiler to use, defualt is clang')
    parser.add_argument('--binder', default='', help='Path to Binder tool. If none is given then download, build and install binder into main/source/build/prefix. Use "--binder-debug" to control which mode of binder (debug/release) is used.')
    parser.add_argument("--binder-debug", action="store_true", help="Run binder tool in debug mode (only relevant if no '--binder' option was specified)")
    parser.add_argument("--binder-config", default="rosetta.config", help="Binder config file. [Default='rosetta.config']")
    parser.add_argument("--print-build-root", action="store_true", help="Print path to where PyRosetta binaries will be located with given options and exit. Use this option to automate package creation.")
    parser.add_argument('--cross-compile', action="store_true", help='Specify for cross-compile build')
    parser.add_argument('--pybind11', default='', help='Path to pybind11 source tree')
    parser.add_argument('--annotate-includes', action="store_true", help='Annotate includes in generated PyRosetta source files')
    parser.add_argument('--trace', action="store_true", help='Binder will add trace output to to generated PyRosetta source files')

    parser.add_argument('-p', '--create-package', default='', help='Create PyRosetta Python package at specified path (default is to skip creating package)')
    parser.add_argument('--external-link', default=None, choices=["debug", "release"], help="Optional, link externally compiled rosetta libraries from the given cmake build directory rather than rebuilding in extension modoule.")

    parser.add_argument('--python-include-dir', default=None, help='Path to python C headers. Use this if CMake fails to autodetect it')
    parser.add_argument('--python-lib', default=None, help='Path to python library. Use this if CMake fails to autodetect it')

    parser.add_argument('--gcc-install-prefix', default=None, help='Path to GCC install prefix which will be used to determent location of libstdc++ for Binder build. Default is: auto-detected. Use this option if you would like to build Binder with compiler that was side-installed and which LLVM build system failed to identify. To see what path Binder uses for libstdc++ run `binder -- -xc++ -E -v`.')

    parser.add_argument('--serialization', action="store_true", help="Build PyRosetta with serialization enabled (off by default)")

    parser.add_argument('--binder-extra-options', default=None, help='Specify Binder extra (LLVM) options. Use this to point Binder to additional include dir or add/change LLVM flags. For example on Mac with no Xcode install set this to "-isystem /Library/Developer/CommandLineTools/usr/include/c++/v1" to point Binder to correct location of system includes.')

    parser.add_argument('--documentation', help='Generate PyRosetta documentation at specified path (default is to skip documentation creation). Note that Sphinx package need to be installed for this to work correctly.')

    parser.add_argument('--version', help='Supply JSON version file to be used for during package creation and documentation building. File must be in the same format as standard Rosetta .release.json used to mark release versions. If no file is supplied script will fallback to use main/.version.json.')

    global Options
    Options = parser.parse_args()

    rosetta_source_path = os.path.abspath('./../../../')

    binding_build_root = get_binding_build_root(rosetta_source_path)

    if Options.print_build_root: print(binding_build_root, end=('\n' if sys.stdout.isatty() else '') ); sys.exit(0)

    print('Creating PyRosetta in "{}" mode in: {}'.format(Options.type, binding_build_root))

    copy_supplemental_files(rosetta_source_path)

    if Options.skip_generation_phase:
        print('Option --skip-generation-phase is supplied, skipping generation phase...')

    else:
        if not Options.pybind11: Options.pybind11 = install_pybind11(rosetta_source_path + '/build/prefix')
        if not Options.binder:
            execute('Updating Binder and other Git submodules...', 'cd {}/.. && git submodule update --init --recursive -- source/src/python/PyRosetta/binder'.format(rosetta_source_path) )
            output = execute('Checking if Binder submodule present...',  'cd {}/.. && git submodule status'.format(rosetta_source_path), return_='output', silent=True)
            if 'source/src/python/PyRosetta/binder' not in output: print('ERROR: Binder submodule is not found... terminating...'); sys.exit(1)
            Options.binder = install_llvm_tool('binder', rosetta_source_path+'/src/python/PyRosetta/binder/source', rosetta_source_path + '/build/prefix', Options.binder_debug)

        generate_bindings(rosetta_source_path)


    if Options.skip_building_phase: print('Option --skip-building-phase is supplied, skipping building phase...')
    else: build_generated_bindings(rosetta_source_path)

    if Options.documentation: generate_documentation(rosetta_source_path, Options.documentation, Options.version)
    if Options.create_package: create_package(rosetta_source_path, Options.create_package)


if __name__ == "__main__":

    # Check if current working dir is where this file is located...
    if os.path.dirname(os.path.realpath(__file__)) != os.getcwd():
        print('This script must be run from within setup directory! Exiting...')
        sys.exit(1)

    main(sys.argv)
