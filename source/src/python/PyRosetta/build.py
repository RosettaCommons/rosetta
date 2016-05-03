#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true: :collapseFolds=10:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   build.py
## @brief  Build PyRosetta
## @author Sergey Lyskov


import os, sys, argparse, platform, subprocess

from collections import OrderedDict

if sys.platform.startswith("linux"): Platform = "linux" # can be linux1, linux2, etc
elif sys.platform == "darwin" : Platform = "macos"
elif sys.platform == "cygwin" : Platform = "cygwin"
elif sys.platform == "win32" : Platform = "windows"
else: Platform = "unknown"
PlatformBits = platform.architecture()[0][:2]


_banned_dirs_ = 'src/utility/pointer src/protocols/jd3'.split()  # src/utility/keys src/utility/options src/basic/options
_banned_headers_ = 'utility/py/PyHelper.hh utility/keys/KeyCount.hh utility/keys/KeyLookup.functors.hh'
_banned_headers_ +=' core/scoring/fiber_diffraction/FiberDiffractionKernelGpu.hh' # GPU code

# moved to config file _banned_namespaces_ = 'utility::options'.split()

def is_dir_banned(directory):
    for b in _banned_dirs_:
        if b in directory: return True
    return False


def get_rosetta_include_directories():
    ''' return list of include directories for compilation '''
    r = 'src external external/include external/boost_1_55_0 external/dbio external/dbio/sqlite3 external/libxml2/include'.split()
    r.append('src/platform/'+Platform)
    return r


def get_defines():
    ''' return list of #defines '''
    return 'PYROSETTA BOOST_ERROR_CODE_HEADER_ONLY BOOST_SYSTEM_NO_DEPRECATED BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS CXX11 PTR_STD'.split() + ([] if Options.debug else ['NDEBUG'])



def execute(message, command_line, return_='status', until_successes=False, terminate_on_failure=True):
    print(message);  print(command_line)
    while True:

        p = subprocess.Popen(command_line, bufsize=0, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, errors = p.communicate()

        output = output.decode('utf-8') + errors.decode('utf-8')
        exit_code = p.returncode

        print(output)

        if exit_code and until_successes: pass  # Thats right - redability COUNT!
        else: break

        print( "Error while executing {}: {}\n".format(message, output) )
        print("Sleeping 60s... then I will retry...")
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



def get_binding_build_root(rosetta_source_path, source=False, build=False):
    ''' Calculate bindings build path using current platform and compilation settings and create dir if it is not there yet '''

    p = os.path.join(rosetta_source_path, 'build/PyRosetta')

    #p = os.path.join(p, 'cross_compile' if Options.cross_compile else (Platform+ '/' + options.compiler) )
    p =  os.path.join(p, Platform+ '/' + get_compiler_family() )

    #p = os.path.join(p, 'monolith' if True else 'namespace' )

    p = os.path.join(p, 'debug' if Options.debug  else 'release')
    #p = os.path.abspath( os.path.join(p, '_build_' if Options.monolith else 'rosetta') )

    source_p = p+'/source'
    build_p  = p+'/build'

    if not os.path.isdir(source_p): os.makedirs(source_p)
    if not os.path.isdir(build_p) : os.makedirs(build_p)

    if source: return source_p
    if build:  return build_p

    return p



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
        t  = 'add_library({}\n{})\n'.format(l, '\n'.join( [ rosetta_source_path + '/external/' + s for s in libs[l]] ))
        t += '\ntarget_compile_options({} PUBLIC -fPIC {})\n'.format(l, ' '.join([ '-D'+d for d in defines[l].split() ] ) )   #  target_compile_definitions
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

        t  = 'add_library({}\n{})\n'.format(lib, '\n'.join( [ rosetta_source_path + '/src/' + s for s in sources] ))
        t += '\ntarget_compile_options({} PUBLIC -fPIC)\n'.format(lib)  # Enable Position Independent Code generation for libraries
        update_source_file(prefix + lib + '.cmake', t)

        libs.append(lib)

    return libs



def generate_cmake_file(rosetta_source_path, extra_sources):
    ''' generate cmake file in bindings source location '''

    prefix = get_binding_build_root(rosetta_source_path, source=True) + '/'

    libs = generate_rosetta_cmake_files(rosetta_source_path, prefix) + generate_rosetta_external_cmake_files(rosetta_source_path, prefix)

    rosetta_cmake =  ''.join( ['include({}.cmake)\n'.format(l) for l in libs] )
    rosetta_cmake += '\ninclude_directories({})\n\n'.format(' '.join(get_rosetta_include_directories()+[Options.pybind11] ) )
    rosetta_cmake += 'add_definitions({})\n'.format(' '.join([ '-D'+d for d in get_defines()] ) )

    cmake = open('template.cmake').read()

    cmake = cmake.replace('#%__Rosetta_cmake_instructions__%#', rosetta_cmake)
    cmake = cmake.replace('#%__PyRosetta_sources__%#', '\n'.join(extra_sources))
    cmake = cmake.replace('#%__Rosetta_libraries__%#', ' '.join(libs))

    update_source_file(prefix + 'CMakeLists.txt', cmake)

    #test = prefix + 'test.cpp'
    #with open(test, 'w') as f: f.write('#include <utility/exit.hh>\nint main(void) { utility_exit_with_message("Wow..."); return 0; }\n')
    # with open(prefix + 'CMakeLists.txt', 'w') as f:
    #     f.write('cmake_minimum_required(VERSION 2.8)\n\nproject(rosetta)\n\n')

    #     libs = generate_rosetta_cmake_files(rosetta_source_path, prefix)

    #     for i in libs: f.write('include({}.cmake)\n'.format(i))

    #     for d in get_rosetta_include_directories(): f.write( 'include_directories({})\n'.format(d) )

    #     for d in get_defines(): f.write( 'add_definitions(-D{})\n'.format(d) )

    #     f.write('set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")\n')

    #     #f.write( 'add_executable(test {})\n'.format(test) )
    #     #f.write( 'target_link_libraries(test {})\n'.format(' '.join(libs) ) )

    #     f.write( 'target_link_libraries(test {})\n'.format(' '.join(libs) ) )
    #     #extra_sources


    # extra_objs = [  # additioanal source for external libs
    #     'dbio/cppdb/atomic_counter.cpp', "dbio/cppdb/conn_manager.cpp", "dbio/cppdb/driver_manager.cpp", "dbio/cppdb/frontend.cpp",
    #     "dbio/cppdb/backend.cpp", "dbio/cppdb/mutex.cpp", "dbio/cppdb/pool.cpp", "dbio/cppdb/shared_object.cpp", "dbio/cppdb/sqlite3_backend.cpp",
    #     "dbio/cppdb/utils.cpp", 'dbio/sqlite3/sqlite3.c', ]


def generate_bindings(rosetta_source_path):
    ''' Generate bindings using binder tools and return list of source files '''
    prefix = get_binding_build_root(rosetta_source_path, source=True) + '/'

    for d in ['src', 'external']:
        s = prefix + d
        if os.path.islink(s): os.unlink(s)
        os.symlink(rosetta_source_path + '/' + d, s)

    # generate include file that contains all headers
    include = prefix + 'all_rosetta_includes.hh'
    with open(include, 'w') as fh:
        #for path in 'ObjexxFCL utility numeric basic core protocols'.split():
        for path in 'ObjexxFCL utility numeric basic core protocols'.split():
            for dir_name, _, files in os.walk(rosetta_source_path + '/src/' + path):
                for f in files:
                    if f.endswith('.hh'):
                        #if f == 'exit.hh':
                        #if dir_name.endswith('utility'):
                        if not is_dir_banned(dir_name):
                        #if 'utility/' in dir_name  and  'src/utility/pointer' not in dir_name:
                            header = dir_name[len(rosetta_source_path+'/src/'):] + '/' + f
                            if header not in _banned_headers_  and  not header.startswith('basic/options/keys/OptionKeys.cc.gen'):
                                #print(header)
                                fh.write( '#include <{}>\n'.format(header) )


    includes = ''.join( [' -I'+i for i in get_rosetta_include_directories()] )
    defines  = ''.join( [' -D'+d for d in get_defines()] )

    execute('Generating bindings...', 'cd {prefix} && {} --config {config} --annotate-includes --root-module rosetta --prefix {prefix} {} -- -std=c++11 {} {}'.format(Options.binder, include, includes, defines, prefix=prefix, config=os.path.abspath('./rosetta.config') ) )

    sources = open(prefix+'rosetta.sources').read().split()

    generate_cmake_file(rosetta_source_path, sources)


def  build_generated_bindings(rosetta_source_path):
    ''' Build generate bindings '''
    prefix = get_binding_build_root(rosetta_source_path, build=True) + '/'

    config = '-DCMAKE_BUILD_TYPE={}'.format('Debug' if Options.debug else 'Release')
    config += ' -DCMAKE_C_COMPILER=`which clang` -DCMAKE_CXX_COMPILER=`which clang++`'if Options.compiler == 'clang' else ' -DCMAKE_C_COMPILER=`which gcc` -DCMAKE_CXX_COMPILER=`which g++`'

    execute('Running CMake...', 'cd {prefix} && cmake -G Ninja {} ../source'.format(config, prefix=prefix))

    execute('Building...', 'cd {prefix} && ninja'.format(prefix=prefix))



def main(args):
    ''' PyRosetta building script '''

    parser = argparse.ArgumentParser()

    parser.add_argument("--debug", action="store_true", help="Build bindings in debug mode")

    parser.add_argument('--compiler', default='gcc', help='Compiler to use, defualt is gcc on Linux and clang on Mac')

    parser.add_argument('--binder', default='binder', help='Path to Binder tool')

    parser.add_argument("--print-build-root", action="store_true", help="Print path to where PyRosetta binaries will be located with given options and exit. Use this option to automate package creation.")

    parser.add_argument('--cross-compile', action="store_true", help='Specify for cross-compile build')

    parser.add_argument('--pybind11', default='', help='Path to pybind11 source tree')

    global Options
    Options = parser.parse_args()

    rosetta_source_path = os.path.abspath('./../../../')

    binding_build_root = get_binding_build_root(rosetta_source_path)

    if Options.print_build_root: print(binding_build_root, end=''); sys.exit(0)

    generate_bindings(rosetta_source_path)

    build_generated_bindings(rosetta_source_path)


if __name__ == "__main__":
    main(sys.argv)
