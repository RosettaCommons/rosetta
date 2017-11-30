#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   update_ui_project.py
## @brief  create/update ui project files
## @author Sergey Lyskov

from __future__ import print_function

import os, sys, argparse, platform, imp, subprocess

from collections import OrderedDict


if sys.platform.startswith("linux"): Platform = "linux" # can be linux1, linux2, etc
elif sys.platform == "darwin" : Platform = "macos"
elif sys.platform == "cygwin" : Platform = "cygwin"
elif sys.platform == "win32" : Platform = "windows"
else: Platform = "unknown"
PlatformBits = platform.architecture()[0][:2]


def get_rosetta_external_libraries(): return 'z'.split()


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
    defines = 'BOOST_ERROR_CODE_HEADER_ONLY BOOST_SYSTEM_NO_DEPRECATED BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS PTR_STD'  # MULTI_THREADED
    if Platform == 'macos': defines += ' UNUSUAL_ALLOCATOR_DECLARATION'
    #if Options.type in 'Release MinSizeRel': defines += ' NDEBUG'
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




library_project_template = '''\
QT -= core gui

CONFIG += object_parallel_to_source no_keywords {config}

TARGET = {name}
TEMPLATE = lib

DEFINES += {defines}

INCLUDEPATH = {includes}

QMAKE_CXXFLAGS += {flags}

QMAKE_CFLAGS += {flags}

SOURCES += \\
\t{sources}

HEADERS += \\
\t{headers}

LIBS +={libs}
'''


def generate_rosetta_external_project_files(rosetta_source_path, prefix):
    # libs = OrderedDict([
    #     ('sqlite3_r', ['dbio/sqlite3/sqlite3.c', 'dbio/sqlite3/sqlite3.h' ] ),
    #     ('cppdb', [ 'dbio/cppdb/' + f for f in 'atomic_counter.cpp atomic_counter.h backend.cpp backend.h connection_specific.h conn_manager.cpp conn_manager.h defs.h driver_manager.cpp driver_manager.h errors.h frontend.cpp frontend.h mutex.cpp mutex.h numeric_util.h pool.cpp pool.h ref_ptr.h shared_object.cpp shared_object.h sqlite3_backend.cpp utils.cpp utils.h'.split() ] ),
    # ])
    # defines = dict(cppdb = r'CPPDB_EXPORTS CPPDB_DISABLE_SHARED_OBJECT_LOADING CPPDB_DISABLE_THREAD_SAFETY CPPDB_WITH_SQLITE3' +
    #                        r' CPPDB_LIBRARY_PREFIX=\\\"lib\\\" CPPDB_LIBRARY_SUFFIX=\\\".dylib\\\" CPPDB_SOVERSION=\\\"0\\\" ' +
    #                        r' CPPDB_MAJOR=0 CPPDB_MINOR=3 CPPDB_PATCH=0 CPPDB_VERSION=\\\"0.3.0\\\"',
    #                sqlite3_r = 'SQLITE_DISABLE_LFS SQLITE_OMIT_LOAD_EXTENSION SQLITE_THREADSAFE=0')

    libs = OrderedDict([
        ('external', ['dbio/sqlite3/sqlite3.c', 'dbio/sqlite3/sqlite3.h' ] + \
                     [ 'dbio/cppdb/' + f for f in 'atomic_counter.cpp atomic_counter.h backend.cpp backend.h connection_specific.h conn_manager.cpp conn_manager.h defs.h driver_manager.cpp driver_manager.h errors.h frontend.cpp frontend.h mutex.cpp mutex.h numeric_util.h pool.cpp pool.h ref_ptr.h shared_object.cpp shared_object.h sqlite3_backend.cpp utils.cpp utils.h'.split() ] ),
    ])

    defines = dict(external = r'CPPDB_EXPORTS CPPDB_DISABLE_SHARED_OBJECT_LOADING CPPDB_DISABLE_THREAD_SAFETY CPPDB_WITH_SQLITE3' +
                              r' CPPDB_LIBRARY_PREFIX=\\\"lib\\\" CPPDB_LIBRARY_SUFFIX=\\\".dylib\\\" CPPDB_SOVERSION=\\\"0\\\" ' +
                              r' CPPDB_MAJOR=0 CPPDB_MINOR=3 CPPDB_PATCH=0 CPPDB_VERSION=\\\"0.3.0\\\"' +
                              r' SQLITE_DISABLE_LFS SQLITE_OMIT_LOAD_EXTENSION SQLITE_THREADSAFE=0')

    scons_file_extension = '.external.settings'
    external_scons_files = [f for f in os.listdir(rosetta_source_path+'/external') if f.endswith(scons_file_extension)]
    for scons_file in external_scons_files:
        G = {}
        sources= []
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

    linking_libs = ''
    for lib in libs:
        lib_prefix = prefix + '/' + lib + '/'
        if not os.path.isdir(lib_prefix): os.makedirs(lib_prefix)

        t = library_project_template.format(name=lib,
                                            config='',
                                            sources = ' \\\n\t'.join( [ 'external/' + s for s in libs[lib] if not s.endswith('h') ] ),
                                            headers = ' \\\n\t'.join( [ 'external/' + h for h in libs[lib] if h.endswith('h') ] ),
                                            #sources = ' \\\n\t'.join( [ os.path.abspath(rosetta_source_path + '/external/' + s) for  if not s.endswith('h') ] ),
                                            #headers = ' \\\n\t'.join( [ os.path.abspath(rosetta_source_path + '/external/' + h) for h in libs[lib] if h.endswith('h') ] ),
                                            defines = defines[lib] + ' '+ ' '.join( get_defines() ),
                                            includes = ' '.join( [ '$$PWD/../../../../' + i for i in get_rosetta_include_directories() ] ),
                                            flags = ' '.join( [ '-isystem$$PWD/../../../../' + i for i in get_rosetta_system_include_directories() ] ),
                                            libs = linking_libs,
        )

        update_source_file(lib_prefix + lib + '.pro', t)

        # for s in libs[lib]:
        #     filename = 'external/' + s
        #     dirname = os.path.dirname(filename)
        #     if not os.path.isdir(lib_prefix+dirname): os.makedirs(lib_prefix+dirname)
        #     if not os.path.islink(lib_prefix+filename): os.link( os.path.abspath(rosetta_source_path + '/' + filename), lib_prefix + filename)


        for l in 'external'.split():
            s = lib_prefix + l
            if not os.path.islink(s): os.symlink('../../../../'+l, s)


        linking_libs += ' \\\n\t-L$$OUT_PWD/../{0}/ -l{0}'.format(lib)


    return list( libs.keys() ), linking_libs


def generate_rosetta_libraries_project_files(rosetta_source_path, prefix, linking_libs):
    scons_file_prefix = rosetta_source_path + '/src/'
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

    all_libs.sort(key=key, reverse=False)

    for lib in all_libs:
        lib_mangled = lib.replace('.', '_')
        #if lib not in 'ObjexxFCL utility numeric basic': continue
        #if lib not in 'ObjexxFCL utility numeric basic core.1 core.2 core.3 core.4': continue

        G = {}
        sources, headers = [], []

        with open(scons_file_prefix + lib + lib_suffix) as f:
               code = compile(f.read(), scons_file_prefix + lib + lib_suffix, 'exec')
               exec(code, G)

        for dir_ in G['sources']:
            for f in G['sources'][dir_]:
                if not f.endswith('.cu'): sources.append( dir_ + '/' + f + '.cc')
                #if f.endswith('.cu'): sources.append( dir_ + '/' + f )
                #else: sources.append( dir_ + '/' + f + '.cc')

            for h in os.listdir(scons_file_prefix + dir_):
                if h.endswith('.hh') or h.endswith('.h'): headers.append(dir_ + '/' + h)
                # # only add headers if .cc could be found...
                # if h.endswith('.hh') or h.endswith('.h'):
                #     base = h.replace('.fwd.hh', '').replace('.hh', '').replace('.h', '')
                #     if os.path.isfile(scons_file_prefix + dir_ + '/' + base + '.cc'): headers.append(dir_ + '/' + h)


        sources.sort()

        lib_prefix = prefix + '/' + lib_mangled + '/'
        if not os.path.isdir(lib_prefix): os.makedirs(lib_prefix)

        t = library_project_template.format(name=lib_mangled,
                                            config='c++11',
                                            sources = ' \\\n\t'.join([s for s in sources]),
                                            headers = ' \\\n\t'.join([h for h in headers]),
                                            defines = ' '.join( get_defines() ),
                                            includes = ' '.join( [ '$$PWD/../../../../' + i for i in get_rosetta_include_directories() ] ),
                                            flags = ' '.join( [ '-isystem$$PWD/../../../../' + i for i in get_rosetta_system_include_directories() ] ),
                                            libs = linking_libs,
        )

        update_source_file(lib_prefix + lib_mangled + '.pro', t)

        # for filename in sources + headers:
        #     dirname = os.path.dirname(filename)
        #     if not os.path.isdir(lib_prefix+dirname): os.makedirs(lib_prefix+dirname)
        #     if not os.path.islink(lib_prefix+filename): os.link( os.path.abspath(rosetta_source_path + '/src/' + filename), lib_prefix + filename)

        for l in 'ObjexxFCL utility numeric basic core protocols'.split():
            s = lib_prefix + l
            if not os.path.islink(s): os.symlink('../../../../src/'+l, s)

        # s = lib_prefix + 'platform'
        # if not os.path.islink(s): os.symlink('../../../../src/platform/' + Platform + '/platform', s)

        libs.append(lib_mangled)
        linking_libs += ' \\\n\t-L$$OUT_PWD/../{0}/ -l{0}'.format(lib_mangled)

    return libs



subdir_project_template = '''\
TEMPLATE = subdirs

SUBDIRS += \\
\t{subdirs}
'''

def generate_rosetta_project_files(rosetta_source_path, prefix):
    external_libs, linking_string = generate_rosetta_external_project_files(rosetta_source_path, prefix)
    rosetta_libs = generate_rosetta_libraries_project_files(rosetta_source_path, prefix, linking_string + ''.join( [' -l' + l for l in get_rosetta_external_libraries()] ) )

    libs = external_libs + rosetta_libs

    project = subdir_project_template.format(subdirs = ' \\\n\t'.join(libs) )

    for i in range(1, len(libs)): project += '\n{}.depends\t{}= {}\n'.format(libs[i], '\t' if len( libs[i]) < 8  else '', ' '.join(libs[:i]) )

    update_source_file(prefix + '/rosetta.pro',  project)

    return libs



def generate_app_project_files(rosetta_source_path, prefix):
    ''' create link to ui project files and return list of projects
    '''
    apps_project_root = prefix

    #if not os.path.isdir(apps_project_root): os.makedirs(apps_project_root)
    #projects = [ p for p in os.listdir(rosetta_source_path+'/src/ui/apps/public') if os.path.isdir(rosetta_source_path+'/src/ui/apps/public/'+p) ]

    projects = ['workbench', 'tests']  # WARNING: DO NOT ADD EXTRA APPS to this list! We intentionally want to keep only one UI app public to enforce merge-ability of apps. If you feel that you need to change this - please write to devel list first!


    pilot_filename = 'config.py'
    if os.path.isfile(pilot_filename):
        pilot = imp.load_source('pilot', pilot_filename)

        if hasattr(pilot, 'pilot_apps'):
            for a in pilot.pilot_apps:
                pth, name = a
                print('Adding Pilot app {}/{}.pro...'.format(pth, name) )

                if not os.path.islink(apps_project_root + '/' + name): os.symlink('../../src/ui/apps/pilot/'+pth, apps_project_root + '/' + name)
                projects.append(name)


        # if hasattr(pilot, 'public_apps'):
        #     print( 'Warning local config for public_apps detected, setting public_app={}'.format(pilot.public_apps) )
        #     projects = pilot.public_apps

    for p in projects:
        if not os.path.islink(apps_project_root + '/' + p): os.symlink('../../src/ui/apps/'+p, apps_project_root + '/' + p)




    #project = 'TEMPLATE = subdirs\n\nSUBDIRS += {}\n'.format( ' '.join(projects) )
    #for p in projects: project +='\n{}.depends = rosetta\n'.format(p)
    #update_source_file(apps_project_root + '/apps.pro', project)

    return projects


""" version for putting apps into qt/apps subproject, but qmake does not allow to specify project-level dependency for project on different hierarchy leveles
def generate_app_project_files(rosetta_source_path, prefix):
    ''' create link to ui project files and return list of projects
    '''
    apps_project_root = prefix + '/apps'

    if not os.path.isdir(apps_project_root): os.makedirs(apps_project_root)

    #projects = [ p[:-len('.pro')] for p in os.listdir(rosetta_source_path+'/src/ui') if p.endswith('.pro') ]
    projects = [ p for p in os.listdir(rosetta_source_path+'/src/ui/apps') if os.path.isdir(rosetta_source_path+'/src/ui/apps/'+p) ]

    for p in projects:
        if not os.path.islink(apps_project_root + '/' + p): os.symlink('../../../src/ui/apps/'+p, apps_project_root + '/' + p)

    project = 'TEMPLATE = subdirs\n\nSUBDIRS += {}\n'.format( ' '.join(projects) )
    # we can only set dependency for projects on the same level, so line below is commented out
    #for p in projects: project +='\n{}.depends = rosetta\n'.format(p)
    update_source_file(apps_project_root + '/apps.pro', project)

    return projects
"""


app_run_shell_script_template = '''\
#!/bin/bash

# Shell script to build/run app {name} from command line
# To use thi script copy it into build location (if you using Qt Creator default will be: main/source/src/build/build-qt-Desktop_Qt_<version>)
# Then uncomment desirable build commands below

# build all Rosetta libs. Unless you developeing Rosetta code this step will only need to be done once and could be executed from Qt Creator
# so it is commented by default
# cd rosetta && make $@ && cd ..

# build single Rosetta lib (utility in this example). Uncomment this if you working on particular Rosetta lib.
# cd rosetta/utility && make $@ && cd ../..

# build ui lib
cd ui && make $@ && cd ..

# build {name} app
cd {name} && make $@ && cd ..

export DYLD_LIBRARY_PATH={rosetta_libs}:`pwd`/ui:${{DYLD_LIBRARY_PATH+:$DYLD_LIBRARY_PATH}}

./{name}/{name}.app/Contents/MacOS/{name}
'''

def generate_app_shell_scripts(prefix, rosetta_libs, apps):
    import stat
    for app in apps:
        shell_script = prefix + '/' + app + '.sh'
        print('Creating app build/run script: {}'.format(shell_script))
        with open(shell_script, 'w') as f:
            f.write(app_run_shell_script_template.format(name=app,
                                                         rosetta_libs=':'.join( ['`pwd`/rosetta/'+l for l in rosetta_libs]))
            )

        os.chmod(shell_script, os.stat(shell_script).st_mode | stat.S_IEXEC)


def main(args):
    ''' UI project building script '''

    parser = argparse.ArgumentParser()

    global Options
    Options = parser.parse_args()

    rosetta_source_path = os.path.abspath('./../../')

    print( 'Creating/Updating Qt project files and assuming that Rosetta source dir is at {}...'.format(rosetta_source_path) )

    execute('Updating version, options and residue-type-enum files...', 'cd {} && ./version.py && ./update_options.sh && ./update_ResidueType_enum_files.sh'.format(rosetta_source_path) )

    ui_project_root = os.path.abspath(rosetta_source_path + '/build/qt')

    if not os.path.isdir(ui_project_root): os.makedirs(ui_project_root)

    rosetta_libs = generate_rosetta_project_files(rosetta_source_path, ui_project_root + '/rosetta')

    projects = generate_app_project_files(rosetta_source_path, ui_project_root)

    # ui lib
    if not os.path.islink(ui_project_root + '/ui'): os.symlink('../../src/ui', ui_project_root + '/ui')
    projects.append('ui')

    # creating root project, this is just an umbrella that hold other projects
    root_project = 'TEMPLATE = subdirs\n\nSUBDIRS += {} rosetta\n'.format( ' '.join(projects) )
    for p in projects: root_project +='\n{}.depends = rosetta {}\n'.format(p, 'ui' if p!='ui' else '')
    update_source_file(ui_project_root + '/qt.pro', root_project)

    # root_project = 'TEMPLATE = subdirs\n\nSUBDIRS += apps rosetta\napps.depends = rosetta\n'
    # update_source_file(ui_project_root + '/qt.pro', root_project)

    generate_app_shell_scripts(ui_project_root, rosetta_libs, [a for a in projects if a != 'ui'])


if __name__ == "__main__":

    # Check if current working dir is where this file is located...
    if os.path.dirname(os.path.realpath(__file__)) != os.getcwd():
        print('This script must be run from within source/src/ui directory! Exiting...')
        sys.exit(1)

    main(sys.argv)
