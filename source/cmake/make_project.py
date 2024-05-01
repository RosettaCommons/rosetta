#!/usr/bin/env python3

# make_project.py - modified from mini-interactive by tex.

from __future__ import print_function

import os, sys, os.path, subprocess
import timeit
import build_util

# Note that CMake encodes absolute path information in the compile,
# so we don't actually lose anything by using the absolute path
PATH_TO_SOURCE_DIR = os.path.dirname( os.path.dirname( os.path.abspath(sys.argv[0]) ) ) + "/"
PATH_TO_SOURCE_DIR = PATH_TO_SOURCE_DIR.replace('\\', '/')


def execute_rosetta_script_from_source_dir(script_name, parser='bash'):
    if parser == 'python':
        parser = sys.executable
    command = parser + ' ' + script_name
    subprocess.run(command, cwd="..", check=True, stderr=subprocess.STDOUT, shell=True)


def update_version(): execute_rosetta_script_from_source_dir('version.py','python')

def update_options(): execute_rosetta_script_from_source_dir('update_options.sh','bash')

def update_submodules(): execute_rosetta_script_from_source_dir('update_submodules.sh','bash')

def update_ResidueType_enum_files(): execute_rosetta_script_from_source_dir('update_ResidueType_enum_files.sh','bash')


def project_callback(project, project_path, project_files):
    print('making project files for project ' + project + ' ...', end='')

    project_path = project_path.replace("\\",'/')
    platform = sys.platform

    if platform == 'darwin':
        platform = 'macos'
    elif platform.startswith('linux'):
        platform = 'linux'

    if project == 'apps' or project == 'pilot_apps':
        app_files = {}
        #print('making project files for ' + project + ' from ' + project_path)
        apps_output = open( 'build/' + project + '.all.cmake', 'w' )

        if not os.path.exists( 'build/apps' ):
            os.mkdir( 'build/apps' )
        for entry in project_files:
            file_path = '/'.join( [ '../../src/apps', entry[0] ] )
            for file in entry[ 1 ]:
                full_path = file_path + '/' + file
                split_file = file.split( '.' )
                extension = split_file.pop()

                if extension == 'cc':
                    tag = ''.join( split_file )
                    app_files[ os.path.basename( tag ) ] = full_path
        for key in app_files:
            output = ''
            output += 'ADD_EXECUTABLE( %s %s )' % (key, app_files[ key ] ) + '\n'
            output += 'TARGET_LINK_LIBRARIES( %s  ${LINK_ALL_LIBS} )' % ( key ) + '\n'
            output += 'SET_TARGET_PROPERTIES( %s PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}" )' % ( key ) + '\n'
            output += 'SET_TARGET_PROPERTIES( %s PROPERTIES LINK_FLAGS "${LINK_FLAGS}" )' % ( key ) + '\n'

            symlink_var = key + '_symlink'
            output += 'ADD_CUSTOM_TARGET( %s  ALL)\n' % ( symlink_var )
            output += "# cmake -E create_symlink won't choke if the symlink already exists\n"
            output += "ADD_CUSTOM_COMMAND( TARGET %s POST_BUILD\n" % ( symlink_var )
            output += "\t\tCOMMAND ${CMAKE_COMMAND} -E make_directory ../../bin/ \n"
            output += "\t\tCOMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_BINARY_DIR}/%s" % ( key ) # no endline
            output += " ../../bin/%s\n" % ( key )
            output += "\t\tCOMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_BINARY_DIR}/%s" % ( key ) # no endline
            output += " ../../bin/%s.%s${COMPILER}${MODE}\n" % ( key, platform )
            output += "\t\tCOMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_BINARY_DIR}/%s" % ( key ) # no endline
            output += " ../../bin/%s.default.%s${COMPILER}${MODE}\n" % ( key, platform )
            output += ")\n"
            output += 'ADD_DEPENDENCIES( %s %s )\n' % ( symlink_var, key )
            output += 'ADD_DEPENDENCIES( %s BUILD_ROSETTA_LIBS )\n' % ( key )
            output += 'ADD_DEPENDENCIES( %s %s )\n' % ( project, symlink_var )
            output += 'install(TARGETS %(key)s RUNTIME DESTINATION bin OPTIONAL)' % dict(key=key)
            apps_file = key + '.cmake'
            open( 'build/apps/' + apps_file, 'w').write( output )
            apps_output.write( 'INCLUDE( ../build/apps/%s )\n' % apps_file )
        apps_output.close()
    else:
        cmake_files = ''
        for dir, files in project_files:
            cmake_files += '\n\t' + '\n\t'.join([project_path + dir + file for file in files])

        #print('making project files for ' + project + ' from ' + project_path)

        output = 'SET(' + project + '_files' + cmake_files + '\n)\n'

        open('build/' + project + '.cmake', 'w').write(output)
    print('done.')

def project_external_callback(project, project_path, project_files, other_settings):
    print('making external project files for project ' + project + ' ...', end='')
    project_path = project_path.replace("\\",'/')

    cmake_files = ''
    for dir, files in project_files:
        cmake_files += '\n\t' + '\n\t'.join([project_path + dir + file for file in files])

    cmake_defines = ''
    cmake_compileflags = ''
    cmake_linkflags = ''

    if 'defines' in other_settings:
        cmake_defines = ';'.join( other_settings['defines'] )

    compile_flags = []
    if 'ccflags' in other_settings:
        compile_flags.extend( other_settings['ccflags'] )
    #Unfortunately, CMAKE doesn't have language-specific flags settings
    if 'cflags' in other_settings:
        compile_flags.extend( other_settings['cflags'] )
    if 'cxxflags' in other_settings:
        compile_flags.extend( other_settings['cxxflags'] )
    if compile_flags:
        cmake_compileflags = ' '.join( compile_flags )

    if 'link_flags' in other_settings:
        cmake_linkflags = ';'.join( other_settings['link_flags'] )

    output = ''
    output += 'SET(' + project + '_files' + cmake_files + '\n)\n'
    output += 'SET(' + project + '_defines ' + cmake_defines + ')\n'
    output += 'SET(' + project + '_compileflags ' + cmake_compileflags + ')\n'
    output += 'SET(' + project + '_linkflags ' + cmake_linkflags + ')\n'

    open('build/external_' + project + '.cmake', 'w').write(output)
    print('done.')

    return True

def project_test_callback(test, project_path, test_path, test_files, test_inputs):
    print('making test files for test ' + test + ' ...', end='')
    project_path = project_path.replace("\\",'/')

    output = ''
    cmake_headers = []
    cmake_testfiles = []
    cmake_testinputs = []
    cmake_directories = set()

    for dir, files in test_files:
            cmake_directories.add(dir)
            cmake_directories.update( [dir + cxx.rsplit('/',1)[0] for (header, cxx, root) in files if cxx is not None and '/' in cxx] )
            cmake_testfiles.extend( [dir + header[:-3] + ".cc" for (header, cxx, root) in files if root != True and cxx != None] )
            cmake_headers.extend( [dir + header for (header, cxx, root) in files if root != True and cxx == None] )

    for dir, files in test_inputs:
            cmake_directories.add(dir)
            cmake_directories.update( [dir + file.rsplit('/',1)[0] for file in files if '/' in file] )
            cmake_testinputs.extend( [dir + file for file in files] )


    #Make directories alphabetical, so parent directories are created before child directories
    cmake_directories = list(cmake_directories)
    cmake_directories.sort()

    output += 'SET(' + test + '_testfiles\n\t' + '\n\t'.join(cmake_testfiles) + '\n)\n\n'
    output += 'SET(' + test + '_testdirectories\n\t' + '\n\t'.join(cmake_directories) + '\n)\n\n'
    output += 'SET(' + test + '_testheaders\n\t' + '\n\t'.join(cmake_headers) + '\n)\n\n'
    output += 'SET(' + test + '_testinputs\n\t' + '\n\t'.join(cmake_testinputs) + '\n)\n\n'

    open(test_path + "test_" + test + '.cmake', 'w').write(output)
    print('done.')

# get in the script directory
rosetta_cmake_directory = os.path.dirname( sys.argv[0] ) # where this script is located
os.chdir( rosetta_cmake_directory or "./" )

update_version()
update_options()
update_submodules()  # We do have the update in the build process, but if we don't trigger that properly we can error out before then.
update_ResidueType_enum_files()

build_util.external_main(PATH_TO_SOURCE_DIR, sys.argv, project_external_callback)

build_util.project_main(PATH_TO_SOURCE_DIR, sys.argv, project_callback)

build_util.test_main(PATH_TO_SOURCE_DIR, sys.argv, project_test_callback)
