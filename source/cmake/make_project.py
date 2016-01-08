#!/usr/bin/env python2.7

# make_project.py - modified from mini-interactive by tex.

import os, string, sys, os.path
import timeit
import build_util

# Note that CMake encodes absolute path information in the compile,
# so we don't actually lose anything by using the absolute path
PATH_TO_SOURCE_DIR = os.path.dirname( os.path.dirname( os.path.abspath(sys.argv[0]) ) ) + "/"

def update_version():
	cmd = 'cd ..; python2.7 version.py'
	os.system(cmd)

def update_options():
	cmd = 'cd ..; ./update_options.sh'
	os.system(cmd)

def update_ResidueType_enum_files():
	cmd = 'cd ..; ./update_ResidueType_enum_files.sh'
	os.system(cmd)


def project_callback(project, project_path, project_files):
	print 'making project files for project ' + project + ' ...',

	if project == 'apps' or project == 'pilot_apps':
		app_files = {}
		#print 'making project files for ' + project + ' from ' + project_path
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
			output += 'ADD_CUSTOM_COMMAND( TARGET %s POST_BUILD COMMAND python2.7 ../smart_symlink.py ${COMPILER} ${MODE} %s )\n' % ( symlink_var, key )
			output += 'ADD_DEPENDENCIES( %s %s )\n' % ( symlink_var, key )
			output += 'ADD_DEPENDENCIES( %s BUILD_ROSETTA_LIBS )\n' % ( key )
			output += 'ADD_DEPENDENCIES( %s %s )\n' % ( project, symlink_var )
			apps_file = key + '.cmake'
			open( 'build/apps/' + apps_file, 'w').write( output )
			apps_output.write( 'INCLUDE( ../build/apps/%s )\n' % apps_file )
		apps_output.close()
	else:
		cmake_files = ''
		for dir, files in project_files:
			cmake_files += '\n\t' + string.join([project_path + dir + file for file in files], '\n\t')

		#print 'making project files for ' + project + ' from ' + project_path

		output = ''
		output += 'SET(' + project + '_files' + cmake_files + '\n)\n'

		open('build/' + project + '.cmake', 'w').write(output)
	print 'done.'

def project_external_callback(project, project_path, project_files, defines):
	print 'making external project files for project ' + project + ' ...',

	cmake_files = ''
	for dir, files in project_files:
		cmake_files += '\n\t' + string.join([project_path + dir + file for file in files], '\n\t')

        cmake_defines = ''

        if defines :
            cmake_defines = ';'.join( defines )

	#print 'making project files for ' + project + ' from ' + project_path

	output = ''
	output += 'SET(' + project + '_files' + cmake_files + '\n)\n'
        output += 'SET(' + project + '_defines ' + cmake_defines + ')\n'

	open('build/external_' + project + '.cmake', 'w').write(output)
	print 'done.'

def project_test_callback(test, project_path, test_path, test_files, test_inputs):
	print 'making test files for test ' + test + ' ...',

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
	print 'done.'

# get in the script directory
rosetta_cmake_directory = os.path.dirname( sys.argv[0] ) # where this script is located
os.chdir( rosetta_cmake_directory or "./" )

update_version()
update_options()
update_ResidueType_enum_files()

build_util.external_main(PATH_TO_SOURCE_DIR, sys.argv, project_external_callback)

build_util.project_main(PATH_TO_SOURCE_DIR, sys.argv, project_callback)

build_util.test_main(PATH_TO_SOURCE_DIR, sys.argv, project_test_callback)

