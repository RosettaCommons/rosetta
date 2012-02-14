#!/usr/bin/python -u

# make_project.py - modified from mini-interactive by tex.

import os, string, sys
import timeit
import build_util

PATH_TO_ROOT = '../../'
MINI_DIR = (os.path.abspath(sys.argv[0])).split("/")[-3]

def update_svn_version():
	cmd = 'cd ..; python svn_version.py'
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
			if ( project == 'pilot_apps' ):
				project = 'apps'
			file_path = '/'.join( [ '../../src', project, entry[0] ] )
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
			output += 'ADD_CUSTOM_COMMAND( TARGET %s POST_BUILD COMMAND python ../smart_symlink.py ${COMPILER} ${MODE} %s )\n' % ( symlink_var, key ) 
			output += 'ADD_DEPENDENCIES( %s %s )\n' % ( symlink_var, key )
			output += 'ADD_DEPENDENCIES( %s BUILD_ROSETTA_LIBS )\n' % ( key )
			apps_file = key + '.cmake'
			open( 'build/apps/' + apps_file, 'w').write( output )
			apps_output.write( 'INCLUDE( ../build/apps/%s )\n' % apps_file )
		apps_output.close()
	else:
		cmake_files = ''
		for dir, files in project_files:
			cmake_files += '\n\t' + string.join(['../' + project_path + dir + file for file in files], '\n\t')

		#print 'making project files for ' + project + ' from ' + project_path

		output = ''
		output += 'SET(' + project + '_files' + cmake_files + '\n)\n'

		open('build/' + project + '.cmake', 'w').write(output)
	print 'done.'


update_svn_version()

build_util.project_main(PATH_TO_ROOT + MINI_DIR + "/", sys.argv, project_callback)

