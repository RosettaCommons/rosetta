PATH_TO_ROOT = '../../'

import hashlib, os, string, sys, shutil
#sys.path.append(PATH_TO_ROOT + 'script/')
import build_util, xcode_util

sys.path.append( '..' )
import version

# generate new svn_version file
#os.popen( 'cd ..; python svn_version.py' )
starting_directory = os.path.basename( os.getcwd() )
os.chdir( '..' )
version.svn_version()
os.chdir( starting_directory )

# (re)generate options files
os.system( 'cd ..; ./update_options.sh' )

# (re)generate residue property files
os.system( 'cd ..; ./update_residue_properties.sh' )

def project_callback(project, project_path, project_files):

	# get keys

	group_key = xcode_util.PROJECT_KEYS[project][0]
	sources_key = xcode_util.PROJECT_KEYS[project][1]

	# read file

	xcode_filename = 'Rosetta.xcodeproj/project.pbxproj'

	if not os.path.exists( xcode_filename ):
		shutil.copyfile( xcode_filename + '.template', xcode_filename )

	lines = open(xcode_filename, 'r').readlines()

	# find the relevant sections

	build_file_begin_ln = xcode_util.find_line('/* Begin PBXBuildFile section */', lines, 0)
	build_file_end_ln = xcode_util.find_line('/* End PBXBuildFile section */', lines, build_file_begin_ln)

	file_ref_begin_ln = xcode_util.find_line('/* Begin PBXFileReference section */', lines, build_file_end_ln)
	file_ref_end_ln = xcode_util.find_line('/* End PBXFileReference section */', lines, file_ref_begin_ln)

	group_begin_ln = xcode_util.find_line('/* Begin PBXGroup section */', lines, file_ref_end_ln)
	group_end_ln = xcode_util.find_line('/* End PBXGroup section */', lines, group_begin_ln)

	proj_sources_ln = xcode_util.find_line('\t\t' + sources_key + ' /* Sources */ = {', lines, group_end_ln)
	source_files_begin_ln = xcode_util.find_line('\t\t\tfiles = (', lines, proj_sources_ln)
	source_files_end_ln = xcode_util.find_line('\t\t\t);', lines, source_files_begin_ln)

	# separate sections

	a_lines = lines[0:build_file_begin_ln + 1]
	build_file_lines = lines[build_file_begin_ln + 1:build_file_end_ln]
	b_lines = lines[build_file_end_ln:file_ref_begin_ln + 1]
	file_ref_lines = lines[file_ref_begin_ln + 1:file_ref_end_ln]
	c_lines = lines[file_ref_end_ln:group_begin_ln + 1]
	group_lines = lines[group_begin_ln + 1:group_end_ln]
	d_lines = lines[group_end_ln:source_files_begin_ln + 1]
	source_files_lines = lines[source_files_begin_ln + 1:source_files_end_ln]
	e_lines = lines[source_files_end_ln:]

	# remove current information

	xcode_util.remove_groups_and_file_refs(group_key, project, group_lines, file_ref_lines)
	xcode_util.remove_build_files_and_sources(source_files_lines, build_file_lines)

	# generate new information

	root = xcode_util.make_groups_and_file_refs(group_key, project, project, project_files)
	xcode_util.add_new_lines('../src/', root, group_lines, file_ref_lines, build_file_lines, source_files_lines)

	#write new sections

	outfile = open(xcode_filename, 'w')

	outfile.writelines(a_lines)
	outfile.writelines(build_file_lines)
	outfile.writelines(b_lines)
	outfile.writelines(file_ref_lines)
	outfile.writelines(c_lines)
	outfile.writelines(group_lines)
	outfile.writelines(d_lines)
	outfile.writelines(source_files_lines)
	outfile.writelines(e_lines)

build_util.project_main(PATH_TO_ROOT + 'source/', sys.argv, project_callback)
