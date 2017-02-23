#!/usr/bin/env python

PATH_TO_ROOT = '../../'

import os, string, sys
sys.path.append(PATH_TO_ROOT + 'script/')
import build_util, xcode_util

def project_callback(test, test_path, test_files):
	# set up data

	xcode_test = test + '.test'

	xcode_files = []
	for dir, files in test_files:
		header_files = []
		source_files = []
		for file in files:
			header_files.append(file[0])
			if (file[1]):
				source_files.append(file[1])

		print dir
		xcode_files.append(('test/' + dir, header_files))
		xcode_files.append(('xcode/test_area/' + dir, source_files))

	print xcode_files

	# get keys

	group_key = xcode_util.TEST_KEYS[xcode_test][0]
	sources_key = xcode_util.TEST_KEYS[xcode_test][1]

	# read file

	#xcode_filename = 'RosettaTests.xcodeproj/project.pbxproj'
	xcode_filename = 'Rosetta.xcodeproj/project.pbxproj'
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

	xcode_util.remove_groups_and_file_refs(group_key, xcode_test, group_lines, file_ref_lines)
	xcode_util.remove_build_files_and_sources(source_files_lines, build_file_lines)

	# generate new information

	root = xcode_util.make_groups_and_file_refs(group_key, xcode_test, '', xcode_files)
	xcode_util.add_new_lines('../', root, group_lines, file_ref_lines, build_file_lines, source_files_lines)

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

build_util.test_main(PATH_TO_ROOT + 'source/', sys.argv, project_callback)
