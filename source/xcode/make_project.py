#!/usr/bin/env python

#Original Author: Unknown
#Last Major Edits: Jared Adolf-Bryfogle (jadolfbr@gmail.com)



import os
work_dir = os.getcwd()
if os.path.dirname(__file__): os.chdir(os.path.dirname(__file__))
script_dir = os.getcwd()

import os, string, sys, shutil
import build_util, xcode_util

KNOWN_PROJECTS = [
    'basic',
    'utility',
    'numeric',
    'ObjexxFCL',
    'core.1',
    'core.2',
    'core.3',
    'core.4',
    'core.5',
    #'interactive',
    #'game',
    ]

if 'cygwin' in sys.platform:
    KNOWN_PROJECTS.append('protocols.1')
    KNOWN_PROJECTS.append('protocols.2')
else:
    KNOWN_PROJECTS += [
    'protocols.1',
    'protocols_a.2',
    'protocols_b.2',
    'protocols.3',
    'protocols.4',
    'protocols_a.5',
    'protocols_b.5',
    'protocols_c.5',
    'protocols_d.5',
    'protocols_e.5',
    'protocols_f.5',
    'protocols_a.6',
    'protocols_b.6',
    'protocols_c.6',
    'protocols_d.6',
    'protocols_e.6',
    'protocols.7',
    'protocols.8',
        ]

#Needed by scons, but have no source files:
"""
    'protocols_g.5',
    'protocols_h.5',
"""


KNOWN_TESTS = [
    'basic',
    'apps',
    'core',
    'demo',
    'devel',
    'numeric',
    'protocols',
    'utility',
    ]

def project_callback(project, project_path, project_files):

	# get keys
	print project
	#print project_files
	group_key = xcode_util.PROJECT_KEYS[project][0]
	sources_key = xcode_util.PROJECT_KEYS[project][1]

	# read file

	xcode_filename = 'Rosetta.xcodeproj/project.pbxproj'

	#if not os.path.exists( xcode_filename ): - JAB - new libraries will fail
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

def project_main(path_to_mini, argv, project_callback, known_projects):
    if len(argv) < 2:
        print('usage: %s [project]...' % argv[0])
        print('  known projects:')
        for p in known_projects + ['all']:
            print('    %s' % p)
        sys.exit(-1)

    projects = argv[1:]
    if projects == ['all']:
        projects = known_projects

    for project in projects:
        if project not in known_projects:
            print('unknown project: ' + project)
            sys.exit(-1)

        project_path, project_files = build_util.list_project_files(path_to_mini, project)
        project_callback(project, project_path, project_files)


if __name__ == "__main__":

	main_dir = '../../'
	source_dir = '../'

	sys.path.append('..')

	# generate new version file
	os.system('cd {source_dir}; python version.py'.format(source_dir=source_dir))

	# (re)generate options files
	os.system('cd {source_dir}; ./update_options.sh'.format(source_dir=source_dir))

	# (re)generate ResidueType enum files
	os.system('cd {source_dir}; ./update_ResidueType_enum_files.sh'.format(source_dir=source_dir))

	#print repr(KNOWN_PROJECTS)

	test_projects=['protocols.3', 'protocols.4']

	project_main(main_dir + 'source/', sys.argv, project_callback, KNOWN_PROJECTS)
