from __future__ import print_function
import os, stat, shutil, sys, os.path
from functools import reduce

def execfile(fname, settings):
    with open(fname) as f:
        code = compile(f.read(), fname, 'exec')
        exec(code, settings)


# Load projects.settings file and put resultant data into PROJECT_SETTINGS
# This script (build_util.py) is in rosetta_source/cmake/, and want rosetta_source/projects.settings
PROJECTS_SETTINGS = {}
execfile(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../projects.settings") ,PROJECTS_SETTINGS)
print(PROJECTS_SETTINGS["projects"]["external"])
KNOWN_PROJECTS = PROJECTS_SETTINGS["projects"]["src"]
KNOWN_TESTS =PROJECTS_SETTINGS["projects"]["test"]
KNOWN_EXTERNAL =PROJECTS_SETTINGS["projects"]["external"]

def list_project_files(path_to_source_dir, project_name, look_for_extension = 'all'):
    path_to_project = path_to_source_dir + 'src/'
    settings_file_name = path_to_project + project_name + '.src.settings'

    # allow user to override 'normal' .settings file with ,settings.my file.
    if len( look_for_extension ) > 0:
        settings_file_name_with_extension = settings_file_name + '.'+ look_for_extension
        if os.path.exists( settings_file_name_with_extension ):
            settings_file_name = settings_file_name_with_extension
            print('Found {} and going to use it!'.format(settings_file_name))

    # brutal hack for pilot_apps.src.settings.all
    if project_name == 'pilot_apps' and not os.path.exists( settings_file_name ): settings_file_name += '.all'

    settings_dict = {}
    execfile(settings_file_name, settings_dict)

    project_files = []
    for dir, srcfiles in settings_dict['sources'].items():
        if len(dir) != 0 and dir[-1] != '/':
            dir += '/'

        hash_index = project_name.find('#')
        if hash_index == -1:
            project_name_nohash = project_name
        else:
            project_name_nohash = project_name[:hash_index]
        #full_path = path_to_project + '/' + project_name_nohash + '/' + dir
        full_path = path_to_project + '/' + dir

        # another brutal hack for pilot_apps.src.settings.all
        if project_name == 'pilot_apps' or project_name == 'apps':
            full_path = path_to_project + '/apps/' + '/' + dir

        old_srcfiles = srcfiles
        srcfiles = []
        for srcfile in old_srcfiles:
            if len(srcfile) > 3 and srcfile[-3:] == ".cu" :
                continue # APL NOTE: cannot currently compile cuda files with cmake
            if os.path.exists(full_path + srcfile + '.cc'):
                srcfiles.append(srcfile + '.cc')
            elif os.path.exists(full_path + srcfile + '.c'):
                srcfiles.append(srcfile + '.c')
            else:
                raise RuntimeError('Nonexistent source file: ' + full_path + srcfile)

        if not os.path.exists(full_path):
            # Git ignores empty directories, so if there aren't any headers or source files the listdir() below may choke.
            # A missing directory where something was expected should error out above.
            continue

        hdrfiles = os.listdir(full_path)
        hdrfiles = [hdrfile for hdrfile in hdrfiles if hdrfile.endswith('.hh') or hdrfile.endswith('.h')]

        dirfiles = sorted(hdrfiles + srcfiles)
        project_files.append((dir, dirfiles))

    #return (path_to_project + project_name_nohash + '/', project_files)
    return (path_to_project + '/', project_files)

def list_test_files(path_to_source_dir, test_name):
    path_to_test = path_to_source_dir + 'test/'
    settings_file_name = path_to_test + test_name + '.test.settings'
    settings_dict = {}
    execfile(settings_file_name, settings_dict)

    test_files = []
    for dir, srcfiles in settings_dict['sources'].items():
        if len(dir) != 0 and dir[-1] != '/':
            dir += '/'
        dir = test_name + '/' + dir

        dirfiles = []
        for file in srcfiles:
            dirfiles.append((file + '.cxxtest.hh', file + '.cxxtest.cc', False))

        hdrfiles = os.listdir(path_to_test + '/' + dir)
        hdrfiles = [hdrfile for hdrfile in hdrfiles if hdrfile.endswith('.hh') and not hdrfile.endswith('.cxxtest.hh') ]
        for file in hdrfiles:
            dirfiles.append((file, None, None))

        test_files.append((dir, dirfiles))

    dirfiles = []
    dirfiles.append((test_name + '.cxxtest.hh', test_name + '.cxxtest.cc', True))

    hdrfiles = os.listdir(path_to_test + '/' + test_name)
    hdrfiles = [hdrfile for hdrfile in hdrfiles if hdrfile.endswith('.hh') and not hdrfile.endswith('.cxxtest.hh') ]
    for file in hdrfiles:
        dirfiles.append((file, None, None))

    test_files.append((test_name + '/', dirfiles))

    test_inputs_dict = {}
    for input in settings_dict['testinputfiles']:
        filename = os.path.basename(input)
        dir = os.path.dirname(input)
        if len(dir) != 0 and dir[-1] != '/':
            dir += '/'
        dir = test_name + '/' + dir

        if dir in test_inputs_dict:
            test_inputs_dict[dir].append(filename)
        else:
            test_inputs_dict[dir] = [filename]

    test_inputs = []
    for dir, inputfiles in test_inputs_dict.items():
        test_inputs.append((dir, inputfiles))

    return (path_to_test, test_files, test_inputs)

def list_external_files(path_to_mini, external_name):
    path_to_external = path_to_mini + '/external/'
    settings_file_name = path_to_external + external_name + '.external.settings'
    if not os.path.exists( settings_file_name ):
        return ( "", [] , [] )
    settings_dict = {}
    execfile(settings_file_name, settings_dict)

    project_files = []
    for dir, srcfiles in settings_dict['sources'].items():
        if len(dir) != 0 and dir[-1] != '/':
            dir += '/'

        full_path = path_to_external + dir

        if not os.path.exists(full_path) and 'only_with_extras' in settings_dict and settings_dict['only_with_extras']:
            # This may be part of a submodule which isn't currently loaded, but isn't needed for the current only_with_extras build
            # Don't crash (yet).
            continue

        old_srcfiles = srcfiles
        srcfiles = []
        for srcfile in old_srcfiles:
            if srcfile.endswith( ".cu" ):
                continue # APL NOTE: cannot currently compile cuda files with cmake
            if os.path.exists( full_path + srcfile):
                    srcfiles.append(srcfile)
            elif os.path.exists(full_path + srcfile + '.cc'):
                srcfiles.append(srcfile + '.cc')
            elif os.path.exists(full_path + srcfile + '.c'):
                srcfiles.append(srcfile + '.c')
            else:
                raise RuntimeError('Nonexistent source file: ' + full_path + srcfile)

        if not os.path.exists(full_path):
            # Git ignores empty directories, so if there aren't any headers or source files the listdir() below may choke.
            # A missing directory where something was expected should error out above.
            continue

        hdrfiles = os.listdir(full_path)
        hdrfiles = [hdrfile for hdrfile in hdrfiles if hdrfile.endswith('.hh') or hdrfile.endswith('.h')]

        dirfiles = sorted(hdrfiles + srcfiles)
        project_files.append((dir, dirfiles))

    return (path_to_external, project_files, settings_dict)


def update_libraries_list(projects):
    exclude = ['apps', 'pilot_apps']
    out = open( 'build/libraries.cmake', 'w')
    out.write('SET( LIBRARIES ' + reduce(lambda x,y: x + ' ' + y, filter(lambda x: x not in exclude, projects)) + ')')

def update_test_list(projects, test_path):
    out = open('build/test_libraries.cmake', 'w')
    try:
        out.write('SET( TEST_LIBRARIES ' + ' '.join( projects) + ')\n')
        out.write('SET( SRCDIR ' + test_path + ' )\n' )
    finally:
        out.close()

def update_externals_list(projects):

    with open('build/external_libraries.cmake', 'w') as out:

        for project, other_settings in projects.items():
            project_entry = \
                'SET(EXTERNAL_LIBRARIES ${EXTERNAL_LIBRARIES} %s)\n' % project

            if other_settings.get("only_with_extras", False):
                # Add ifdef guards around the external module definition to check
                # for each required extra.
                guard = ""
                endguard = ""

                for e in other_settings["only_with_extras"]:
                    guard += 'if( ${EXTRAS} MATCHES "%s" )\n' % e
                    endguard += 'endif()\n'

                project_entry = "\n" + guard + project_entry + endguard + "\n"

            out.write(project_entry)

def test_main(path_to_source_dir, argv, project_test_callback = None):

    tests =  [ "all" ]

    if tests == [ "all" ]:
        tests = KNOWN_TESTS

    for test in tests:

        if test not in KNOWN_TESTS:
            print('unknown test: ' + test)
            sys.exit(-1)

        test_path, test_files, test_inputs = list_test_files(path_to_source_dir, test)
        project_test_callback(test, path_to_source_dir, "build/", test_files, test_inputs)

    update_test_list(tests, path_to_source_dir + "test/")

def external_main(path_to_mini, argv, project_external_callback = None):
    externals =  [ "all" ]

    if externals == [ "all" ]:
        externals = KNOWN_EXTERNAL

        buildable_externals = dict()
    for external in externals:
        print("Examining external" + external)
        if external not in KNOWN_EXTERNAL:
            print('unknown external project: ' + external)
            sys.exit(-1)

        project_path, project_files, other_settings = list_external_files(path_to_mini, external)
        if project_files:
            can_build = project_external_callback(external, project_path, project_files, other_settings)
            if can_build:
                buildable_externals[external] = other_settings

    update_externals_list(buildable_externals)


def project_main(path_to_source_dir, argv, project_callback):
    if len(argv) < 2:
        print('usage: {} [project]...'.format(argv[0]))
        print('  known projects:')
        for p in KNOWN_PROJECTS + ['all','my']:
            print('    {}'.format(p))
        sys.exit(-1)
    update_libraries_list(KNOWN_PROJECTS)


    projects = argv[1:]
    if projects == ['all']:
        projects = KNOWN_PROJECTS

    look_for_extension = 'all'
    if projects == ['my']:
        look_for_extension = 'my'
        projects = KNOWN_PROJECTS

    for project in projects:
        if project not in KNOWN_PROJECTS:
            print('unknown project: ' + project)
            sys.exit(-1)

        project_path, project_files = list_project_files(path_to_source_dir, project, look_for_extension)
        project_callback(project, project_path, project_files)


