from __future__ import print_function
import os, stat, shutil, sys



TEST_DIR = 'test_area/'
CXX_TEST_TEMPLATE_FORMAT = '%stest/cxxtest_main.tpl'
CXX_TEST_GEN_PART_FORMAT = 'external/cxxtest/cxxtestgen.py --part --have-eh --abort-on-fail -o %s %s'
CXX_TEST_GEN_ROOT_FORMAT = 'external/cxxtest/cxxtestgen.py --root --have-eh --abort-on-fail --error-printer --template=' + CXX_TEST_TEMPLATE_FORMAT + ' -o %s %s'



def list_project_files(path_to_mini, project_name):
    path_to_project = path_to_mini + 'src/'
    settings_file_name = path_to_project + project_name + '.src.settings'
    settings_dict = {}

    with open(settings_file_name) as f:
        code = compile(f.read(), settings_file_name, 'exec')
        exec(code, settings_dict)

    # execfile(settings_file_name, settings_dict)

    project_files = []
    for dir, srcfiles in settings_dict['sources'].items():
        if len(dir) != 0 and dir[-1] != '/':
            dir += '/'

        full_path = path_to_project + dir

        hdrfiles = os.listdir(full_path)
        hdrfiles = [hdrfile for hdrfile in hdrfiles if hdrfile.endswith('.hh') or hdrfile.endswith('.h')]

        old_srcfiles = srcfiles
        srcfiles = []
        for srcfile in old_srcfiles:
            if os.path.exists(full_path + srcfile + '.cc'):
                srcfiles.append(srcfile + '.cc')
            elif os.path.exists(full_path + srcfile + '.c'):
                srcfiles.append(srcfile + '.c')
            elif os.path.exists(full_path + srcfile ) and srcfile[-3:] == '.cu':
                continue
            else:
                raise RuntimeError('Nonexistant source file: ' + full_path + srcfile)

        dirfiles = sorted(hdrfiles + srcfiles)
        project_files.append((dir, dirfiles))

    return (path_to_project, project_files)


def list_test_files(path_to_mini, test_name):
    path_to_test = path_to_mini + 'test/'
    print(path_to_test)
    settings_file_name = path_to_test + test_name + '.test.settings'
    print(settings_file_name)
    settings_dict = {}
    execfile(settings_file_name, settings_dict)

    test_files = []
    for dir, srcfiles in settings_dict['sources'].iteritems():
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

        if test_inputs_dict.has_key(dir):
            test_inputs_dict[dir].append(filename)
        else:
            test_inputs_dict[dir] = [filename]

    test_inputs = []
    for dir, inputfiles in test_inputs_dict.iteritems():
        test_inputs.append((dir, inputfiles))

    return (path_to_test, test_files, test_inputs)

def _need_update(src, dst):
    if not os.path.exists(dst):
        return True
    else:
        src_stat = os.stat(src)
        dst_stat = os.stat(dst)

        if src_stat[stat.ST_MTIME] >= dst_stat[stat.ST_MTIME]:
            return True

    return False

def update_test_files(path_to_mini, test_path, test_files):
    for dir, files in test_files:
        if not os.path.exists(TEST_DIR + dir):
            os.makedirs(TEST_DIR + dir)

        for file in files:
            if not file[1]:
                continue

            header = test_path + dir + file[0]
            source = TEST_DIR + dir + file[1]
            root = file[2]

            if _need_update(header, source) or _need_update(CXX_TEST_TEMPLATE_FORMAT % path_to_mini, source):
                if root:
                    cmd = path_to_mini + CXX_TEST_GEN_ROOT_FORMAT % (path_to_mini, source, header)
                else:
                    cmd = path_to_mini + CXX_TEST_GEN_PART_FORMAT % (source, header)

                print(cmd)
                os.system(cmd)

def update_test_inputs(path_to_mini, test_path, test_inputs):
    for dir, inputs in test_inputs:
        if not os.path.exists(TEST_DIR + dir):
            os.makedirs(TEST_DIR + dir)

        for input in inputs:
            src = test_path + dir + input
            dst = TEST_DIR + dir + input

            if _need_update(src, dst):
                print(src + ' -> ' + dst)
                shutil.copyfile(src, dst)

def clean_test(path_to_mini, test_path, test_files, test_inputs):
    to_remove = []

    for dir, files in test_files:
        for file in files:
            if file[1]:
                to_remove.append(TEST_DIR + dir + file[1])

    for dir, inputs in test_inputs:
        for input in inputs:
            to_remove.append(TEST_DIR + dir + input)

    for rem in to_remove:
        if os.path.exists(rem):
            print('removing: ' + rem)
            os.remove(rem)



def test_main(path_to_mini, argv, project_callback = None):
    if len(argv) < 2:
        print('usage: %s [test] [action]...' % argv[0])
        print('  known tests:')
        for t in KNOWN_TESTS:
            print('    %s' % t)
        print('  known actions:')
        print('     "project": update .vsproj from .src.settings "sources"')
        print('     "files": update .cxxtest.cc files from .cxxtest.hh files')
        print('     "inputs": update input files from .src.settings "inputs"')
        print('     "clean": remove .cxxtest.cc files and inputs')
        sys.exit(-1)

    tests = argv[1:]

    #JAB - wtf is whats?
    #whats = argv[2:]
    if tests == ['all']:
        tests = KNOWN_TESTS


    for test in tests:
        print("test "+test)
        if test not in KNOWN_TESTS:

            print('unknown test: ' + test)
            sys.exit(-1)

        test_path, test_files, test_inputs = list_test_files(path_to_mini, test)
        project_callback(test, TEST_DIR, test_files)


#test_main(sys.argv[1], sys.argv[2:])
