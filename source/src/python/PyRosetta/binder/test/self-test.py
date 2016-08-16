#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true: :collapseFolds=10:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   self-test.py
## @brief  Script to run Binder tool self-test
## @author Sergey Lyskov


import os, sys, argparse, platform, subprocess, shutil

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
    else: return exit_code


def get_test_files(dir_):
        return [dir_ + '/' + f for f in os.listdir(dir_) if f.startswith('T') and f.endswith('.hpp')]


def run_test(test_path, build_dir):
    if not os.path.isfile(test_path): print('Could not find test file: {} ... Exiting...'.format(test_path)); sys.exit(1)

    test = os.path.basename(test_path)
    test_name = test[:-len('.hpp')]
    source_dir = os.path.dirname(test_path)

    source_include = build_dir + '/' + test + '.include'
    with open(source_include, 'w') as f: f.write( '#include <{}>\n'.format(test) )

    root_module = test[:-len('.hpp')].replace('.', '_')
    config = test_path[:-len('.hpp')] + '.config'
    if not os.path.isfile(config): config = ''

    command_line = '{binder} --bind "" --root-module {root_module} --prefix {build_dir} --single-file --annotate-includes {config} {source} -- -x c++ -std=c++11 -I {source_dir} -I {source_dir}/..'.format(binder=Options.binder, root_module=root_module, build_dir=build_dir, source_dir=source_dir, source=source_include,
                                                                                                                                                                                         config='--config {}'.format(config) if config else '')

    execute('{} Running test...'.format(test), command_line);

    command_line = 'cd {build_dir} && clang++ -O3 -shared -std=c++11 -I {pybind11} -I/usr/include/python2.7 -I./.. -I./../.. {root_module}.cpp -o {root_module}.so -fPIC'.format(pybind11=Options.pybind11, root_module=root_module, build_dir=build_dir)
    execute('{} Compiling binding results...'.format(test), command_line);

    command_line = "cd {build_dir} && python -c 'import {root_module}'".format(root_module=root_module, build_dir=build_dir)
    execute('{} Testing imports...'.format(test), command_line);

    ref = source_dir + '/' + test_name + '.ref'
    new = build_dir + '/' + root_module + '.cpp'
    command_line = 'diff {ref} {new}'.format(**vars())

    r = execute('Comparing results for test {}...'.format(test), command_line, terminate_on_failure=not Options.accept);

    not_binded = [l for l in open(new).read().split('\n') if 'not_binded' in l]
    if not_binded: print('{}\n"not_binded" string found in results for test {}!!!\n'.format('\n'.join(not_binded), test)); sys.exit(1)

    if r  and  Options.accept:
        a = input( 'Accept new results from test {} as reference? (yes or no) '.format(test) )
        if a == 'yes': shutil.copyfile(new, ref)



def main():
    ''' Binder tool self-test script '''

    parser = argparse.ArgumentParser()

    parser.add_argument('--binder', default='binder', help='Path to Binder tool')

    parser.add_argument('--pybind11', default='', help='Path to pybind11 source tree')

    parser.add_argument("--accept", action="store_true", help="Run tests and accept new tests results as reference")

    parser.add_argument('args', nargs=argparse.REMAINDER)

    global Options
    Options = parser.parse_args()

    test_source_dir = os.path.dirname( os.path.abspath(__file__) )

    tests = Options.args or sorted( get_test_files(test_source_dir) )

    build_dir = test_source_dir + '/build/'
    if os.path.isdir(build_dir): print('Removing old test dir {0} ...'.format(build_dir));  shutil.rmtree(build_dir)  # remove old dir if any
    os.makedirs(build_dir)
    #if not os.path.isdir(build_dir): os.makedirs(build_dir)

    for t in tests: run_test(t, build_dir)

    print('Done!')


if __name__ == "__main__":
    main()
