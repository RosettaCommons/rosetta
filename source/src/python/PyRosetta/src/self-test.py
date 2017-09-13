#!/usr/bin/env python3.5
# was !/usr/bin/env python
# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   self-test.py
## @brief  Run bindings test and demo scrips
## @author Sergey Lyskov

from __future__ import print_function

import os, os.path, sys, json, datetime, subprocess, argparse, time, glob, signal, shutil

def execute(message, command_line, return_='status', until_successes=False, terminate_on_failure=True, silent=False):
    print(message);  print(command_line)
    while True:

        p = subprocess.Popen(command_line, bufsize=0, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, errors = p.communicate()

        output = output.decode('utf-8', errors="replace") + errors.decode('utf-8', errors="replace")
        exit_code = p.returncode

        if exit_code  or  not silent: print(output)

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
    else: return False

_jobs_ = []
def mfork():
    ''' Check if number of child process is below Options.jobs. And if it is - fork the new pocees and return its pid.
    '''
    while len(_jobs_) >= Options.jobs :
        for p in _jobs_[:] :
            r = os.waitpid(p, os.WNOHANG)
            if r == (p, 0):  # process have ended without error
                _jobs_.remove(p)
            elif r[0] == p :  # process ended but with error, special case we will have to wait for all process to terminate and call system exit.
                for p in _jobs_: os.waitpid(p, 0)
                print('Some of the unit test suite terminate abnormally!')
                sys.exit(1)

        if len(_jobs_) >= Options.jobs: time.sleep(.5)

    sys.stdout.flush()
    pid = os.fork()
    if pid: _jobs_.append(pid) # We are parent!
    return pid


_test_output_ = '.test.output/'  # dir to store test results

def test_name(test): return os.path.basename(test)[:-3]

def json_file_name(test): return _test_output_+ '.test.' + test_name(test) + '.json'


def run_test(test):
    name = test_name(test)
    json_file = json_file_name(test)

    if os.path.isfile(json_file): os.remove(json_file)

    started = datetime.datetime.today()

    command_line = 'SET PYTHONPATH=%CD%;%PYTHONPATH%' if sys.platform == "win32" else 'export PYTHONPATH=`pwd`:$PYTHONPATH && ulimit -t 4096'
    command_line += ' && {0} {1} '.format(sys.executable, test)

    res, output = execute('\nExecuting %s...' % test, command_line, return_='tuple')
    run_time = '\nFinished {0} in {1}'.format(name, datetime.datetime.today() - started)
    print(run_time)

    output += run_time

    with open(json_file, 'w') as f: json.dump({ name : dict(log=output, state='failed' if res else 'passed') }, f, sort_keys=True, indent=2)
    sys.stdout.flush()


def main(args):
    parser = argparse.ArgumentParser(usage="Generate pyrosetta distribution.")

    parser.add_argument('-j', '--jobs',
      default=1, type=int,
      help="Number of processors to use on when building. (default: use all avaliable memory)",
    )

    parser.add_argument("--delete-tests-output", action="store_true", default=False, help="Do not run tests, instead delete tests output files and exit. [off by default]" )

    parser.add_argument("--enable-gui",
        help = "Also run GUI tests",
        default = False,
        action = "store_true")

    parser.add_argument("--gui-only",
        help = "Run only GUI tests",
        default = False,
        action = "store_true")

    parser.add_argument('args', nargs=argparse.REMAINDER)

    global Options;
    Options = parser.parse_args(args=args[1:])


    def get_py_files(dir_):
        return [dir_ + '/' + f for f in os.listdir(dir_) if f.endswith('.py')  and  f!='__init__.py'  and  \
                ( f.startswith('G') if Options.gui_only else ( True if Options.enable_gui else not f.startswith('G')) )  ]

    tests = Options.args or sorted( get_py_files('test') + get_py_files('demo') )

    print('Preparing to run:\n{}'.format('\n'.join(tests) ) )

    # removing previous output files if any...
    if os.path.isdir(_test_output_): print( 'Removing old test dir {0}...'.format(_test_output_) );  shutil.rmtree(_test_output_)  # remove old dir if any
    if Options.delete_tests_output: return
    os.mkdir(_test_output_)

    json_file = _test_output_ + '.test.results.json'

    if Options.jobs>1:
        def signal_handler(signal_, f):
            print('Ctrl-C pressed... killing child jobs...')
            for j in _jobs_:
                os.killpg(os.getpgid(j), signal.SIGKILL)

        signal.signal(signal.SIGINT, signal_handler)


    for t in tests:
        if Options.jobs > 1:
            pid = mfork()
            if not pid:  # we are child process
                run_test(t)
                sys.exit(0)

        else:
            run_test(t)

    for p in _jobs_: os.waitpid(p, 0)  # waiting for all child process to termintate...


    results = dict(tests={})
    state = 'passed'
    for t in tests:
        test_results = json.load( open( json_file_name(t) ) )
        if 'state' not in test_results[ test_name(t) ]  or  test_results[ test_name(t) ]['state']!= 'passed': state = 'failed'
        results['tests'].update(test_results)

    with open(json_file, 'w') as f: json.dump(dict(state=state, results=results, log=''), f, sort_keys=True, indent=2)

    for t in sorted(results['tests']):
        if results['tests'][t]['state'] == 'passed': print(t, '- passed')

    failed_tests = [ t for t in sorted(results['tests']) if results['tests'][t]['state'] != 'passed' ]

    if failed_tests:
        print('\nFollowing PyRosetta Tests FAILED:')
        for t in failed_tests: print(t, '- FAILED')

    if state == 'passed': print('\nAll PyRosetta Tests passed!\n')
    else: sys.exit(1)


if __name__ == "__main__":
    print('{}::__main__, started at {}...'.format(sys.argv[0], datetime.datetime.now()) )
    main(sys.argv)
