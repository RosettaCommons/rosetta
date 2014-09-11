#!/usr/bin/env python
# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   TestBindings.py
## @brief  Run bindings test and demo scrips
## @author Sergey Lyskov


import os, os.path, sys, json, commands, datetime, subprocess, argparse, time, glob, signal


def execute(message, commandline, return_=False, untilSuccesses=False):
    print message, commandline

    while True:
        if sys.platform == "win32":
            po = subprocess.Popen(commandline, bufsize=0, shell=True, stdout=sys.stdout, stderr=sys.stderr)
            while po.returncode is None: po.wait()
            res = po.returncode

        else:
            (res, output) = commands.getstatusoutput(commandline)
            print output

        if res and untilSuccesses: pass  # Thats right - redability COUNT!
        else: break

        print "Error while executing %s: %s\n" % (message, output)
        print "Sleeping 60s... then I will retry..."
        time.sleep(60)

    if return_ == 'tuple': return res, output

    if res:
        print "\nEncounter error while executing: " + commandline
        if return_==True: return True
        else: sys.exit(1)

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
                print 'Some of the unit test suite terminate abnormally!'
                sys.exit(1)

        if len(_jobs_) >= Options.jobs: time.sleep(.5)
    pid = os.fork()
    if pid: _jobs_.append(pid) # We are parent!
    return pid


def test_name(test): return os.path.basename(test)[:-3]

def json_file_name(test): return '.test.' + test_name(test) + '.json'


def run_test(test):
    name = test_name(test)
    json_file = json_file_name(test)

    if os.path.isfile(json_file): os.remove(json_file)

    started = datetime.datetime.today()

    command_line = 'SET PYTHONPATH=%CD%;%PYTHONPATH%' if sys.platform == "win32" else 'export PYTHONPATH=`pwd`:$PYTHONPATH'
    command_line += ' && {0} {1} '.format(sys.executable, test)

    res, output = execute('\nExecuting %s...' % test, command_line, return_='tuple')
    run_time = '\nFinished {0} in {1}'.format(name, datetime.datetime.today() - started)
    print run_time

    output += run_time

    with file(json_file, 'w') as f: json.dump({ name : dict(log=output, state='failed' if res else 'finished') }, f, sort_keys=True, indent=2)
    sys.stdout.flush()


def main(args):
    parser = argparse.ArgumentParser(usage="Generate pyrosetta distribution.")

    parser.add_argument('-j', '--jobs',
      default=1, type=int,
      help="Number of processors to use on when building. (default: use all avaliable memory)",
    )

    parser.add_argument('args', nargs=argparse.REMAINDER)

    global Options;
    Options = parser.parse_args(args=args[1:])

    def get_py_files(dir_):
	if not os.path.exists(dir_): return []
	return [dir_ + '/' + f for f in os.listdir(dir_) if f.endswith('.py')  and  f!='__init__.py']

    tests = Options.args or sorted( get_py_files('test') + get_py_files('demo') )

    print 'Preparing to run:\n%s\n' % '\n'.join(tests)

    # removing previous output files if any...
    for f in '.test_kic.pdb _.pdb ddG_out_1.txt dna_output.fasc dna_output_1.pdb dock_output.fasc dock_output_1.pdb fold_output.fasc fold_output_1.pdb loop10.pdb loop13.pdb loop16.pdb loop23.pdb loop26.pdb loop30.pdb loop_output.fasc loop_output_1.pdb poly-A_final.pdb poly-A_low.pdb refine_output.fasc refine_output_1.pdb sample_resfile'.split():
        if os.path.isfile(f): os.remove(f)

    for f in glob.iglob('.test.*.json'):
        if os.path.isfile(f): os.remove(f)

    json_file = '.test_bindings.json'
    if os.path.isfile(json_file): os.remove(json_file)


    if Options.jobs>1:
        def signal_handler(signal_, f):
            print 'Ctrl-C pressed... killing child jobs...'
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
    state = 'finished'
    for t in tests:
        test_results = json.load( file( json_file_name(t) ) )
        if 'state' not in test_results[ test_name(t) ]  or  test_results[ test_name(t) ]['state']!= 'finished': state = 'failed'
        results['tests'].update(test_results)

    with file(json_file, 'w') as f: json.dump(dict(state=state, results=results, log=''), f, sort_keys=True, indent=2)


    if state == 'finished': print '\nAll PyRosetta Tests passed!\n'
    else:
        for t in results['tests']:
            if results['tests'][t]['state'] == 'finished': print t, '- passed'
        print '\nFollowing PyRosetta Tests FAILED:'
        for t in results['tests']:
            if results['tests'][t]['state'] != 'finished': print t, '- FAILED'



if __name__ == "__main__":
    print '%s::__main__, started at %s...' % (sys.argv[0], datetime.datetime.now())
    main(sys.argv)


#tests = execute('getting list of tests...', 'ls test/T*.py', return_= 'output').split()
#tests = filter(lambda x: x.startswith('T') and x.endswith('.py'), sorted( os.listdir('test/') ) )

#tests = filter(lambda x: x.endswith('.py'), sorted( os.listdir('test/') + os.listdir('demos/') ) )

#for t in map(lambda x: x.replace('/', '.'), args[1:]) or tests:
    #print '\nRunning %s...' % t
    #__import__( t[:-3] )
    #__import__( 'test.' + t[:-3] )
    #execute('Executing %s...' % t, 'export PYTHONPATH=`pwd`:$PYTHONPATH && python %s' % t)
    #execute('Executing %s...' % t, 'source SetPyRosettaEnvironment.sh && python %s' % t)
