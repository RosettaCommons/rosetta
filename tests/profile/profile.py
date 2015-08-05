#!/usr/bin/env python
# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
#
# @author Sergey Lyskov
#

import os, sys, os.path, time, commands, subprocess, datetime
from optparse import OptionParser
from shutil import move
import json

# Create global 'Platform' that will hold info of current system
if sys.platform.startswith("linux"): Platform = "linux" # can be linux1, linux2, etc
elif sys.platform == "darwin" : Platform = "macos"
else: Platform = "_unknown_"
#PlatformBits = platform.architecture()[0][:2]


def getSubprocessMemoryInfo(parentPID):
    size = 0
    for line in commands.getoutput('ps --ppid %s -o vsize,pid' % parentPID).split('\n')[1:]:
        sz, pid = map(int, line.split())
        size += sz + getSubprocessMemoryInfo(pid)

    return size




def run(test, options):
    workdir = os.path.abspath( os.path.join("tests", test) )

    # Running tests
    platform = Platform

    minidir = os.path.join(options.main, "source")

    bindir = os.path.join(options.main, "source", "bin")
    compiler = 'gcc'
    mode = 'release'
    binext = platform+compiler+mode

    print 'Running test %s...' % test
    print '  Test working dir is: %s' % workdir
    print '  Rosetta home dir is: %s' % minidir
    print '  Rosetta bin dir is: %s' % bindir
    print '  Rosetta database dir is: %s' % options.database

    additional_flags = options.additional_flags

    templates = dict(minidir=minidir, database=options.database, workdir=workdir, platform=platform, bin=bindir, compiler=compiler, mode=mode, binext=binext, additional_flags=additional_flags)

    fname = os.path.join(workdir, 'command')
    cmd = file(fname).read().strip()
    cmd = cmd % templates  # variable substitution using Python printf style
    # Writing result to .sh file for future refference.
    f = file(fname+'.sh', 'w');  f.write(cmd);  f.close()

    output_dir = os.path.join(workdir, 'output')
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    # Now we actualy run test...
    proc = subprocess.Popen(["bash", fname+'.sh'], preexec_fn=os.setpgrp)
    retcode=None
    start = time.time() # refined start time
    #if self.timeout == 0:
    # retcode = proc.wait() # does this block all threads?
    timeout = options.timeout
    if timeout == 0 : timeout = 99999999

    memory = [(0, 0)]
    time.sleep(1)

    while time.time() - start <= timeout:
        retcode = proc.poll()
        if retcode is not None: break
        memory.append( (int(time.time() - start), getSubprocessMemoryInfo(os.getpid())/1000. ) )
        #print 'Time: %s, Memory: %s' % memory[-1]
        time.sleep(1)

    # writing results to a file and generating memory graph
    memory_data_fn = output_dir + '/memory.data'
    f = file(memory_data_fn, 'w')
    f.write('%-10s %-10s\n' % ('time', 'memory'));  map(lambda x: f.write('%-10s %-10s\n' % x), memory);  f.close()

    print commands.getoutput("""echo 'set terminal png small;set xlabel \"time sec\";set ylabel \"memory MB\";plot "%s" u 1:2 notitle' | gnuplot > %s/memory.png""" % (memory_data_fn, output_dir))

    max_memory_allocated = max( map(lambda x: x[1], memory) )
    execution_time = memory[-1][0]

    filename = output_dir + '/.results.yaml'
    old_filename = output_dir+'/.old_results.yaml'
    if os.path.isfile(filename): move(filename, output_dir+'/.old_results.yaml')
    if options.daemon  and  os.path.isfile(old_filename): os.remove(old_filename)

    yaml_data = { 'execution_time_s' : execution_time, 'max_memory_allocated_MB': max_memory_allocated }
    f = file(filename, 'w');  json.dump( yaml_data, f );  f.close()


    failed = False

    if retcode is None:
        import signal
        print "*** Test %s exceeded the timeout and will be killed!" % test
        os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
        failed = True

    if retcode != 0 and retcode is not None: failed = True

    if failed:
        error_string = "*** Test %s did not run!  Check your --mode flag and paths. [%s]\n" % (test, datetime.datetime.now())
        print error_string,
        with file(os.path.join(workdir, ".test_did_not_run.log"), 'w') as f: f.write(error_string)  # Writing error_string to a file, so regression test should fail for sure

    return dict(state='failed' if failed else 'finished', max_memory_allocated=max_memory_allocated, execution_time=execution_time, memory_usage=[ dict(time=m[0], memory=m[1]) for m in memory])




def main(argv):
    '''A simple system for running protocol profile end-to-end tests in Mini.
    '''

    json_results_file = '.profile_test_results.json'
    if os.path.isfile(json_results_file): os.remove(json_results_file)


    # Attempt to resolve checkout root directory.
    try:
        main_dir = subprocess.check_output(["git", "rev-parse", "--show-toplevel"]).strip()
    except (subprocess.CalledProcessError, OSError) as e:
        main_dir = os.path.abspath('./../../')

    parser = OptionParser(usage="usage: %prog [OPTIONS] [TESTS]")
    parser.set_description(main.__doc__)

    parser.add_option("-m", "--main",
      default=main_dir,
      help="Rosetta main directory. (default: %s)" % (main_dir if main_dir else "None"))

    parser.add_option("-d", "--database",
      default=None,
      help="Directory where rosetta database is located. (default: <main>/database)",
    )

    parser.add_option("--additional_flags",
      default="",
      help="Add additional flags to porfile tests. (default: None)",
    )

    parser.add_option("-t", "--timeout",
      default=0,
      type="int",
      help="Maximum runtime for each test, in minutes (default: no limit)",
      metavar="MINUTES",
    )

    parser.add_option("--daemon",
      action="store_true",
      dest="daemon",
      help="Run tests in deamon mode. This will skip some steps and void generation of file with old results.",
    )


    (options, args) = parser.parse_args(argv)

    if not options.main:
        print "Unable to resolve rosetta main repository root, run profile.py from within rosetta main checkout or specify root with --main."
        return 1

    if not os.path.isdir(options.main):
        print "Invalid rosetta main repository root %s, run profile.py from within rosetta main checkout or specify root with --main." % options.main
        return 1

    if not options.database:
        options.database = os.path.join(options.main, "database")

    if not os.path.isdir(options.database):
        print "Can't find database at %s; please use -d" % options.database
        return 1

    # Normalize path before we change directories!
    options.main = os.path.abspath(options.main)
    options.database = os.path.abspath(options.database)

    # Switch to test directory
    os.chdir(os.path.join(options.main, "tests", "profile"))

    # Each test consists of a directory with a "command" file in it.
    if len(args) > 0:
        tests = args
    else:
        tests = [ d for d in os.listdir("tests") if not d.startswith(".") and os.path.isdir(os.path.join("tests", d)) ]


    tests_results = {}
    # Now actually running the tests...
    for test in tests:
        tests_results[test] = run(test, options)

    with file(json_results_file, 'w') as f: json.dump(tests_results, f, sort_keys=True, indent=2)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
