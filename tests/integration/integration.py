#!/usr/bin/env python
# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

# For help, run with -h.  Requires Python 2.4+, like the unit tests.

import sys, thread, commands
if not hasattr(sys, "version_info") or sys.version_info < (2,4):
    raise ValueError("Script requires Python 2.4 or higher!")

import os, shutil, threading, subprocess, signal, time, re, random, datetime
import json
from os import path
from optparse import OptionParser, IndentedHelpFormatter


class NT:  # named tuple
    def __init__(self, **entries): self.__dict__.update(entries)
    def __repr__(self):
        r = '|'
        for i in dir(self):
            if not i.startswith('__'): r += '%s --> %s, ' % (i, getattr(self, i))
        return r[:-2]+'|'

Jobs = []  # Global list of NameTuples  (pid, tag, start_time, out_dir,...)

def write_runtimes(runtimes, dir):
    try:
      time_file = open(dir+'/runtimes.yaml', 'w')
      json.dump(runtimes, time_file, sort_keys=True, indent=2)
      time_file.close()
    except:
      pass # if no JSON, just forget this step



def main(argv):
    '''
A simple system for running regression tests for Rosetta.

Each test has its own subdirectory in tests/, which serves as its name.
Each test has a file named "command", which should have the command to be executed.
Variable substitution is possible using Python printf format; see examples.
To create a new test, just make a new directory with a "command" file and other needed files.
See tests/HOW_TO_MAKE_TESTS for more information.

When the tests are run, one of two things happens.
(1) If the "ref" directory does not exist, it is created.
Each subdirectory of tests/ is copied to ref/, and the test is executed there.
(2) Otherwise, the "new" directory is wiped out and re-created, and tests are run in new/ instead.
Afterwards, "diff" is used to compare each subdirectory of ref/ with the corresponding subdirectory in new/.
A test is considered passed if there are no differences, and failed otherwise.

Intended usage is to run the tests once with a clean checkout, and again after changes.
The exact results of many runs are hardware-dependent, so "expected" results cannot be kept in SVN.
If tests fail in expected ways, the appropriate subdirectories can be copied from new/ to ref/ to update the expected results.

EXAMPLES:

rm -r ref/; ./integration.py    # create reference results using only default settings

./integration.py    # test reference results against new results

./integration.py -d ~/minidb -j2    # again, but using 2 processors and custom database location

./integration.py ligand_dock_7cpa    # only run the "ligand_dock_7cpa" test
    '''
    parser = OptionParser(usage="usage: %prog [OPTIONS] [TESTS]", formatter=ParagraphHelpFormatter())
    parser.set_description(main.__doc__)
    # parser.add_option("-short", ["--long"],
    #   action="store|store_true|store_false|append",
    #   default=True|False|...
    #   type="string|int|float",
    #   dest="opt_name",
    #   help="store value in PLACE",
    #   metavar="PLACE",
    # )
    parser.add_option("-d", "--database",
      default="", # processed below
      help="Path to Rosetta database. (default: ../../database, $ROSETTA_DB, ~/rosetta_database)",
    )

    parser.add_option("-m", "--mini_home",
      #default=path.join( path.expanduser("~"), "mini"),
      default= path.join( path.dirname( path.dirname( path.dirname(path.abspath(sys.argv[0])) ) ), 'source'),
      help="Directory where Mini is found (default: ../../source/)",
    )
    parser.add_option("-j", "--jobs",
      default=1,
      type="int",
      help="number of processors to use on local machine (default: 1)",
    )

    parser.add_option("--fork", action="store_true", dest="fork", default=True,
      help="Use Unix fork() instead of Python subprocess module. (off by default)"
    )

    parser.add_option("--disable-fork", action="store_false", dest="fork",
      help="Use Unix fork() instead of Python subprocess module. (off by default)"
    )

    parser.add_option("--host",
      default=[],
      action="append",
      help="ssh to HOST to run some tests. Assumes shared filesystem. May specify multiple times. May specify HOST or USER@HOST. For 5 jobs per host, HOST/5",
    )
    parser.add_option("--digs",
      default=0,
      type="int",
      help="Baker lab shortcut: use NUM processors from the digs, selected at random",
      metavar="NUM",
    )
    parser.add_option("-t", "--timeout",
      default=0,
      type="int",
      help="Maximum runtime for each test, in seconds (default: no limit)",
      metavar="SECONDS",
    )
    parser.add_option("-c", "--compiler",
      default="gcc",
      help="In selecting binaries, which compiler was used? (default: gcc)",
    )
    parser.add_option("--mode",
      default="release",
      help="In selecting binaries, which mode was used? (default: release)",
    )
    parser.add_option("--extras",
                      default="",
                      dest="extras",
                      help="in selecting binaries, which options were specified? (default: '')",
    )
    parser.add_option("--daemon", action="store_true", dest="daemon", default=False,
      help="generate daemon friendly output (off by default)"
    )
    parser.add_option("--fulldiff", action="store_true", dest="fulldiff", default=False,
      help="Include diff results for each files in to results. (off by default)"
    )
    parser.add_option("--compareonly", action="store_true", dest="compareonly", default=False,
      help="Do not run test themself, just compare results in new and ref. (off by default)"
    )
    parser.add_option("--skip-comparison", action="store_true", dest="skip_comparison", default=False,
      help="Just run test themself but do not compare results (off by default)"
    )
    parser.add_option("--yaml",
      default=None,
      help="Save results to specified file in YAML format. (default: None)",
    )

    parser.add_option("--additional_flags",
      default="",
      help="Add additional flags to integration tests. (default: None)",
    )



    parser.add_option("--dbms_database_name",
      default="rosetta_tests",
      help="For testing relational databases: the name of the database. " +
        "The variable 'test' is substituted into the parameter value. (default: rosetta_tests)",
    )

    parser.add_option("--dbms_pq_schema",
      default="integration_test_%(test)s",
      help="For testing relational databases: the name of the postgres schema. " +
        "The variable 'test' is substituted into the parameter value. (default: integration_test_%(test)s)",
    )

    parser.add_option("--dbms_host",
      default="",
      help="For testing relational databases: the name of the host. (default: '')",
    )

    parser.add_option("--dbms_user",
      default="",
      help="For testing relational databases: the name of the user, Note the password can be set with .pgpass or .my.cnf. (default: '')",
    )

    parser.add_option("--dbms_port",
      default="",
      help="For testing relational databases: the name of the user. (default: '')",
    )
    parser.add_option("--unordered",
      action="store_true",
      help="Do not order tests by expected runtime prior to launching. (default is to order)",
    )

    parser.add_option("--valgrind",
      default=False, action="store_true",
      help="Enable valgrind checking mode, instead of regular integration testing.",
    )
    parser.add_option("--valgrind_path",
      default=None,
      help="The path to the valgrind executable. (default: find in path)",
    )
    parser.add_option("--trackorigins",
      action="store_true",
      help="Tell Valgrind to output information about where uninitialized variables came from (default is no tracking, which is faster)",
    )
    parser.add_option("--leakcheck",
      action="store_true",
      help="Tell Valgrind to output details about memory leaks in addition to finding uninitialized variables. (default is no info, which is faster)",
    )

    (options, remaining_args) = parser.parse_args(args=argv)

    # Strip off whitespace and blank remaining arguments
    args = [a.strip() for a in remaining_args if a.strip()]

    options.mini_home = path.abspath( options.mini_home )

    if options.digs > 0:
        options.num_procs = 0 # don't use local processors too, b/c local *is* a dig
        # From the "freewhip" script on the whips themselves:
        # digs = ["dig01","dig02","dig03","dig04","dig05","dig06","dig07","dig08", "dig09","dig11","dig12","dig13","dig14","dig15","dig16", "dig17","dig19","dig20","dig21","dig22","dig23", "dig25","dig26"] * 2
        digs = [ "dig"+str(x) for x in range(1,33) ] * 4 # we use up to 4 processors per box
        assert options.digs <= len(digs)
        random.shuffle(digs)
        options.host.extend( digs[:options.digs] )

    #Parse database (will update options object)
    if parse_database( options, parser ) != 0: # Returns 1 on error
        return 1

    #Parse valgrind options (will update options object)
    if not options.valgrind and ( options.trackorigins or options.leakcheck or options.valgrind_path is not None ):
        print "Valgrind specific options set, but Valgrind mode not enabled. Enabling it for you."
        options.valgrind = True
    if options.valgrind:
        parse_valgrind_options(options)

    print 'Using Rosetta source dir at: ', options.mini_home
    print 'Using Rosetta database dir at:', options.database
    if options.valgrind:
        print 'Using Valgrind at:', options.valgrind_path

    global Options;  Options = options
    Options.num_procs = Options.jobs

    # Make sure the current directory is the script directory:
    # Using argv[] here causes problems when people try to run the script as "python integration.py ..."
    #os.chdir( path.dirname(sys.argv[0]) ) # argv[0] is the script name
    if not path.isdir("tests"):
        print "You must run this script from rosetta/tests/integration/"
        return 2

    # Where are we running the tests?
    if options.valgrind:
        outdir = "valgrind"
        rename_to_ref = False
    else:
        # Always put the output in the 'new' directory, and move to ref later, if appropriate
        outdir = "new"
        if not path.isdir("ref"): rename_to_ref = True
        else: rename_to_ref = False

    # Each test consists of a directory with a "command" file in it.
    if len(args) > 0:
        tests = args
    else:
        tests = [ d for d in os.listdir("tests") if not d.startswith(".") and path.isdir(path.join("tests", d)) ]

    if not options.unordered:
        tests = order_tests(tests)

    if not options.compareonly:
        if len(args)>0 and path.isdir(outdir): #we have individual tests... only remove single test directories
            for test in tests:
                testdir=outdir+'/'+test
                if path.isdir(testdir): shutil.rmtree(testdir)
        else:
        # Remove everything in the current outdir, then re-create it empty
            if path.isdir(outdir): shutil.rmtree(outdir)
            os.mkdir(outdir)

    runtimes={}
    if not options.compareonly:
        queue = Queue()
        queue.TotalNumberOfTasks = len(tests)

        # Write substitution parameters to result directory
        with open(path.join( outdir, "test_parameters.json"), "w") as parameters_file:
            json.dump(generateIntegrationTestGlobalSubstitutionParameters(), parameters_file, sort_keys=True, indent=2)

        for test in tests:
            queue.put(test)
            #shutil.copytree( path.join("tests", test), path.join(outdir, test) )
            copytree( path.join("tests", test), path.join(outdir, test),
                accept=lambda src, dst: path.basename(src) != '.svn' )

        if options.fork or options.jobs==1:
            simple_job_running( generateTestCommandline, queue, outdir, runtimes, options )
        else:
            parallel_job_running(  generateTestCommandline, queue, outdir, runtimes, options )

    # removing absolute paths to root Rosetta checkout from tests results and replacing it with 'ROSETTA'
    rosetta_dir = os.path.abspath('../..')
    for test in tests:
        for dir_, _, files in os.walk( path.join(outdir, test) ):
            for f in files:
                if f == 'command.sh': continue
                fname = dir_ + '/' + f
                data = file(fname).read()
                if rosetta_dir in data:
                    with file(fname, 'w') as f: f.write( data.replace(rosetta_dir, 'ROSETTA_MAIN') )

    # Analyze results
    print

    #if outdir == "ref":
    if rename_to_ref:
        os.renames(outdir, 'ref')

        print "Just generated 'ref' results [renamed '%s' to 'ref'];  run again after making changes." % outdir
        if options.daemon:
            print "SUMMARY: TOTAL:%i PASSED:%i FAILED:%i." % (len(tests), len(tests), 0)

        if options.yaml:
            f = file(options.yaml, 'w');  f.write("{total : %s, failed : 0, details : {}}" % len(tests));  f.close()
        write_runtimes(runtimes, 'ref')

    else:
        if options.skip_comparison:
            print 'Skipping comparison/analysis phase because command line option "--skip-comparison" was specified...'

        else:
            errors = 0
            results = {}
            full_log = ''
            for test in tests:
                if options.valgrind:
                    errors += analyze_valgrind_test(test, outdir, results, full_log )
                else:
                    errors += analyze_integration_test(test, outdir, results, full_log)

            if options.daemon:
                print "SUMMARY: TOTAL:%i PASSED:%i FAILED:%i." % (len(tests), len(tests)-errors, errors)
            else:
                if errors:
                    if options.valgrind:
                        print "%i test(s) failed.  Examine respective valgrind.out file(s) for details." % errors
                    else:
                        print "%i test(s) failed.  Use 'diff' to compare results." % errors
                else:
                    print "All tests passed."

            if options.yaml:
                try:
                  data = dict(total=len(tests), failed=errors, details=results, brief=makeBriefResults(full_log).decode('utf8', 'replace'))
                  f = file(options.yaml, 'w')
                  json.dump(data, f, sort_keys=True, indent=2)
                  f.close()
                  '''
                  f = file(options.yaml, 'w')
                  brief = makeBriefResults(full_log)
                  brief = brief.replace('"', '\\"')
                  brief = '"' + brief.replace('\n', '\\n') + '"'
                  f.write("{total : %s, failed : %s, details : %s, brief : %s}" % (len(tests), errors, results, brief) )
                  f.close()
                  '''
                except:
                  pass

        if not options.compareonly: write_runtimes(runtimes, outdir)
        if not options.valgrind:
            #compare_times has hardcoded new/ref dependancies
            from compare_times import compare_times
            compare_times(verbose=False)
    return 0

# Parse the path to the database in the options object
def parse_database( options, option_parser ):
    if not path.isdir( options.database ):
        if options.database == option_parser.get_default_values().database:
            if os.environ.get('ROSETTA3_DB') is not None and \
                    path.isdir(os.environ.get('ROSETTA3_DB')):
                options.database = os.environ.get('ROSETTA3_DB')
            else:  options.database = path.join( path.dirname( path.dirname( path.dirname(path.abspath(sys.argv[0])) ) ), 'database')

            if not path.isdir( options.database ):
                options.database = path.join( path.expanduser("~"), "rosetta_database")

            if not path.isdir( options.database ):
                print "Can't find database at %s; please set $ROSETTA3_DB or use -d" % options.database
                return 1

    # Normalize path before we change directories!
    options.database = path.abspath(options.database)
    return 0

def parse_valgrind_options(options):
    ## By default timeout is off for integration tests - observe set value
    #if options.timeout > 0 :
    #    print "Turning off timeout with valgrind." # Valgrind runs take a long time.
    #    options.timeout = 0
    if options.valgrind_path is None:
        import distutils.spawn
        options.valgrind_path = distutils.spawn.find_executable('valgrind')
        if options.valgrind_path is None:
            print "Unable to find valgrind - install or specify the path with the --valgrind_path option."
            sys.exit(1)
    options.valgrind_path = path.abspath( options.valgrind_path )
    if not os.path.exists( options.valgrind_path ):
        print "Cannot find Valgrind at", options.valgrind_path, "install or use the --valgrind_path option."
        sys.exit(1)

#
# Generate brief version of results, only ~20 first lines of difference will be shown.
#
def makeBriefResults(results):
    def replace_fun(match):  # return only first 20 lines from match object
        s = match.group()
        #return '[' + s + ']'
        lines = s.split('\n')
        if len(lines) > 40:
            lines = lines[:40] + ['---- diff output truncated ----\n']
        return '\n'.join( lines )

    r = re.compile( r'FAIL .*?^$', re.DOTALL | re.MULTILINE)
    res = r.sub(replace_fun, results)
    return res


# Wrap new line simbols inside given strings by using '\' character
def wrapNewLine(s):
    r = ''

#
# Order tests based on decreasing expected runtime. Unknown tests get run first.
# Expected runtime is taken from the runtimes.yaml file in the ref directory,
# if it exists, else from the file "runtimes" in the current directory.
# (Which is just an old version of the integration.py output.)
#
def order_tests(tests):
    times = {}
    try:
        times = json.load( open('ref/runtimes.yaml') )
    except Exception:
        pass
    if not times:
        try:
            f = open('runtimes','r')
        except Exception:
            return tests #Any problems? Just skip ordering
        try:
            for line in f:
                line = line.split()
                if len(line) > 5 and line[0] == "Finished" and line[4] == "seconds":
                    times[ line[1] ] = int(line[3])
        finally:
            f.close()
    #Decorate, sort, undecorate (A side effect is we'll alphabetize any missing tests)
    ordered = [ ( times.get(test, 9999), test ) for test in tests ]
    ordered.sort(reverse=True)
    return [ test for (time, test) in ordered ]

# -------------------------------------
def execute(message, command_line, return_=False, untilSuccesses=False, print_output=True, verbose=True):
    if verbose:
        print message
        print command_line

    while True:
        #(res, output) = commands.getstatusoutput(commandline)

        po = subprocess.Popen(command_line+ ' 1>&2', bufsize=0, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        #po = subprocess.Popen(command_line+ ' 1>&2', bufsize=0, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        f = po.stderr
        output = ''
        for line in f:
            #po.poll()
            if print_output: print line,
            output += line
            sys.stdout.flush()
        f.close()
        while po.returncode is None: po.wait()
        res = po.returncode
        #print '_____________________', res

        if res and untilSuccesses: pass  # Thats right - redability COUNT!
        else: break

        print "Error while executing %s: %s\n" % (message, output)
        print "Sleeping 60s... then I will retry..."
        time.sleep(60)

    if res:
        if print_output: print "\nEncounter error while executing: " + command_line
        if not return_: sys.exit(1)

    if return_ == 'output': return output
    else: return res

def print_(msg, color=None, background=None, bright=False, blink=False, action='print', endline=True):
    ''' print string with color and background. Avoid printing and return results str instead if action is 'return'. Also check for 'Options.no_color'
    '''
    colors = dict(black=0, red=1, green=2, yellow=3, blue=4, magenta=5, cyan=6, white=7)  # standard ASCII colors

    if 'Options' in globals()  and  hasattr(Options, 'color')  and  not Options.color: s = str(msg)
    else:
        s  = ['3%s' % colors[color] ] if color else []
        s += ['4%s' % colors[background] ] if background else []
        s += ['1'] if bright else []
        s += ['5'] if blink else []

        s = '\033[%sm%s\033[0m%s' % (';'.join(s), msg, '\n' if endline else '')

    if action == 'print': sys.stdout.write(s)
    else: return s


def mFork(times, tag=None, overhead=0, **args):
    ''' Check if number of child process is below Options.jobs. And if it is - fork the new pocees and return its pid.
    '''
    #print_('Groups:%s' % os.getgroups(), color='cyan')
    while len(Jobs) >= Options.jobs + overhead:
        for j in Jobs[:] :
            r = os.waitpid(j.pid, os.WNOHANG)
            if r == (j.pid, 0):  # process have ended without error
                normal_finish = getattr(j, 'normal_finish', lambda x: None)
                normal_finish(j, times)
                Jobs.remove(j)

            elif r[0] == j.pid :
                error_finish = getattr(j, 'error_finish', lambda x: None)
                error_finish(j, times)
                Jobs.remove(j)

            else:
                #pass
                if j.timeout:
                    if time.time() - j.start_time > j.timeout :
                        #print '~~~~~~~~~~~ pids:', j.pid, os.getpid(), os.getppid()
                        #print '~~~~~~~~~~~ groups:', os.getpgid(j.pid), os.getpgrp()
                        os.kill(j.pid, signal.SIGKILL)  #
                        #os.killpg(os.getpgid(j.pid), signal.SIGKILL)
                        timeout_finish = getattr(j, 'timeout_finish', lambda x: None)
                        timeout_finish(j, times)
                        Jobs.remove(j)
                        break

        if len(Jobs) >= Options.jobs + overhead: time.sleep(.2)

    sys.stdout.flush();  sys.stderr.flush();
    pid = os.fork()
    if pid: pass # We are parent!
    Jobs.append( NT(times=times, pid=pid, tag=tag, start_time=time.time(), **args) )
    return pid, Jobs[-1]


def mWait(tag=None, all_=False, timeout=0):
    ''' Wait for process tagged with 'tag' for completion
    '''
    while True :
        #print 'Waiting for %s: ' % tag, Jobs

        for j in [ x for x in Jobs if x.tag==tag or all_==True]:
            #print 'Waiting2: ', Jobs
            #try:
            r = os.waitpid(j.pid, os.WNOHANG)
            times = j.times
            if r == (j.pid, 0):  # process have ended without error
                normal_finish = getattr(j, 'normal_finish', lambda x: None)
                normal_finish(j, times)
                Jobs.remove(j)

            elif r[0] == j.pid :  # process ended but with error, special case we will have to wait for all process to terminate and call system exit.
                error_finish = getattr(j, 'error_finish', lambda x: None)
                error_finish(j, times)
                Jobs.remove(j)

            elif j.timeout:
                if time.time() - j.start_time > j.timeout :
                    os.kill(j.pid, signal.SIGKILL)  #os.killpg(os.getpgid(j.pid), signal.SIGKILL)
                    timeout_finish = getattr(j, 'timeout_finish', lambda x: None)
                    timeout_finish(j, times)
                    Jobs.remove(j)

        time.sleep(.2)
        if not Jobs: return


def generateIntegrationTestGlobalSubstitutionParameters(host=None):
    # Variables that may be referrenced in the cmd string:
    python = sys.executable
    minidir = Options.mini_home
    database = Options.database

    bin = path.join(minidir, "bin")
    pyapps = path.join(minidir, "src", "python", "apps")

    if sys.platform.startswith("linux"):
        platform = "linux" # can be linux1, linux2, etc
    elif sys.platform == "darwin":
        platform = "macos"
    elif sys.platform == "cygwin":
        platform = "cygwin"
    else:
        platform = "_unknown_"

    compiler = Options.compiler
    mode = Options.mode
    extras = Options.extras if Options.extras else 'default'
    binext = extras+"."+platform+compiler+mode
    #print 'binext: %s, extras: %s, Options.extras: %s' % (binext, repr(extras), repr(Options.extras) )

    additional_flags = Options.additional_flags
    dbms_host = Options.dbms_host
    dbms_user = Options.dbms_user
    dbms_port = Options.dbms_port

    return dict(locals())

def generateIntegrationTestSubstitutionParameters(test, outdir, host=None):
    """ Generate substitution parameters for integration command generation."""

    params = generateIntegrationTestGlobalSubstitutionParameters(host)
    params["dbms_database_name"] = Options.dbms_database_name % { 'test': test }
    params["dbms_pq_schema"] = Options.dbms_pq_schema % { 'test': test }
    params["workdir"] = path.abspath( path.join(outdir, test) )

    return params

def generateTestCommandline(test, outdir, options=None, host=None):
    ''' Generate and write command.sh and return command line that will run given integration test
    '''
    # Read the command from the file "command"
    params = generateIntegrationTestSubstitutionParameters(test, outdir, host)
    workdir = params["workdir"]

    if options.valgrind:
        # We need to adjust the "bin" variable to use valgrind instead
        preamble = options.valgrind_path
        if( options.trackorigins ):
            preamble = preamble + " --track-origins=yes"
        if( options.leakcheck ):
            preamble = preamble + " --leak-check=full"
        params["bin"] = preamble + " " + params["bin"]

    cmd=''
    # A horrible hack b/c SSH doesn't honor login scripts like .bash_profile
    # when executing specific remote commands.
    # This causes problems with e.g. the custom Python install on the Whips.
    # So we replace the default remote PATH with the current local one.
    if host is not None:
      cmd = 'PATH="%s"\n%s' % (os.environ["PATH"], cmd)
    cmd += '\n'
    cmd += file(path.join(workdir, "command")).read().strip()
    cmd = cmd % params # variable substitution using Python printf style

    if options.valgrind:
        # We need to remove the existance testing commands from the commandfile
        # Substitute all occurances of the "[ -x blah blah blah ]" pattern with the "true" command
        cmd = re.sub( r'\[ -x[^\]]*\]', 'true', cmd )

        # We also don't need the standard output restrictions (as we're not doing comparisons
        cmd = re.sub( r'egrep -vf ../../ignore_list', 'cat', cmd )

    cmd_line_sh = path.join(workdir, "command.sh")
    f = file(cmd_line_sh, 'w');  f.write(cmd);  f.close() # writing back so test can be easily re-run by user lately...
    #if "'" in cmd: raise ValueError("Can't use single quotes in command strings!")
    #print cmd; print

    return cmd_line_sh, workdir

def analyze_integration_test( test, outdir, results, full_log ):
    """Look at the specific integration test, and check for errors"""
    dir_before = path.join("ref", test)
    dir_after = path.join(outdir, test)
    # diff returns 0 on no differences, 1 if there are differences

    flags = ["-rq"]
    if Options.fulldiff: flags = ["-r"]
    flags += ["--exclude=command.sh"]

    proc = subprocess.Popen(["diff"] + flags + [dir_before, dir_after], stdout=subprocess.PIPE)

    full_log_msg = "FAIL %s\n" % test

    if Options.daemon:
        msg = "FAIL %s #####" % test

        if Options.fulldiff:
            for diff in proc.stdout.readlines():
                msg += "~ %s\n" % diff.strip()
                full_log_msg += "     %s\n" % diff.strip()
        else :
            lines = proc.stdout.readlines()
            for line in lines :
                cols = line.split()
                if len(cols) < 4 :
                    msg += "~ %s\n" % line.strip()
                    full_log_msg += "     %s\n" % line.strip()
                else:
                    msg += "~ nonempty diff %s %s\n" % ( cols[1], cols[3].strip() )
                    full_log_msg += "     nonempty diff %s %s\n" %  ( cols[1], cols[3].strip() )
        msg += "#####"

    else:
        msg = "FAIL %s\n" % test
        #for diff in proc.stdout.readlines():
        #    msg += "     %s\n" % diff.strip()
        #    full_log_msg += "     %s\n" % diff.strip()

        if Options.fulldiff:
            for diff in proc.stdout.readlines():
                msg += "~ %s\n" % diff.strip()
                full_log_msg += "     %s\n" % diff.strip()
        else :
            lines = proc.stdout.readlines()
            for line in lines :
                cols = line.split()
                if len(cols) < 4 or not( cols[2] == "and" ) :
                    msg += "    %s\n" %  line.strip()
                    full_log_msg += "     %s\n" % line.strip()
                else:
                    msg += "    nonempty diff %s %s\n" % ( cols[1], cols[3].strip() )
                    full_log_msg += "     nonempty diff %s %s\n" %  ( cols[1], cols[3].strip() )
    full_log_msg += '\n'

    result = proc.wait()
    results[test] = result

    if result == 0:
        print "ok   %s" % test
        full_log += "ok   %s\n" % test
        return 0
    else:
        #runtimes[test] = float('nan')
        print msg
        full_log += full_log_msg
        return 1

def analyze_valgrind_test( test, outdir, results, full_log ):
    valgrind_output = []
    dir = path.join(outdir, test)

    # Find all the Valgrind formatted lines in the log files
    def recurse( dirname, outlines ):
        for fn in os.listdir( dirname ): #Looking for all the files
            fn = os.path.join( dirname, fn )
            if os.path.isdir( fn ):
                recurse( fn, outlines )
                continue
            if (not os.path.isfile( fn )) or fn.endswith("valgrind.out"):
                continue
            f = open( fn )
            try:
                for line in f:
                    if line.startswith("=="):
                        outlines.append(line)
            finally:
                f.close()

    recurse(dir, valgrind_output)

    # Save the valgrind output specifically
    f = open( os.path.join(dir, "valgrind.out" ), 'w' )
    try:
        f.writelines( valgrind_output )
    finally:
        f.close()

    # Check that the number of log files with valgrind output matches the number we expect from
    # the number of commands we ran (i.e. there isn't a missing log somewhere)
    number_expected = 0
    f = open( os.path.join( dir, "command.sh" ) )
    try:
        for line in f:
            # Lines for valgrind runs start in column 0.
            # To turn off counting a line (e.g. because it's surrounded in an if clause)
            # Simply indent it.
            if line.startswith( Options.valgrind_path ):
                number_expected += 1
    finally:
        f.close()

    # Count the summary error lines, and sum the number of errors
    number_logs = 0
    total_errors = 0
    for line in valgrind_output:
        if line.find("ERROR SUMMARY") == -1:
            continue
        number_logs += 1
        total_errors += int(line.split()[3])

    msg = ''
    if number_logs != number_expected:
        msg += " --Log file(s) missing (%d of %d) -- " % (number_logs,number_expected)
    if total_errors != 0:
        msg += " Found " + str( total_errors ) + " Valgrind error(s)."

    results[test] = total_errors # I think this is what we should be attaching to the results object

    if len(msg)  == 0:
        print "ok   %s" % test
        full_log += "ok   %s\n" % test
        return 0
    else:
        print "FAIL %s: %s" % (test, msg)
        full_log += "FAIL %s: %s" % (test, msg)
        return 1

def simple_job_running( GenerateJob, queue, outdir, runtimes, options ):
    '''Using the function GenerateJob to generate the job commandlines, run the jobs in queue with the given options.

    GenerateJob signature:
    cmd_line_sh, workdir = GenerateJob(test, outdir)
    '''
    def signal_handler(signal_, f):
        print 'Ctrl-C pressed... killing child jobs...'
        for nt in Jobs:
            os.killpg(os.getpgid(nt.pid), signal.SIGKILL)

    signal.signal(signal.SIGINT, signal_handler)

    while not queue.empty():
        test = queue.get()
        if test is None: break

        cmd_line_sh, workdir = GenerateJob(test, outdir, options);

        def run(nt, times):
            #execute('Running Test %s' % test, 'bash ' + cmd_line_sh)
            extra = 'ulimit -t %s && ' % Options.timeout  if Options.timeout else ''
            res = execute('Running Test %s' % nt.test, '%sbash %s' % (extra, cmd_line_sh), return_=True)
            if res:
                error_string = "*** Test %s did not run!  Check your --mode flag and paths. [%s]\n" % (test, datetime.datetime.now())
                file(path.join(nt.workdir, ".test_did_not_run.log"), 'w').write(error_string)
                print error_string,
                times[test] = float('nan')

            #execute('Just sleeeping %s...' % test, 'ulimit -t%s && sleep 60 && echo "Done!"' % Options.timeout)
            #execute('Just echo %s...' % test, 'echo "%s Done!"' % test)
            #print 'Not even echo %s... ' % test

        def normal_finish(nt, times):
            queue.task_done()
            percent = (100* (queue.TotalNumberOfTasks-queue.qsize())) / queue.TotalNumberOfTasks
            elapse_time = time.time() - nt.start_time
            print "Finished %-40s in %3i seconds\t [~%4s test (%s%%) started, %4s in queue, %4d running]" % (nt.test, elapse_time, queue.TotalNumberOfTasks-queue.qsize(), percent, queue.qsize(), queue.unfinished_tasks-queue.qsize() )
            if nt.test not in times:
                times[nt.test] = elapse_time


        def error_finish(nt, times):
            error_string = "*** Test %s did not run!  Check your --mode flag and paths. [%s]\n" % (nt.test, datetime.datetime.now())
            file(path.join(nt.workdir, ".test_did_not_run.log"), 'w').write(error_string)
            print error_string,
            times[nt.test] = float('nan')
            normal_finish(nt, times)

        def timeout_finish(nt, times):
            error_string = "*** Test %s exceeded the timeout=%s  and will be killed! [%s]\n" % (nt.test, Options.timeout, datetime.datetime.now())
            file(path.join(nt.workdir, ".test_got_timeout_kill.log"), 'w').write(error_string)
            print error_string,
            times[nt.test] = float('inf')
            normal_finish(nt, times)

        if options.jobs > 1:
            pid, nt = mFork(times=runtimes, test=test, workdir=workdir, queue=queue, timeout=Options.timeout, normal_finish=normal_finish, error_finish=error_finish, timeout_finish=timeout_finish)
            if not pid:  # we are child process
                signal.signal(signal.SIGINT, signal.SIG_DFL)
                run(nt,runtimes)
                sys.exit(0)
        else:
            nt = NT(times=runtimes, test=test, workdir=workdir, queue=queue, start_time=time.time(), timeout=Options.timeout, normal_finish=normal_finish, error_finish=error_finish, timeout_finish=timeout_finish)
            run(nt,runtimes)
            if nt.timeout and (time.time() - nt.start_time > nt.timeout): nt.timeout_finish(nt,runtimes)
            else: normal_finish(nt,runtimes)

    mWait(all_=True)  # waiting for all jobs to finish before movinf in to next phase

def parallel_job_running(GenerateJob, queue, outdir, runtimes, options ):
    # Start worker thread(s)
    for i in range(options.num_procs):
        worker = Worker(queue, outdir, options, times=runtimes, timeout_minutes=options.timeout, GenerateJob=GenerateJob)
        thread = threading.Thread(target=worker.work)
        #thread.setDaemon(True) # shouldn't be necessary here
        thread.start()
    for host in options.host:
        if host.count('/') > 1:
          sys.exit("only one forward slash per host specification")
        parts = host.split('/')
        nodes=None
        if len(parts) == 1:
          nodes=1
        if len(parts) == 2:
          host= parts[0]
          nodes= int(parts[1])
        for node in range(nodes):
          worker = Worker(queue, outdir, options, times=runtimes, host=host, timeout_minutes=options.timeout,GenerateJob=GenerateJob)
          thread = threading.Thread(target=worker.work)
          #thread.setDaemon(True) # shouldn't be necessary here
          thread.start()

    # Wait for them to finish
    queue.join()

class Worker:
    def __init__(self, queue, outdir, opts, times, host=None, timeout_minutes=0, GenerateJob=generateTestCommandline):
        self.queue = queue
        self.outdir = outdir
        self.opts = opts
        self.host = host
        self.timeout = timeout_minutes * 60
        self.times = times
        self.GenerateJob = GenerateJob

    def work(self):
        running=0
        try:
            while True:
                test = self.queue.get_nowait()
                try: # Actually catch exception and ignore it.  Python 2.4 can't use "except" and "finally" together.
                    start = time.time() # initial guess at start time, in case of exception
                    try: # Make sure job is marked done even if we throw an exception
                        cmd_line_sh, workdir = self.GenerateJob(test, self.outdir, options=self.opts, host=self.host)

                        if self.host is None:
                            print "Running  %-40s on localhost ..." % test
                            proc = subprocess.Popen(["bash",  cmd_line_sh], preexec_fn=os.setpgrp)
                        # Can't use cwd=workdir b/c it modifies *local* dir, not remote dir.
                        else:
                            print "Running  %-40s on %20s ..." % (test, self.host)
                            bash_cmd='bash '+cmd_line_sh
                            proc = subprocess.Popen(["ssh", self.host, bash_cmd], preexec_fn=os.setpgrp)#, cwd=workdir)
                            #os._exit(os.EX_IOERR)
                        start = time.time() # refined start time
                        if self.timeout == 0:
                            retcode = proc.wait() # does this block all threads?
                        else:
                            while time.time() - start <= self.timeout:
                                retcode = proc.poll()
                                if retcode is not None: break
                                time.sleep(1)
                            if retcode is None:
                                print "*** Test %s exceeded the timeout and will be killed! [%s]\n" % (test, datetime.datetime.now())
                                self.times[test] = float('inf')
                                #os.kill(proc.pid, signal.SIGTERM)
                                os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
                        if retcode != 0 and retcode is not None:
                            self.times[test] = float('nan')
                            if self.host is None:
                              error_string = "*** Test %s did not run on host %s!  Check your --mode flag and paths. [%s]\n" % (test, 'local_host', datetime.datetime.now())
                            else:
                              error_string = "*** Test %s did not run on host %s!  Check your --mode flag and paths. [%s]\n" % (test, self.host, datetime.datetime.now())
                            print error_string,

                            # Writing error_string to a file, so integration test should fail for sure
                            file(path.join(workdir, ".test_did_not_run.log"), 'w').write(error_string)

                    finally: # inner try
                        percent = (100* (self.queue.TotalNumberOfTasks-self.queue.qsize())) / self.queue.TotalNumberOfTasks
                        elapse_time = time.time() - start
                        print "Finished %-40s in %3i seconds\t [~%4s test (%s%%) started, %4s in queue, %4d running]" % (test, elapse_time, self.queue.TotalNumberOfTasks-self.queue.qsize(), percent, self.queue.qsize(), self.queue.unfinished_tasks-self.queue.qsize() )
                        if test not in self.times: self.times[test] = elapse_time
                        self.queue.task_done()

                except Exception, e: # middle try
                    print e
        except Empty: # outer try
            pass # we're done, just return


class ParagraphHelpFormatter(IndentedHelpFormatter):
    '''
    A help formatter that respects paragraph breaks (blank lines) in usage strings.
    '''
    def _format_text(self, text):
        paragraphs = re.split('\n([ \t]*\n)+', text)
        paragraphs = [ IndentedHelpFormatter._format_text(self, p.strip()) for p in paragraphs ]
        return '\n'.join(paragraphs) # each already ends in a newline

def copytree(src, dst, symlinks=False, accept=lambda srcname, dstname: True):
    """Recursively copy a directory tree using copy2(), with filtering.
    Copied from shutil so I could filter out .svn entries.
    """
    names = os.listdir(src)
    os.makedirs(dst)
    errors = []
    for name in names:
        srcname = os.path.join(src, name)
        dstname = os.path.join(dst, name)
        if not accept(srcname, dstname): continue
        try:
            if symlinks and os.path.islink(srcname):
                linkto = os.readlink(srcname)
                os.symlink(linkto, dstname)
            elif os.path.isdir(srcname):
                copytree(srcname, dstname, symlinks, accept)
            else:
                shutil.copy2(srcname, dstname)
            # XXX What about devices, sockets etc.?
        except (IOError, os.error), why:
            errors.append((srcname, dstname, str(why)))
        # catch the Error from the recursive copytree so that we can
        # continue with other files
        except shutil.Error, err:
            errors.extend(err.args[0])
    try:
        shutil.copystat(src, dst)
    #except WindowsError:
        # can't copy file access times on Windows
        pass
    except OSError, why:
        errors.extend((src, dst, str(why)))
    if errors:
        raise shutil.Error, errors

################################################################################
# Python 2.4 lacks support for join() / task_done() in the Queue class,
# so I pasted the 2.5 implementation here.
# With 2.5+, you can just do "from Queue import *" instead.

from time import time as _time
from collections import deque

class Empty(Exception):
    "Exception raised by Queue.get(block=0)/get_nowait()."
    pass

class Full(Exception):
    "Exception raised by Queue.put(block=0)/put_nowait()."
    pass

class Queue:
    """Create a queue object with a given maximum size.

    If maxsize is <= 0, the queue size is infinite.
    """
    def __init__(self, maxsize=0):
        try:
            import threading
        except ImportError:
            import dummy_threading as threading
        self._init(maxsize)
        # mutex must be held whenever the queue is mutating.  All methods
        # that acquire mutex must release it before returning.  mutex
        # is shared between the three conditions, so acquiring and
        # releasing the conditions also acquires and releases mutex.
        self.mutex = threading.Lock()
        # Notify not_empty whenever an item is added to the queue; a
        # thread waiting to get is notified then.
        self.not_empty = threading.Condition(self.mutex)
        # Notify not_full whenever an item is removed from the queue;
        # a thread waiting to put is notified then.
        self.not_full = threading.Condition(self.mutex)
        # Notify all_tasks_done whenever the number of unfinished tasks
        # drops to zero; thread waiting to join() is notified to resume
        self.all_tasks_done = threading.Condition(self.mutex)
        self.unfinished_tasks = 0

    def task_done(self):
        """Indicate that a formerly enqueued task is complete.

        Used by Queue consumer threads.  For each get() used to fetch a task,
        a subsequent call to task_done() tells the queue that the processing
        on the task is complete.

        If a join() is currently blocking, it will resume when all items
        have been processed (meaning that a task_done() call was received
        for every item that had been put() into the queue).

        Raises a ValueError if called more times than there were items
        placed in the queue.
        """
        self.all_tasks_done.acquire()
        try:
            unfinished = self.unfinished_tasks - 1
            if unfinished <= 0:
                if unfinished < 0:
                    raise ValueError('task_done() called too many times')
                self.all_tasks_done.notifyAll()
            self.unfinished_tasks = unfinished
        finally:
            self.all_tasks_done.release()

    def join(self):
        """Blocks until all items in the Queue have been gotten and processed.

        The count of unfinished tasks goes up whenever an item is added to the
        queue. The count goes down whenever a consumer thread calls task_done()
        to indicate the item was retrieved and all work on it is complete.

        When the count of unfinished tasks drops to zero, join() unblocks.
        """
        self.all_tasks_done.acquire()
        try:
            while self.unfinished_tasks:
                self.all_tasks_done.wait()
        finally:
            self.all_tasks_done.release()

    def qsize(self):
        """Return the approximate size of the queue (not reliable!)."""
        self.mutex.acquire()
        n = self._qsize()
        self.mutex.release()
        return n

    def empty(self):
        """Return True if the queue is empty, False otherwise (not reliable!)."""
        self.mutex.acquire()
        n = self._empty()
        self.mutex.release()
        return n

    def full(self):
        """Return True if the queue is full, False otherwise (not reliable!)."""
        self.mutex.acquire()
        n = self._full()
        self.mutex.release()
        return n

    def put(self, item, block=True, timeout=None):
        """Put an item into the queue.

        If optional args 'block' is true and 'timeout' is None (the default),
        block if necessary until a free slot is available. If 'timeout' is
        a positive number, it blocks at most 'timeout' seconds and raises
        the Full exception if no free slot was available within that time.
        Otherwise ('block' is false), put an item on the queue if a free slot
        is immediately available, else raise the Full exception ('timeout'
        is ignored in that case).
        """
        self.not_full.acquire()
        try:
            if not block:
                if self._full():
                    raise Full
            elif timeout is None:
                while self._full():
                    self.not_full.wait()
            else:
                if timeout < 0:
                    raise ValueError("'timeout' must be a positive number")
                endtime = _time() + timeout
                while self._full():
                    remaining = endtime - _time()
                    if remaining <= 0.0:
                        raise Full
                    self.not_full.wait(remaining)
            self._put(item)
            self.unfinished_tasks += 1
            self.not_empty.notify()
        finally:
            self.not_full.release()

    def put_nowait(self, item):
        """Put an item into the queue without blocking.

        Only enqueue the item if a free slot is immediately available.
        Otherwise raise the Full exception.
        """
        return self.put(item, False)

    def get(self, block=True, timeout=None):
        """Remove and return an item from the queue.

        If optional args 'block' is true and 'timeout' is None (the default),
        block if necessary until an item is available. If 'timeout' is
        a positive number, it blocks at most 'timeout' seconds and raises
        the Empty exception if no item was available within that time.
        Otherwise ('block' is false), return an item if one is immediately
        available, else raise the Empty exception ('timeout' is ignored
        in that case).
        """
        self.not_empty.acquire()
        try:
            if not block:
                if self._empty():
                    raise Empty
            elif timeout is None:
                while self._empty():
                    self.not_empty.wait()
            else:
                if timeout < 0:
                    raise ValueError("'timeout' must be a positive number")
                endtime = _time() + timeout
                while self._empty():
                    remaining = endtime - _time()
                    if remaining <= 0.0:
                        raise Empty
                    self.not_empty.wait(remaining)
            item = self._get()
            self.not_full.notify()
            return item
        finally:
            self.not_empty.release()

    def get_nowait(self):
        """Remove and return an item from the queue without blocking.

        Only get an item if one is immediately available. Otherwise
        raise the Empty exception.
        """
        return self.get(False)

    # Override these methods to implement other queue organizations
    # (e.g. stack or priority queue).
    # These will only be called with appropriate locks held

    # Initialize the queue representation
    def _init(self, maxsize):
        self.maxsize = maxsize
        self.queue = deque()

    def _qsize(self):
        return len(self.queue)

    # Check whether the queue is empty
    def _empty(self):
        return not self.queue

    # Check whether the queue is full
    def _full(self):
        return self.maxsize > 0 and len(self.queue) == self.maxsize

    # Put a new item in the queue
    def _put(self, item):
        self.queue.append(item)

    # Get an item from the queue
    def _get(self):
        return self.queue.popleft()
################################################################################


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
