#!/usr/bin/env python
# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.
#
# For help, run with -h.  Requires Python 2.4+, like the unit tests.
# Author: Sergey Lyskov
# Author: Jared Adof-Bryfogle (demo extension)
from __future__ import print_function

import sys
if not hasattr(sys, "version_info") or sys.version_info < (2,4):
    raise ValueError("Script requires Python 2.4 or higher!")

import os, shutil, threading, subprocess, signal, time, re, random, datetime, copy, traceback
import io
import json
import glob
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
platforms = ['linux', 'macos', 'windows']
compilers = ['gcc', 'clang', 'icc']
modes     = ['release', 'debug', 'release_debug']
this_file_path = os.path.realpath(__file__)
this_file_dir = os.path.dirname(this_file_path)

#We assume a typical Rosetta install for testing purposes

# old, before submodules layout
# root_rosetta_dir = os.path.abspath(os.path.join(this_file_dir,"../../.."))
# root_main_dir =  os.path.realpath(os.path.join(root_rosetta_dir, "main"))
# root_demos_dir = os.path.realpath(os.path.join(root_rosetta_dir, "demos"))
# root_tools_dir = os.path.realpath(os.path.join(root_rosetta_dir, "tools"))

# new directory layout when demos and tools are submodules to main
root_rosetta_dir = os.path.abspath(os.path.join(this_file_dir,"../.."))
root_main_dir =  root_rosetta_dir
root_demos_dir = os.path.realpath(os.path.join(root_rosetta_dir, "demos"))
root_tools_dir = os.path.realpath(os.path.join(root_rosetta_dir, "tools"))

# ----------------------------------------------------------------

timeout_factors = {}
with open("timeout_factors.json") as f:
    timeout_factors = json.loads(f.read())

pwd = os.getcwd()

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

A special case are MPI integration tests, which are defined by command.mpi files.  These can be run with the --mpi_tests
flag (which implies extras=mpi).

EXAMPLES:

rm -r ref/; ./integration.py    # create reference results using only default settings

./integration.py    # test reference results against new results

./integration.py -d ~/minidb -j2    # again, but using 2 processors and custom database location

./integration.py ligand_dock_7cpa    # only run the "ligand_dock_7cpa" test

./integration.py tests/ligand_dock_7cpa # only run the "ligand_dock_7cpa" test

./integration.py -d ~/database -j5 --mpi-tests # run the special MPI integration tests using a custom database location and 5 processors.





EXAMPLES For Running Demos/Tutorials

./integration.py --demos

./integration.py --tutorials

./integration.py --tutorials Tutorial_4_relax

./integration.py --tutorials tutorials/Tutorial_4_relax


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
      default="gcc" if sys.platform != "darwin" else "clang",
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
      default=[], action="append",
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


    parser.add_option("--suffix",
      default=None,
      help="Specify `command` suffix. When used only tests that have `command.<suffix>` is run, with corresponding extras build. Default is None: run only plain test (command)",
    )

    parser.add_option("--mpi-tests",
      action="store_const", const='mpi', dest="suffix",
      help="[DEPRECATED, USE `--suffix mpi` INSTEAD] Run only those tests defined in command.mpi files.  This option implies --extras=mpi.",
    )



    parser.add_option("--demos",
      help="Signifiy we are testing the [public] demos.",
      default=False,
      action="store_true"
    )
    parser.add_option("--tutorials",
      help = "Signify we are testing the tutorials",
      default = False,
      action = "store_true",
    )
    parser.add_option("--demos_root_dir",
      help = "If you are running the demos or tutorials test using a non-typical Rosetta directory structure, pass the dir here.",
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
        print("Valgrind specific options set, but Valgrind mode not enabled. Enabling it for you.")
        options.valgrind = True
    if options.valgrind:
        parse_valgrind_options(options)

    global Options;  Options = options

    if Options.suffix == 'thread' and Options.jobs != 1:
        print('WARNING: `--suffix thread` is only compatible with `--jobs 1`! Setting `jobs` to `1`...')
        Options.jobs = 1

    Options.num_procs = Options.jobs

    print('Using Rosetta source dir at: ' + Options.mini_home)
    print('Using Rosetta database dir at:' + Options.database)
    if Options.demos or Options.tutorials:
        print('\nUsing Rosetta tools at: ' + root_tools_dir)
        print('Using Rosetta demos at: ' + root_demos_dir)
    if options.valgrind:
        print('Using Valgrind at:' + options.valgrind_path)

    # Make sure the current directory is the script directory:
    # Using argv[] here causes problems when people try to run the script as "python integration.py ..."
    #os.chdir( path.dirname(sys.argv[0]) ) # argv[0] is the script name
    #if not path.isdir("tests") and not Options.demos:
    #    print "You must run this script from rosetta/tests/integration/"
    #    return 2


    # Each test consists of a directory with a "command" file in it.

    #Print the current SHA1 for demos, main, and tools.

    print("\nCurrent Versions Tested:")
    print("MAIN:  " + get_git_sha1(Options.mini_home))
    if get_git_sha1(Options.mini_home) != get_git_sha1(Options.database):
        print("\n====== WARNING WARNING WARNING ========")
        print("\t Rosetta database version doesn't match source version!")
        print("\t DATABASE: " + get_git_sha1(Options.database))
        print("====== WARNING WARNING WARNING ========\n")
    print("TOOLS: " + get_git_sha1(root_tools_dir))
    print("DEMOS: " + get_git_sha1(root_demos_dir))
    print("\n")

    globalparams = generateIntegrationTestGlobalSubstitutionParameters()
    print("Python: `"  + globalparams.get("python","")  +"`")
    print("Python2: `" + globalparams.get("python2","") +"`")
    print("Python3: `" + globalparams.get("python3","") +"`")
    print("\n")

    #All tests are in a subdirectory.  We set these up here.
    test_subdir = "tests"
    demo_subdir = "public"
    tutorial_subdir = "tutorials"


    tests = []

    if Options.demos or Options.tutorials:
        os.chdir(root_demos_dir)
    else:
        os.chdir(this_file_dir)

    if len(args) > 0:
        if Options.demos:
            tests = get_tests(demo_subdir, args)
        elif Options.tutorials:
            tests = get_tests(tutorial_subdir, args)
        else:
            tests = get_tests(test_subdir, args)

    elif Options.demos or Options.tutorials:
        subdir = tutorial_subdir if Options.tutorials else demo_subdir
        tests = [f for f in glob.glob(subdir+"/*") if os.path.isdir(f)]

    else:
        tests = [path.join(test_subdir,d) for d in os.listdir(test_subdir) if not d.startswith(".") and path.isdir(path.join(test_subdir, d)) ]


    #Setup the Demos or Tutorials
    if Options.demos or Options.tutorials:
        for test_dir in copy.deepcopy(tests):
            if os.path.isdir(test_dir):
                print("Setting up: "+test_dir)
                if not setup_demo_command_file(test_dir):
                    tests.remove(test_dir)

    if len(tests) == 0:
        sys.exit("No tests present.  Exiting.")


    # Where are we running the tests?
    if options.valgrind:
        outdir = "valgrind"
        rename_to_ref = False
    else:
        # Always put the output in the 'new' directory, and move to ref later, if appropriate
        outdir = "new"
        if not path.isdir("ref"): rename_to_ref = True
        else: rename_to_ref = False

    if not options.compareonly:
        if len(args)>0 and path.isdir(outdir): #we have individual tests... only remove single test directories
            for test in tests:
                testbase = os.path.basename(test)
                testdir=os.path.join(outdir,testbase)
                if path.isdir(testdir): shutil.rmtree(testdir)
        else:
            # Remove everything in the current outdir, then re-create it empty
            if path.isdir(outdir): shutil.rmtree(outdir)
            os.mkdir(outdir)

    runtimes={}

    if not options.unordered:
        tests = order_tests(tests)


    if not options.compareonly:
        queue = Queue()
        queue.TotalNumberOfTasks = len(tests)

        # Write substitution parameters to result directory
        print("Outdir: "+outdir)
        with open(path.join( outdir, "test_parameters.json"), "w") as parameters_file:
            json.dump(globalparams, parameters_file, sort_keys=True, indent=2)

        for test in tests:
            testbase = os.path.basename(test)
            #shutil.copytree( path.join("tests", test), path.join(outdir, test) )
            if Options.demos:
                #print "Copying demo dir from "+test+" to "+ path.join(outdir, test)
                local_copytree( test , path.join(outdir, testbase), accept=lambda src, dst: path.basename(src) != '.svn')
                queue.put(testbase)
            elif (not Options.suffix)  or  os.path.isfile(path.join( test , 'command.' + Options.suffix)):
                #print "Copying " + test + " to " + outdir + "."
                local_copytree( test , path.join(outdir, testbase), accept=lambda src, dst: path.basename(src) != '.svn')
                queue.put(testbase)
            else:
                continue

        if options.fork or options.jobs==1:
            simple_job_running( generateTestCommandline, queue, outdir, runtimes, options, globalparams )
        else:
            parallel_job_running( generateTestCommandline, queue, outdir, runtimes, options, globalparams )

    # removing absolute paths to root Rosetta checkout from tests results and replacing it with 'ROSETTA'
    for test in tests:
        testbase = os.path.basename(test)
        for dir_, _, files in os.walk( path.join(outdir, testbase) ):
            for f in files:
                if f == 'command.sh': continue
                if f == 'command.mpi.sh': continue
                fname = dir_ + '/' + f
                try:
                    with io.open(fname, 'r', encoding="UTF-8", errors='backslashreplace') as f:
                        data = f.read()
                except UnicodeDecodeError:
                    # binary files will not work with the replacements
                    # below
                    continue
                except TypeError as e:
                    # Python2 and early versions of Python3 use a different error class
                    if len(e.args) > 0 and type( e.args[0] ) == str and 'UnicodeDecodeError' in e.args[0]:
                        continue

                mod = False
                if root_rosetta_dir in data:
                    data = data.replace(root_rosetta_dir, 'ROSETTA')
                    mod = True
                # In case commandline specification is used - normalize to standard install directories.
                params = globalparams
                if params['minidir'] in data:
                    data = data.replace( params['minidir'], "ROSETTA/source")
                    mod = True
                if params['database'] in data:
                    data = data.replace( params['database'], "ROSETTA/database")
                    mod = True
                if params['rosetta_tools'] in data:
                    data = data.replace( params['rosetta_tools'], "ROSETTA/tools")
                    mod = True
                if params['rosetta_demos'] in data:
                    data = data.replace( params['rosetta_demos'], "ROSETTA/demos")
                    mod = True
                if mod:
                    with io.open(fname, 'w', encoding="UTF-8", errors='backslashreplace') as f:
                        f.write( data )

    # Analyze results
    print()

    refdir = os.path.join( os.path.dirname(outdir), 'ref' )
    if rename_to_ref:
        os.renames( outdir, refdir )

        print("Just generated 'ref' results [renamed '%s' to '%s'];  run again after making changes." % (outdir, refdir))
        if options.daemon:
            print("SUMMARY: TOTAL:%i PASSED:%i FAILED:%i." % (len(tests), len(tests), 0))

        if options.yaml:
            with open(options.yaml, 'w') as f:
                f.write("{total : %s, failed : 0, details : {}}" % len(tests))
        write_runtimes(runtimes, refdir)

    else:
        if options.skip_comparison:
            print('Skipping comparison/analysis phase because command line option "--skip-comparison" was specified...')

        else:
            errors = 0
            results = {}
            full_log = []
            for test in tests:
                if Options.suffix:
                    test_valid = os.path.isfile( path.join( test, "command." + Options.suffix) )
                else:
                    test_valid = os.path.isfile( path.join( test, "command") )

                if test_valid:
                    testbase = os.path.basename(test)
                    if options.valgrind:
                        errors += analyze_valgrind_test(testbase, outdir, results, full_log )
                    else:
                        errors += analyze_integration_test(testbase, outdir, refdir, results, full_log)
            full_log = ''.join(full_log)

            if options.daemon:
                print("SUMMARY: TOTAL:%i PASSED:%i FAILED:%i." % (len(tests), len(tests)-errors, errors))
            else:
                if errors:
                    if options.valgrind:
                        print("%i test(s) failed.  Examine respective valgrind.out file(s) for details." % errors)
                    else:
                        print("%i test(s) failed.  Use 'diff' to compare results." % errors)
                else:
                    print("All tests passed.")

            if options.yaml:
                try:
                  data = dict(total=len(tests), failed=errors, details=results, brief = makeBriefResults(full_log) )
                  with open(options.yaml, 'w') as f:
                      json.dump(data, f, sort_keys=True, indent=2)
                  '''
                  f = file(options.yaml, 'w')
                  brief = makeBriefResults(full_log)
                  brief = brief.replace('"', '\\"')
                  brief = '"' + brief.replace('\n', '\\n') + '"'
                  f.write("{total : %s, failed : %s, details : %s, brief : %s}" % (len(tests), errors, results, brief) )
                  f.close()
                  '''
                except Exception as e:
                    trace = traceback.format_exc()
                    print( 'Integration script failed with exception while writing results:{}\n{}... Skipping results writing...'.format(trace, e) )


        if not options.compareonly: write_runtimes(runtimes, outdir)
        if not options.valgrind:
            #compare_times has hardcoded new/ref dependancies
            from compare_times import compare_times
            compare_times(verbose=False)

    os.chdir(pwd)
    return 0


def get_tests(subdir, args):
    """
    Get the list of tests in the subdirectory from the list of args.  Args can be full paths, relative, or just the name of the paths

    :param subdir: str
    :param args: str
    :rtype:[str]
    """
    return [path.join(subdir, os.path.basename(a)) for a in args]


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
                print("Can't find database at %s; please set $ROSETTA3_DB or use -d" % options.database)
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
            print("Unable to find valgrind - install or specify the path with the --valgrind_path option.")
            sys.exit(1)
    options.valgrind_path = path.abspath( options.valgrind_path )
    if not os.path.exists( options.valgrind_path ):
        print("Cannot find Valgrind at", options.valgrind_path, "install or use the --valgrind_path option.")
        sys.exit(1)


def setup_demo_command_file(demo_subdir):
    """
    Setup a command.sh file in the demo subdirectory if not already present.  This command.sh is written from the
     demo .md file.  A parsable command should start with the '$' charactor.

    Returns a boolean for if a command file is already present or a new one is written.

    Author: Jared Adolf-Bryfogle and the XRW 2016 Team.

    :param demo_dir: str
    :rtype: bool
    """

    files = glob.glob(demo_subdir+'/*')
    md_files = sorted([f for f in files if f.endswith(".md")])

    if "command" in [os.path.basename(f) for f in files]: return True

    commands = []
    rosetta_binaries = []
    total_commands = 0
    for md_file in md_files:

        with io.open(md_file, 'r', encoding="UTF-8", errors='backslashreplace') as f:
            md_file_data = f.read()

            if os.path.basename(md_file) == "README.md":
                demo_name = os.path.dirname(md_file).split('/')[-1]
            else:
                demo_name = os.path.basename(md_file).replace(".md", "")

            commands.append("# "+demo_name+"\n\n")
            for line in md_file_data.split('\n'):
                #print line
                line = line.strip()
                line = line.replace("```", "") #Replace code line.
                line = line.replace("<code>", "")
                line = line.replace("</code>","")


                #$ Charactor is the line for an actual command that will be tested.
                if not line or not line.startswith("$>"): continue

                line = line.replace("`", "")

                #Replace rosetta binaries with that with which they are being run on.
                #Repalce flags files with any short version present.
                line, line_exe = format_demo_line(line, demo_subdir)
                if len(line_exe) > 0:
                    total_commands+=1
                    line = line+ " -database %(database)s -run:constant_seed -nodelay 2>&1 " \
                           + "| egrep -vf "+os.path.dirname(this_file_path)+"/ignore_list > log"+str(total_commands)

                else:
                    total_commands+=1
                    if not re.search('cd', line):
                        #Parenthesis allow proper IO direction around command (for example piping the symdef files into a file)
                        # They are required.  Thanks to Rocco Moretti for figuring out how to accomplish this.

                        line = "("+line+")"+ " 2>&1 | egrep -vf "+os.path.dirname(this_file_path)+"/ignore_list > log"+str(total_commands)


                rosetta_binaries.extend(line_exe)
                commands.append(line)

    if total_commands == 0:
        return False

    test_cmd = '''test "${PIPESTATUS[0]}" != '0' && exit 1 || true  # Check if the first executable in pipe line return error and exit with error code if so'''

    OUTFILE = open(demo_subdir+"/command.new", "w")

    #These cause spurious failures with multiple redirection and other places.
    #OUTFILE.write("set -e\n")
    #OUTFILE.write("set -o pipefail\n")

    OUTFILE.write("# This is a test script of the demo, created by parsing testable commands given in the demo's README.md file.\n")
    OUTFILE.write("# Commands that need to be tested should start with '$<'  to indicate that the line will be written to this file and tested.\n")
    OUTFILE.write("# - Documentation XRW 2016 Team - \n\n")

    OUTFILE.write("cd %(workdir)s\n\n")

    for exe in sorted(set(rosetta_binaries)):
        chars = "~!@#$%^&*()`+=[]{}\|;:',<.>?" #People do some wierd stuff.
        for c in chars:
            exe = exe.replace(c, "")

        OUTFILE.write("[ -x %(bin)s/{exe}.%(binext)s ] || exit 1\n".format(exe=exe))

    OUTFILE.write("\n")

    for command in commands:
        OUTFILE.write(command+"\n")
        if not command.startswith("#"):
            OUTFILE.write(test_cmd+"\n\n")

    OUTFILE.write("\n")
    OUTFILE.close()

    return True
    #print demo_subdir+" command written"


def format_demo_line(command, demo_subdir):
    """
    Formats a demo line into that which can be run through the integration test framework.
    Returns the new line and a list of Rosetta executables found.

    :param command: str
    :rtype: str, [str]
    """
    executables = []
    short_flags = []

    new_command = command.replace("$>", "")

    #THese have to be in this order, or $ROSETTA3 will be replaced in $ROSETTA3_DB.
    new_command = new_command.replace("$ROSETTA3_DB", "%(database)s")
    new_command = new_command.replace("$ROSETTA3", "%(minidir)s")

    new_command = new_command.replace("$ROSETTA_MAIN", "%(rosetta_main)s")
    new_command = new_command.replace("$ROSETTA_TOOLS", "%(rosetta_tools)s")
    new_command = new_command.replace("$ROSETTA_DEMOS", "%(rosetta_demos)s")
    new_command = new_command.replace("$ROSETTA_BINEXT","%(binext)s")
    commandSP = new_command.split()

    for word in commandSP:
        if '.' in word:
            wordSP = os.path.basename(word).split(".")
            ext = wordSP[-1]
            if match_patterns(ext, platforms) and match_patterns(ext, compilers) and match_patterns(ext, modes):
                executables.append(wordSP[0])

                #People can give path/to/rosetta_scripts.linuxgccrelease in their tutorials and we need to match it.
                new_command = new_command.replace(word, "%(bin)s/{exe}.%(binext)s".format(exe=wordSP[0]))


            else: continue

        if word.startswith("@"):
            shorty = word[1:]+".short"

            #print "Should be adding shorty "+os.path.join(demo_subdir,shorty)
            if os.path.exists(os.path.join(demo_subdir,shorty)):
                #print "Adding Shorty: "+os.path.join(demo_subdir, shorty)
                short_flags.append('@'+shorty)

    new_command = new_command+" "+" ".join(short_flags) #These are the settings to make the demos run short.

    return new_command, executables



def to_string(b):
    ''' Conver bytes to string and handle the errors. If argument is already string - do nothing
    '''
    return b if type(b) == str else str(b, 'utf-8', errors='backslashreplace')


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
    return to_string( res )


# Wrap new line simbols inside given strings by using '\' character
def wrapNewLine(s):
    r = ''

def get_git_sha1(dir):
    """
    Get the SHA1 of the git directory.
    :param dir:
    :return:
    """
    if not os.path.exists(dir):
        return "<<Directory '{dir}' not found.>>".format(dir=dir)
    cmd = "cd {dir}; git rev-parse HEAD; cd {workdir}".format(dir=dir, workdir=os.getcwd())
    # Let's return as a string, which means that for python2 we must decode here?
    bytes_sha = subprocess.check_output(cmd, shell=True)
    return bytes_sha.decode('utf-8', errors="backslashreplace").strip()

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
    ordered = []
    # Anything extremely long in debug mode should not be run on the test server.
    tabooed = ["tests/remodel", "tests/remodel_disulfides", "tests/inverse_rotamer_remodel", "tests/pepspec", "tests/bridge_chains",
              "tests/mp_relax_w_ligand", "tests/membrane_relax", "tests/membrane_relax_hbond", "tests/continuous_sewing_hasher",
              "tests/discontinuous_sewing_hasher"]

    for test in tests:
        # skip tabooed-for-debug tests
        if test in tabooed and Options.mode != "release": continue

        testbase = os.path.basename(test)
        if test in times:
            ordered.append( (times[test], test) )
        elif testbase in times:
            ordered.append( (times[testbase], test) )
        else:
            ordered.append( (9999, test) )
    ordered.sort(reverse=True)
    #print([ test for (time, test) in ordered ])
    return [ test for (time, test) in ordered ]

# -------------------------------------
def execute(message, command_line, return_=False, untilSuccesses=False, print_output=True, verbose=True):
    if verbose:
        print(message)
        print(command_line)

    while True:
        #(res, output) = commands.getstatusoutput(commandline)

        po = subprocess.Popen(command_line, bufsize=0, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (stdout, stderr) = po.communicate()
        stdout = stdout.decode('utf-8', errors="backslashreplace")
        stderr = stderr.decode('utf-8', errors="backslashreplace")
        if print_output:
            print( stdout )
            print( stderr )

        while po.returncode is None: po.wait()
        res = po.returncode
        #print('_____________________' + res)

        if res and untilSuccesses: pass  # Thats right - redability COUNT!
        else: break

        print("Error while executing %s: %s\n" % (message, output))
        print("Sleeping 60s... then I will retry...")
        time.sleep(60)

    if res:
        if print_output: print("\nEncounter error while executing: " + command_line)
        if not return_: sys.exit(1)

    if return_ == 'output': return stdout + stderr
    elif return_ == 'stdout': return stdout
    elif return_ == 'stderr': return stderr
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
                    # Temporarily don't time out on particular tests.
                    if time.time() - j.start_time > j.timeout:
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
                if time.time() - j.start_time > j.timeout:
                    os.kill(j.pid, signal.SIGKILL)  #os.killpg(os.getpgid(j.pid), signal.SIGKILL)
                    timeout_finish = getattr(j, 'timeout_finish', lambda x: None)
                    timeout_finish(j, times)
                    Jobs.remove(j)

        time.sleep(.2)
        if not Jobs: return

def get_binext():
    """
    Get the extension for the binaries, and a dictionary of the local variables.
    :rtype: str
    """
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
    if Options.extras:
        extras = Options.extras
    elif Options.suffix:
        # Use .get() to default to whatever the suffix is for the compiler tag
        # Only need to list those which differ
        extras =  dict(thread='cxx11thread',).get( Options.suffix, Options.suffix )
    else :
        extras='default'

    binext = extras.replace(',', '')+"."+platform+compiler+mode
    return binext, dict(locals())


def verify_python_version(executable, version):
    command_line = "{executable} -c 'import sys; sys.exit(sys.version_info[0] != {version})'".format(**vars())
    execute('Verifying Python {executable} to have version {version}.*...'.format(**vars()), command_line)

def generateIntegrationTestGlobalSubstitutionParameters():
    # Variables that may be referenced in the cmd string:
    python = sys.executable
    if sys.version_info[0] == 2:
        python2 = sys.executable
        python3 = execute("Find python3","which python3",return_="stdout",print_output=False,verbose=False).strip()
    elif sys.version_info[0] == 3:
        python3 = sys.executable
        python2 = execute("Find python2","which python2",return_="stdout",print_output=False,verbose=False).strip()
        if python2 == "":
            # On the Mac test servers, there isn't a `python2`, but the regular `python` should be python2
            # however if we inside Python virtual environemnt then `which python2` migh give nothing so we will look for python2.7
            python2 = execute("Find python (assume python2.7)", "which python2.7",return_="stdout",print_output=False,verbose=False).strip()
    else:
        print("ERROR: Unrecognized Python version!")
        sys.exit(-1)

    if python2 == "":
        print("ERROR: Unable to find Python2 executable -- some integration tests may fail on that basis alone.")
        python2 = "PYTHON2_NOT_FOUND"
    else:
        verify_python_version(python2, 2)

    if python3 == "":
        print("ERROR: Unable to find Python3 executable -- some integration tests may fail on that basis alone.")
        python3 = "PYTHON3_NOT_FOUND"
    else:
        verify_python_version(python3, 3)

    minidir = Options.mini_home
    database = Options.database

    bin = path.join(minidir, "bin")
    raw_bin_dir = bin  # Duplicate so that we can avoid valgrind-related expansion when we want to.
    pyapps = path.join(minidir, "scripts", "python")
    template_dir = path.join(minidir, "code_templates")

    binext, local_vars = get_binext()
    #print 'binext: %s, extras: %s, Options.extras: %s' % (binext, repr(extras), repr(Options.extras) )

    additional_flags = ' '.join( Options.additional_flags )
    dbms_host = Options.dbms_host
    dbms_user = Options.dbms_user
    dbms_port = Options.dbms_port

    #No idea why this won't work in the function below.
    rosetta_demos = root_demos_dir
    rosetta_tools = root_tools_dir
    rosetta_main  = root_main_dir

    ##Merge locals
    local_vars.update(dict(locals()))
    del local_vars['local_vars']
    return local_vars

def generateIntegrationTestSubstitutionParameters(test, outdir, host=None, globalparams=None):
    """
    Generate substitution parameters for integration command generation.
    """

    if globalparams is None:
        params = generateIntegrationTestGlobalSubstitutionParameters()
    else:
        params = globalparams.copy()
    params["host"] = host
    params["dbms_database_name"] = Options.dbms_database_name % { 'test': test }
    params["dbms_pq_schema"] = Options.dbms_pq_schema % { 'test': test }
    params["workdir"] = path.abspath( path.join(outdir, test) )

    return params

def generateTestCommandline(test, outdir, globalparams=None, options=None, host=None):
    """
    Generate and write command.sh and return command line that will run given integration test
    """

    params = generateIntegrationTestSubstitutionParameters(test, outdir, host=host, globalparams=globalparams)
    workdir = params["workdir"]
    if options.valgrind:
        # We need to adjust the "bin" variable to use valgrind instead
        preamble = options.valgrind_path
        if options.trackorigins:
            preamble = preamble + " --track-origins=yes"
        if options.leakcheck:
            preamble = preamble + " --leak-check=full"
        # Certain tests need extended stack sizes
        NEEDS_EXPANDED_STACK = ['fragment_picker','fragmentpicker_integration_demo','ligand_dock_ensemble','repeat_propagate']
        if test in NEEDS_EXPANDED_STACK:
            preamble = preamble + " --main-stacksize=" + str( 64*1024*1024 ) # Default is ~32 MB, try ~64 MB stack
        params["bin"] = preamble + " " + params["bin"]

    cmd=u''
    # A horrible hack b/c SSH doesn't honor login scripts like .bash_profile
    # when executing specific remote commands.
    # This causes problems with e.g. the custom Python install on the Whips.
    # So we replace the default remote PATH with the current local one.
    if host is not None:
      cmd = 'PATH="%s"\n%s' % (os.environ["PATH"], cmd)
    cmd += '\n'
    if options.suffix:
        command_file_name = path.join(workdir,"command." + options.suffix)
        if os.path.isfile(command_file_name):
            with open(command_file_name) as f: cmd += f.read().strip()
        else:
            return None, None #If the command.mpi file doesn't exist, then this isn't an MPI integration test and should be skipped.
    else :
        if os.path.isfile( path.join(workdir, "command" ) ):
            with io.open(path.join(workdir, "command"), 'r', encoding="UTF-8", errors='backslashreplace') as f:
                cmd += f.read().strip()
        elif os.path.isfile( path.join(workdir, "command.new") ):
            with io.open(path.join(workdir, "command.new"), 'r', encoding="UTF-8", errors='backslashreplace') as f:
                cmd += f.read().strip()
        else:
          return None, None #If the command file doesn't exist, we skip it.  It may only be an MPI integration test.


    #Error in replacement.  We mark this test as failed and continue.
    try:
        cmd = cmd % params # variable substitution using Python printf style
    except KeyError as e:
        out =  "*** Error in command file replacement: %("+str(e).strip('\'')+")s"

        print("\n*** Test {s} did not run!  Check your --mode flag and paths.".format(s=test))
        print(out+"\n")

        LOG = open(path.join(workdir, ".test_did_not_run.log"), "w")
        LOG.write(out)
        LOG.close()
    except Exception as e:
        out =  "*** Error in command file replacement: "+str(e)

        print("\n*** Test {s} did not run!  Check your --mode flag and paths.".format(s=test))
        print(out+"\n")

        LOG = open(path.join(workdir, ".test_did_not_run.log"), "w")
        LOG.write(out)
        LOG.close()

        return None, None

    if options.valgrind:
        # We need to remove the existance testing commands from the commandfile
        # Substitute all occurances of the "[ -x blah blah blah ]" pattern with the "true" command
        cmd = re.sub( r'\[ -x[^\]]*\]', 'true', cmd )

        # We also don't need the standard output restrictions (as we're not doing comparisons
        cmd = re.sub( r'egrep -vf ../../ignore_list', 'cat', cmd )

    if options.suffix: cmd_line_sh = path.join(workdir, 'command.' + options.suffix + '.sh')
    else : cmd_line_sh = path.join(workdir, 'command.sh')

    #Write the command.sh file
    # writing back so test can be easily re-run by user later...
    with io.open(cmd_line_sh, 'w', encoding="UTF-8", errors='backslashreplace') as f:
        f.write(cmd)

    #if "'" in cmd: raise ValueError("Can't use single quotes in command strings!")
    #print cmd; print

    return cmd_line_sh, workdir

def analyze_integration_test( test, outdir, refdir, results, full_log ):
    """
    Look at the specific integration test, and check for errors
    """

    dir_before = path.join(refdir, test)
    dir_after = path.join(outdir, test)
    # diff returns 0 on no differences, 1 if there are differences

    flags = ["-rq"]
    if Options.fulldiff: flags = ["-r"]
    flags += ['--exclude=command.sh', '--exclude=*.ignore']

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
        print("ok   %s" % test)
        full_log.append("ok   %s\n" % test)
        return 0
    else:
        #runtimes[test] = float('nan')
        print(msg)
        full_log.append( full_log_msg )
        return 1

def analyze_valgrind_test( test, outdir, results, full_log ):
    valgrind_output = []
    dir = path.join(outdir, test)

    # Find all the Valgrind formatted lines in the log files
    def recurse( dirname, outlines ):
        nfiles_with_output = 0
        for fn in os.listdir( dirname ): #Looking for all the files
            fn = os.path.join( dirname, fn )
            if os.path.isdir( fn ):
                nfiles_with_output += recurse( fn, outlines )
                continue
            if (not os.path.isfile( fn )) or fn.endswith("valgrind.out"):
                continue
            has_output = False
            with io.open(fn, encoding="UTF-8", errors='backslashreplace') as f:
                for line in f:
                    if line.startswith("=="):
                        outlines.append(line)
                        has_output = True
            if has_output:
                nfiles_with_output += 1

        return nfiles_with_output

    number_logs = recurse(dir, valgrind_output)

    # Check that the number of log files with valgrind output matches the number we expect from
    # the number of commands we ran (i.e. there isn't a missing log somewhere)
    number_expected = 0
    with open( os.path.join( dir, "command.sh" ) ) as f:
        for line in f:
            # Lines for valgrind runs start in column 0.
            # To turn off counting a line (e.g. because it's surrounded in an if clause)
            # Simply indent it.
            if line.startswith( Options.valgrind_path ):
                number_expected += 1

    # Count the summary error lines, and sum the number of errors
    number_summary = 0
    total_errors = 0
    for line in valgrind_output:
        if line.find("ERROR SUMMARY") == -1:
            continue
        number_summary += 1
        total_errors += int(line.split()[3])

    msg = ''
    # We're fine as long as the number of counted commands matches either the number of summary lines or the number of logfiles.
    if number_summary != number_expected and number_logs != number_expected:
        msg += " --Log file(s) missing (expected %d results, found %d entries in %d log files) -- " % (number_expected, number_summary, number_logs)
        valgrind_output.append('\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
        valgrind_output.append(msg)
        total_errors += 1 # Missing files count as an "error"
    if total_errors != 0:
        msg += " Found " + str( total_errors ) + " Valgrind error(s)."

    # Save the valgrind output specifically
    with open( os.path.join(dir, "valgrind.out" ), 'w' ) as f:
        f.writelines( valgrind_output )

    results[test] = total_errors # I think this is what we should be attaching to the results object

    if len(msg)  == 0:
        print("ok   %s" % test)
        full_log.append( "ok   %s\n" % test )
        return 0
    else:
        print("FAIL %s: %s" % (test, msg))
        full_log.append( "FAIL %s: %s\n" % (test, msg) )

        return 1

def simple_job_running( GenerateJob, queue, outdir, runtimes, options, globalparams ):
    """
    Using the function GenerateJob to generate the job commandlines, run the jobs in queue with the given options.

    GenerateJob signature:
      cmd_line_sh, workdir = GenerateJob(test, outdir, globalparams=None, options=None, host=None)

    """

    def signal_handler(signal_, f):
        print('Ctrl-C pressed... killing child jobs...')
        for nt in Jobs:
            os.killpg(os.getpgid(nt.pid), signal.SIGKILL)

    signal.signal(signal.SIGINT, signal_handler)

    while not queue.empty():
        test = queue.get()
        if test is None: break

        cmd_line_sh, workdir = GenerateJob(test, outdir, globalparams=globalparams, options=options);

        if ( ( cmd_line_sh is None ) or (workdir is None) ) :
            print("No correct command.sh file found for %s.  Skipping." % test)
            queue.task_done() # Mark the task done, for proper counting.
            continue

        def run(nt, times):
            #execute('Running Test %s' % test, 'bash ' + cmd_line_sh)
            test_timeout = Options.timeout
            if nt.test in timeout_factors:
                test_timeout *= timeout_factors[nt.test]
                test_timeout = int(test_timeout)

            extra = 'ulimit -t %s && ' % test_timeout if test_timeout else ''
            res = execute('Running Test %s' % nt.test, '%sbash %s' % (extra, cmd_line_sh), return_=True)
            if res:
                error_string = "*** Test %s did not run!  Check your --mode flag and paths. [%s]\n" % (test, datetime.datetime.now())
                with open(path.join(nt.workdir, ".test_did_not_run.log"), 'w') as f:
                    f.write(error_string)
                print(error_string, end='')
                times[test] = float('nan')

            #execute('Just sleeeping %s...' % test, 'ulimit -t%s && sleep 60 && echo "Done!"' % Options.timeout)
            #execute('Just echo %s...' % test, 'echo "%s Done!"' % test)
            #print 'Not even echo %s... ' % test

        def normal_finish(nt, times):
            queue.task_done()
            percent = (100* (queue.TotalNumberOfTasks-queue.qsize())) / queue.TotalNumberOfTasks
            elapse_time = time.time() - nt.start_time
            print("Finished %-40s in %3i seconds\t [~%4s test (%s%%) started, %4s in queue, %4d running]" % (nt.test, elapse_time, queue.TotalNumberOfTasks-queue.qsize(), percent, queue.qsize(), queue.unfinished_tasks-queue.qsize() ))
            if nt.test not in times:
                times[nt.test] = elapse_time


        def error_finish(nt, times):
            error_string = "*** Test %s did not run!  Check your --mode flag and paths. [%s]\n" % (nt.test, datetime.datetime.now())
            with open(path.join(nt.workdir, ".test_did_not_run.log"), 'w') as f:
                f.write(error_string)
            print(error_string, end='')
            times[nt.test] = float('nan')
            normal_finish(nt, times)

        def timeout_finish(nt, times):
            test_timeout = Options.timeout
            if nt.test in timeout_factors:
                test_timeout *= timeout_factors[nt.test]
            error_string = "*** Test %s exceeded the timeout=%s  and will be killed! [%s]\n" % (nt.test, test_timeout, datetime.datetime.now())
            with open(path.join(nt.workdir, ".test_got_timeout_kill.log"), 'w') as f:
                f.write(error_string)
            print(error_string, end='')
            times[nt.test] = float('inf')
            normal_finish(nt, times)

        if options.jobs > 1:
            test_timeout = Options.timeout
            if test in timeout_factors:
                test_timeout *= timeout_factors[test]
            pid, nt = mFork(times=runtimes, test=test, workdir=workdir, queue=queue, timeout=test_timeout, normal_finish=normal_finish, error_finish=error_finish, timeout_finish=timeout_finish)
            if not pid:  # we are child process
                signal.signal(signal.SIGINT, signal.SIG_DFL)
                run(nt,runtimes)
                sys.exit(0)
        else:
            test_timeout = Options.timeout
            if test in timeout_factors:
                test_timeout *= timeout_factors[test]
            nt = NT(times=runtimes, test=test, workdir=workdir, queue=queue, start_time=time.time(), timeout=test_timeout, normal_finish=normal_finish, error_finish=error_finish, timeout_finish=timeout_finish)
            run(nt,runtimes)
            if nt.timeout and (time.time() - nt.start_time > nt.timeout): nt.timeout_finish(nt,runtimes)
            else: normal_finish(nt,runtimes)

    mWait(all_=True)  # waiting for all jobs to finish before movinf in to next phase

def parallel_job_running(GenerateJob, queue, outdir, runtimes, options, globalparams ):
    # Start worker thread(s)
    for i in range(options.num_procs):
        worker = Worker(queue, outdir, options, times=runtimes, timeout_minutes=options.timeout, GenerateJob=GenerateJob, globalparams=globalparams)
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
          worker = Worker(queue, outdir, options, times=runtimes, host=host, timeout_minutes=options.timeout,GenerateJob=GenerateJob, globalparams=globalparams)
          thread = threading.Thread(target=worker.work)
          #thread.setDaemon(True) # shouldn't be necessary here
          thread.start()

    # Wait for them to finish
    queue.join()

class Worker:
    def __init__(self, queue, outdir, opts, times, host=None, timeout_minutes=0, GenerateJob=generateTestCommandline, globalparams=None):
        self.queue = queue
        self.outdir = outdir
        self.opts = opts
        self.host = host
        self.timeout = timeout_minutes * 60
        self.times = times
        self.GenerateJob = GenerateJob
        self.globalparams = globalparams

    def work(self):
        running=0
        try:
            while True:
                test = self.queue.get_nowait()
                try: # Actually catch exception and ignore it.  Python 2.4 can't use "except" and "finally" together.
                    start = time.time() # initial guess at start time, in case of exception
                    try: # Make sure job is marked done even if we throw an exception
                        cmd_line_sh, workdir = self.GenerateJob(test, self.outdir, options=self.opts, host=self.host, globalparams=self.globalparams)

                        if ( ( cmd_line_sh is None ) or ( workdir is None ) ) :
                            print("No command file found for %s.  Skipping." % test)
                            continue

                        if self.host is None:
                            print("Running  %-40s on localhost ..." % test)
                            proc = subprocess.Popen(["bash",  cmd_line_sh], preexec_fn=os.setpgrp)
                        # Can't use cwd=workdir b/c it modifies *local* dir, not remote dir.
                        else:
                            print("Running  %-40s on %20s ..." % (test, self.host))
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
                                print("*** Test %s exceeded the timeout and will be killed! [%s]\n" % (test, datetime.datetime.now()))
                                self.times[test] = float('inf')
                                #os.kill(proc.pid, signal.SIGTERM)
                                os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
                        if retcode != 0 and retcode is not None:
                            self.times[test] = float('nan')
                            if self.host is None:
                              error_string = "*** Test %s did not run on host %s!  Check your --mode flag and paths. [%s]\n" % (test, 'local_host', datetime.datetime.now())
                            else:
                              error_string = "*** Test %s did not run on host %s!  Check your --mode flag and paths. [%s]\n" % (test, self.host, datetime.datetime.now())
                            print(error_string, end='')

                            # Writing error_string to a file, so integration test should fail for sure
                            with open(path.join(workdir, ".test_did_not_run.log"), 'w') as f:
                                f.write(error_string)

                    finally: # inner try
                        percent = (100* (self.queue.TotalNumberOfTasks-self.queue.qsize())) / self.queue.TotalNumberOfTasks
                        elapse_time = time.time() - start
                        print("Finished %-40s in %3i seconds\t [~%4s test (%s%%) started, %4s in queue, %4d running]" % (test, elapse_time, self.queue.TotalNumberOfTasks-self.queue.qsize(), percent, self.queue.qsize(), self.queue.unfinished_tasks-self.queue.qsize() ))
                        if test not in self.times: self.times[test] = elapse_time
                        self.queue.task_done()

                except Exception as e: # middle try
                    print(e)
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


def local_copytree(src, dst, symlinks=False, accept=lambda srcname, dstname: True):
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
            #NOTE: JAB - symlinks that are hardcoded to random places will fail on shutil and not be able to be copied!
            if symlinks and os.path.islink(srcname):
                linkto = os.readlink(srcname)
                os.symlink(linkto, dstname)
            elif os.path.isdir(srcname):
                local_copytree(srcname, dstname, symlinks, accept)
            else:
                shutil.copy2(srcname, dstname)
            # XXX What about devices, sockets etc.?
        except (IOError, os.error) as why:
            errors.append((srcname, dstname, str(why)))
        # catch the Error from the recursive copytree so that we can
        # continue with other files
        except shutil.Error as err:
            errors.extend(err.args[0])
    try:
        shutil.copystat(src, dst)
    #except WindowsError:
        # can't copy file access times on Windows
        pass
    except OSError as why:
        errors.extend((src, dst, str(why)))
    if errors:
        raise shutil.Error(errors)


def match_patterns(search_string, patterns):
    """
    Uses RE to match multiple patterns.  Returns boolean of success

    :param search_string: str
    :param patterns: [str]
    :rtype: boolean

    """
    match = False

    for pattern in patterns:
        if re.search(pattern, search_string):
            match = True
            break
    return match

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
