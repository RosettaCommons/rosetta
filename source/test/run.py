#!/usr/bin/env python
# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
#
## @file   test/run.py
## @brief  Script to run unit tests in mini
## @author Sergey Lyskov


import sys
if not hasattr(sys, "version_info") or sys.version_info < (2,4):
    raise ValueError("Script requires Python 2.4 or higher!")

import os, commands, re, subprocess, time
from os import path
from optparse import OptionParser

import yaml

UnitTestExecutable = ["protocols.test", "core.test", "basic.test", "ObjexxFCL.test", "numeric.test", "utility.test", "apps.test", "devel.test"]
#UnitTestExecutable = ["numeric.test", "utility.test"]

# Output info class, simple wrapper for named collection of fields.
class OI:
    def __init__(self, **entries): self.__dict__.update(entries)


class Tester:
    def __init__(self):
        #self.db_path = db_path
        self.systemLog = ""  # System log - we store all information here
        self.results = {}
        self.jobs = []  # list of spawned process pid's
        self.platform = None
        self.testpath = None


    # print and log given mesage, return message
    def log(self, s):
        self.systemLog += s
        print s
        return s


    def mfork(self):
        ''' Check if number of child process is below Options.jobs. And if it is - fork the new pocees and return its pid.
        '''
        while len(self.jobs) >= Options.jobs :
            for p in self.jobs[:] :
                r = os.waitpid(p, os.WNOHANG)
                if r == (p, 0):  # process have ended without error
                    self.jobs.remove(p)
                elif r[0] == p :  # process ended but with error, special case we will have to wait for all process to terminate and call system exit.
                    for p in jobs: os.waitpid(p, 0)
                    print 'Some of the unit test suite terminate abnormally!'
                    sys.exit(1)

            if len(self.jobs) >= Options.jobs: time.sleep(.5)
        pid = os.fork()
        if pid: self.jobs.append(pid) # We are parent!
        return pid



    # Try to identity plaform by using scons compiliation feature.
    def getPlatformID(self):
        if Options.CMake:
            self.testpath = Options.testpath
            self.log( "Skipping platform identification. Using explicit directory instead: " + self.testpath )
            return
        self.log( "Identifying platform...\n")
        cmd_str = "./scons.py unit_test_platform_only log=platform"
        if Options.extras:
            cmd_str += " extras=%s" % (Options.extras)
        if Options.mode:
            cmd_str += " mode=%s" % (Options.mode)
        if Options.compiler:
            cmd_str += " cxx=%s" % (Options.compiler)

        pl = commands.getoutput(cmd_str)
        lines = pl.split('\n')
        for s in lines:
            if  len( s.split() ) > 1 and s.split()[0] == 'Platform:':
                platform = s.split()[1]
                self.log( "Platform found: " + platform )
                self.platform = platform
                self.testpath = "build/test/" + platform
                return
        sys.exit("run.py is about to crash because it could not use SCons to detect your platform.  The most likely reason for this is that you are running it from the wrong directory - it must be run from the rosetta_source directory, not the rosetta_source/test directory, even though it lives in the latter.")
        return "PlatformWasNotFound!!!"  # <-- That should not reall happend.


    # extract information regarding how good unit tests run.
    def extractInfo(self, output):
        # default init, in case we can't extract information
        oi = OI(testCount=0, testFailed=0, failedTestsList=[])

        # extracting number of test
        s = re.search(r'Running (\d\d*) test', output)
        if s:
            g = s.groups()
            oi.testCount = int( g[0] )

        # extracting number of test failed
        s = re.search(r'Failed (\d\d*) of \d\d* test', output)
        if s:
            g = s.groups()
            oi.testFailed = int( g[0] )

            # extracting names of the failed tests
            s = re.findall(r'CXXTEST_ERROR: (\w*) Failed!', output)
            for t in s: oi.failedTestsList.append(t)

            #s = re.search(r'CXXTEST_ERROR: (\w*) Failed!', output)
            #if s:
            #    g = s.groups()
            #    for t in g: oi.failedTestsList.append(t)

        else: # all test pussed?
            s = re.search(r'All tests passed!', output)
            if not s: # Something wrong then, count all tests as failed.
                oi.testFailed = oi.testCount
                oi.failedTestsList = ['All']

        #print oi.__dict__
        return oi


    def runOneLibUnitTests(self, lib, yaml_file, log_file):
        if Options.one  and  ( (lib, Options.one.split(':')[0]) not in self.all_test_suites_by_lib[lib] ): return

        #self.unitTestLog += self.log("-------- %s --------\n" % E)
        path = "cd " + self.testpath + " && "
        mute, unmute  = ' ', ' '
        if Options.mute:   mute   = ' -mute '   + ' '.join(Options.mute)
        if Options.unmute: unmute = ' -unmute ' + ' '.join(Options.unmute)
        exe = "./" + lib + ' ' + Options.one + ' --database ' + Options.database + mute + unmute
        #if self.db_path: exe += " " + self.db_path
        print "Command line:: %s%s" % (path, exe)

        if os.path.isfile(yaml_file): os.remove(yaml_file)

        output = "Running %s unit tests..." % lib

        timelimit = sys.executable + ' ' +os.path.abspath('test/timelimit.py') + ' 10 '
        #print "timelimit:", timelimit

        command_line = path + ' ' + timelimit + exe + " 1>&2"
        #print 'Executing:', command_line


        f = subprocess.Popen(command_line, bufsize=0, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).stderr
        for line in f:
            print line,
            output += line
            sys.stdout.flush()
        f.close()

        #else:
        #    print 'Test %s is not in %s, skipping...' % (Options.one, lib)
        #    print 'tests:', self.all_test_suites_by_lib[lib]

        # saving log to a file...
        f = open(log_file, 'w');  f.write(output);  f.close()


    def runOneSuite(self, lib, suite):
        path = "cd " + self.testpath + " && "
        mute, unmute  = ' ', ' '
        if Options.mute:   mute   = ' -mute '   + ' '.join(Options.mute)
        if Options.unmute: unmute = ' -unmute ' + ' '.join(Options.unmute)
        exe = "./" + lib + ' ' + suite + ' --database ' + Options.database + mute + unmute
        #if self.db_path: exe += " " + self.db_path
        #print "Paths:", path, 'command line:', exe

        log_file = self.testpath + '/' + lib + '.' + suite + '.log'
        yaml_file = self.testpath + '/' + lib + '.'+ suite + '.yaml'

        if os.path.isfile(yaml_file): os.remove(yaml_file)
        if os.path.isfile(log_file): os.remove(log_file)

        #output = ''
        output = "Running %s:%s unit tests..." % (lib, suite)
        #print output

        timelimit = sys.executable + ' ' + os.path.abspath('test/timelimit.py') + ' 30 '
        #print "timelimit:", timelimit

        command_line = path + ' ' + timelimit + exe + " 1>&2"
        #print 'Executing:', command_line
        f = subprocess.Popen(command_line, bufsize=0, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).stderr
        for line in f:
            if line != 'All tests passed!\n':
                print line,
                output += line
                sys.stdout.flush()
        f.close()
        #output += '\n\n'; print '\n'

        # saving log to a file...
        f = open(log_file, 'w');  f.write(output);  f.close()



    # Run unit test.
    def runUnitTests(self):
        self.getPlatformID()
        #self.log( "Run unit tests...\n")
        self.unitTestLog = "================================ UnitTest Results ================================\n"

        logs_yamls = {}

        # Getting list of test suite's
        self.all_test_suites = []
        self.all_test_suites_by_lib = {}

        self.all_tests = []

        for lib in UnitTestExecutable:
            tests = set()
            for test in commands.getoutput(self.testpath + '/' + lib + ' _ListAllTests_').split():
                tests.add( (lib, test.split(':')[0]) )
                self.all_tests.append(test)

            self.all_test_suites.extend( tests )
            self.all_test_suites_by_lib[lib] = tests

        #print 'Suites:', self.all_test_suites
        #print 'Tests:', self.all_tests

        if Options.one:  # or Options.jobs < 5:
            if Options.one not in [s for (l,s) in self.all_test_suites] + self.all_tests:
                print 'Test suite %s not found!' % Options.one
                sys.exit(1)

            for lib in UnitTestExecutable:
                log_file = self.testpath + '/' + lib + '.log'
                yaml_file = self.testpath + '/' + lib + '.yaml'

                logs_yamls[lib] = (log_file, yaml_file)

                if Options.jobs > 1:
                    pid = self.mfork()
                    if not pid:  # we are child process
                        self.runOneLibUnitTests(lib, yaml_file, log_file)
                        sys.exit(0)

                else:
                    self.runOneLibUnitTests(lib, yaml_file, log_file)

            for p in self.jobs: os.waitpid(p, 0)  # waiting for all child process to termintate...

            # extracting and aggregation results, but only if running all tests
            if not Options.one:
                all_results_yaml_file = '.unit_test_results.yaml'
                uf = file(all_results_yaml_file, 'w')
                for lib in logs_yamls:
                    (log_file, yaml_file) = logs_yamls[lib]
                    info = self.extractInfo( file(log_file).read() )
                    info.name = lib
                    self.results[lib] = info

                    if not os.path.isfile(yaml_file):
                        print "Unable to read yaml file with test results %s - unit test run aborted!" % yaml_file
                        uf.close()
                        os.remove(all_results_yaml_file)
                        sys.exit(1)

                    # generating summary yaml file for all tests
                    uf.write(lib + ':\n')
                    f = file(yaml_file)
                    for line in f: uf.write("    "+line)
                    f.close()
                    uf.write('\n')

                uf.close()

            #self.log( "Run unit tests... Done.\n")

        else: # running Unit test on multiple CPU's, new style, fully parallel
            for lib, suite in self.all_test_suites:
                pid = self.mfork()
                if not pid:  # we are child process
                    self.runOneSuite(lib, suite)
                    sys.exit(0)

            for p in self.jobs: os.waitpid(p, 0)  # waiting for all child process to termintate...

            # Now all tests should be finished, all we have to do is to create a log file and aggegated yaml file to emulate single CPU out run
            all_yaml = {}
            for lib in UnitTestExecutable:
                log_file = self.testpath + '/' + lib + '.log'
                yaml_file = self.testpath + '/' + lib + '.yaml'

                logs_yamls[lib] = (log_file, yaml_file)

                log_file_h = file(log_file, 'w')
                yaml_file_h = file(yaml_file, 'w')
                yaml_data = {}
                for l, suite in self.all_test_suites_by_lib[lib]:
                    log_file_h.write( file(self.testpath + '/' + lib + '.' + suite + '.log').read() )
                    data = yaml.load( file(self.testpath + '/' + lib + '.' + suite + '.yaml').read() )
                    for k in data:
                        if k in yaml_data: yaml_data[k] = list( set(yaml_data[k] + data[k]) )
                        else: yaml_data[k] = data[k]


                yaml_file_h.write( yaml.dump(yaml_data) )

                log_file_h.close()
                yaml_file_h.close()

                #if 'ALL_TESTS' in yaml_data:
                if yaml_data:
                    self.results[lib] = OI(testCount=len(yaml_data['ALL_TESTS']), testFailed=len(yaml_data['FAILED_TESTS']), failedTestsList=yaml_data['FAILED_TESTS'], name=lib)
                    all_yaml[lib] = yaml_data


                logs_yamls[lib] = (log_file, yaml_file)

            #print 'All_yaml:', all_yaml
            f = file('.unit_test_results.yaml', 'w');  f.write( yaml.dump(all_yaml) );  f.close()


        '''
        error = False
        for p in self.jobs:
            r = os.waitpid(p, 0)  # waiting for all child process to termintate...
            if r[0] == p :  # process ended but with error, special case we will have to wait for all process to terminate and call system exit.
                error = True

        if error:
            print 'Some of the build scripts return an error, PyRosetta build failed!'
            sys.exit(1)
        '''



    def printSummary(self):
        total = 0
        failed = 0
        failedTestsList = []

        for t in self.results:
            total += self.results[t].testCount
            failed += self.results[t].testFailed
            failedTestsList.extend( [ self.results[t].name + ": " + r  for r in self.results[t].failedTestsList] )
            #print "--> ", self.results[t].failedTestsList

        print "-------- Unit test summary --------"
        print "Total number of tests:", total
        print "  number tests passed:", total-failed
        print "  number tests failed:", failed
        if failedTestsList:
            print "  failed tests:"
            for t in failedTestsList:
                print "   ", t
        print "Success rate: %s%%" % ((total-failed)*100/total)
        print "---------- End of Unit test summary"


def main(args):
    ''' Script to run Unit test in  mini.
    '''
    parser = OptionParser(usage="usage: %prog [OPTIONS] [TESTS]")
    parser.set_description(main.__doc__)

    parser.add_option("-d", '--database',
      default="", #processed below
      action="store",
      help="Path to Rosetta database. (default: $ROSETTA3_DB, ../rosetta_database, ~/rosetta_database)",
    )

    parser.add_option("--extras",
      action="store",
      help="Extras option passed to scons to help identify platform"
    )

    parser.add_option("--mode",
      action="store",
      help="mode option passed to scons to help identify platform"
    )


    parser.add_option("-1", "--one",
      default='', action="store",
      help="Run just one unit test or one test suite. To run one unit test-suite just specify it name, ie: '--one PDB_IO'. \
            To run one unit test specify name of test suit, colon, then name of the test it self, like this: '--one PDB_IO:test_pdb_io'.  \
            It is HIGHLY RECOMMENDED to run one unit test suit instead of one unit test! (when running one test all suite have to be initialized, while when running one suit this is not the case).",
    )

    parser.add_option("-j", "--jobs",
      default=1,
      type="int",
      help="Number of processors to use when running testss (default: 1)",
    )

    parser.add_option('--mute',
      default=[],
      action="append",
      help="Mute specified tracer channels.",
    )

    parser.add_option('--unmute',
      default=[],
      action="append",
      help="UnMute specified tracer channels.",
    )

    parser.add_option("-c", '--compiler',
      default='gcc',
      action="store",
      help="Name of the compiler used.",
    )

    parser.add_option("-C", '--CMake',
      default=False,
      action="store_true",
      help="Was CMake used to build the unit tests? (i.e. obey the --testpath option).",
    )

    parser.add_option("-T", '--testpath',
      default="cmake/build_unit/",
      action="store",
      help="The relative directory where the unit tests were built into. (default: 'cmake/build_unit/')",
    )

    (options, args) = parser.parse_args(args=args[1:])

    unknown_args = [a.strip() for a in args if a.strip()]
    if unknown_args:
        raise ValueError("Unknown options: %s", unknown_args)

    if options.database == parser.get_default_values().database:
        if path.isdir( path.join( "..", "database") ):
            options.database = path.abspath(path.join( "..", "database"))

        elif os.environ.get('ROSETTA3_DB') is not None and \
                path.isdir(os.environ.get('ROSETTA3_DB')):
            options.database = os.environ.get('ROSETTA3_DB')

        elif path.isdir( path.join( path.expanduser("~"), "rosetta_database") ):
            options.database = path.join( path.expanduser("~"), "rosetta_database")

        else:
            raise ValueError("Can't find database at %s; please set $ROSETTA3_DB or use -d" % options.database)

    else: options.database = path.abspath( options.database )

    global Options;  Options = options


    #db_path = None
    #if len(args) <= 1: print "Warning: no database path was specified!"
    #else: db_path = " ".join(args[1:])

    T = Tester()
    T.runUnitTests()
    if not Options.one: T.printSummary()

    print "Done!"



if __name__ == "__main__": main(sys.argv)
