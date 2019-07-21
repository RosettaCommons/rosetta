#!/usr/bin/env python
# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.
#
## @file   test/run.py
## @brief  Script to run unit tests in mini
## @author Sergey Lyskov
from __future__ import print_function

import sys
if not hasattr(sys, "version_info") or sys.version_info < (2,4):
    raise ValueError("Script requires Python 2.4 or higher!")

import os, re, subprocess, time, os.path
from os import path
from optparse import OptionParser

import json, codecs


#The factor by which to multiply the single test timeout value for a suite
SUITE_FACTOR = 3

UnitTestExecutable = ["protocols.test", "core.test", "basic.test", "ObjexxFCL.test", "numeric.test", "utility.test", "apps.test", "devel.test"]
#UnitTestExecutable = ["numeric.test", "utility.test"]

# white list of unit tests for which we allow longer timeouts due to historical reasons
# map lib:suite -> timeout-factor
# Default for tests not listed is "1"
_timeout_factors_ = {
    "core:SingleNCAARotamerLibraryTests_equivalency_no_voronoi" : 8,
    "core:SingleNCAARotamerLibraryTests_equivalency_voronoi" : 8,
    "protocols:InterfaceFeaturesTests" : 4,
    "protocols:BridgeChainsMoverTests" : 20,
    "core:AtomTypeDatabaseIOTests" : 20,
    "protocols:PCSEnergyTests" : 20,
    "core:ValidateBetaBinariesTests" : 20,
    "core:ValidateDefaultBinariesTests": 20,
    "core:ValidateTalarisBinariesTests": 20,
    "core:CyclicGeometry_beta_peptoid_TwoChainTests" : 5,
    "core:CyclicGeometry_peptoid_TwoChainTests" : 5,
    "core:CyclicGeometry_betanov16_TwoChainTests" : 2,
    "core:CyclicGeometry_nmethyl_betanov15_TwoChainTests" : 2,
    "core:AtomicDepthTests" : 2,
    "protocols:MembraneUtil" : 2,
    "protocols:AddMembraneMoverTest" : 2,
}

def to_unicode(b):
    ''' Conver bytes to string and handle the errors. If argument is already in string - do nothing
    '''
    if not hasattr(sys, "version_info") or sys.version_info < (3, 0): return b if type(b) == unicode else unicode(b, 'utf-8', errors='backslashreplace')
    else: return b if type(b) == str else str(b, 'utf-8', errors='backslashreplace')


def execute(message, command_line, return_='status', until_successes=False, terminate_on_failure=True, silent=False, silence_output=False, silence_output_on_errors=False, add_message_and_command_line_to_output=False):
    if not silent: print(message);  print(command_line); sys.stdout.flush();
    while True:

        p = subprocess.Popen(command_line, bufsize=0, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, errors = p.communicate()
        output = output + errors
        output = output.decode(encoding='utf-8', errors='backslashreplace')
        exit_code = p.returncode

        #exit_code, output = subprocess.getstatusoutput(command_line)

        if (exit_code  and  not silence_output_on_errors) or  not (silent or silence_output): print(output); sys.stdout.flush();

        if exit_code and until_successes: pass  # Thats right - redability COUNT!
        else: break

        print( "Error while executing {}: {}\n".format(message, output) )
        print("Sleeping 60s... then I will retry...")
        sys.stdout.flush();
        time.sleep(60)

    if add_message_and_command_line_to_output: output = message + '\nCommand line: ' + command_line + '\n' + output

    if return_ == 'tuple': return(exit_code, output)

    if exit_code and terminate_on_failure:
        print("\nEncounter error while executing: " + command_line)
        if return_==True: return True
        else:
            print('\nEncounter error while executing: ' + command_line + '\n' + output);
            raise BenchmarkError('\nEncounter error while executing: ' + command_line + '\n' + output)

    if return_ == 'output': return output
    else: return exit_code



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
        print(s)
        return s


    def mfork(self):
        ''' Check if number of child process is below Options.jobs. And if it is - fork the new process and return its pid.
        '''
        while len(self.jobs) >= Options.jobs :
            for p in self.jobs[:] :
                r = os.waitpid(p, os.WNOHANG)
                if r == (p, 0):  # process have ended without error
                    self.jobs.remove(p)
                elif r[0] == p :  # process ended but with error, special case we will have to wait for all process to terminate and call system exit.
                    for p in self.jobs:
                        try:
                            os.waitpid(p, 0)
                        except OSError: pass
                    print('Some of the unit test suite terminate abnormally!')
                    sys.exit(1)

            if len(self.jobs) >= Options.jobs: time.sleep(.5)
        pid = os.fork()
        if pid: self.jobs.append(pid) # We are parent!
        return pid

    def wait_all(self):
        ''' Wait for all jobs to finish. '''
        oldjobs, self.jobs = self.jobs, [] # Swap to be safe about adding new jobs while waiting.
        for p in oldjobs: os.waitpid(p, 0) # waiting for all child process to terminate...

    # Try to identity plaform by using scons compiliation feature.
    def getPlatformID(self):
        if Options.CMake:
            self.testpath = Options.testpath
            self.log( "Skipping platform identification. Using explicit directory instead: " + self.testpath )
            return
        self.log( "Identifying platform...\n")
        # We use sys.executable here, so we use the same Python as the test script (if it was explicitly specified).
        cmd_str = sys.executable + " ./scons.py unit_test_platform_only log=platform"
        if Options.extras:
            cmd_str += " extras=%s" % (Options.extras)
        if Options.mode:
            cmd_str += " mode=%s" % (Options.mode)
        if Options.compiler:
            cmd_str += " cxx=%s" % (Options.compiler)
        if Options.cxx_ver:
            cmd_str += " cxx_ver=%s" % (Options.cxx_ver)

        pl = subprocess.check_output(cmd_str, shell=True).decode('utf-8', errors="backslashreplace")
        lines = pl.split('\n')
        for s in lines:
            if  len( s.split() ) > 1 and s.split()[0] == 'Platform:':
                platform = s.split()[1]
                self.log( "Platform found: " + platform )
                self.platform = platform
                self.testpath = "build/test/" + platform
                return
        sys.exit("run.py is about to crash because it could not use SCons to detect your platform.  The most likely reason for this is that you are running it from the wrong directory - it must be run from the rosetta_source directory, not the rosetta_source/test directory, even though it lives in the latter.\nScons output for command line {0} is:\n{1}\n".format(cmd_str, pl))
        return "PlatformWasNotFound!!!"  # <-- That should not really happend.


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


    def runOneLibUnitTests(self, lib, test_name, yjson_file, log_file):
        #self.unitTestLog += self.log("-------- %s --------\n" % E)
        path = "cd " + self.testpath + " && "
        exe = self.buildCommandLine(lib, test_name)
        print("Command line: %s%s" % (path, exe))

        if os.path.isfile(yjson_file): os.remove(yjson_file)

        output = [ "Running %s unit tests...\n" % lib ]

        if( Options.timeout > 0 ):
            my_timeout = Options.timeout * _timeout_factors_.get(lib[:-len('.test')]+':'+test_name, 1)
            if ':' not in test_name:
                my_timeout *= SUITE_FACTOR
            timelimit = sys.executable + ' ' +os.path.abspath('test/timelimit.py') + ' ' + str(my_timeout) + ' '
        else:
            timelimit = ""

        #print( "timelimit:", timelimit )

        command_line = path + ' ' + timelimit + exe + " 1>&2"
        #print( ">>>>>>", 'Executing:', command_line)

        f = subprocess.Popen(command_line, bufsize=0, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).stderr
        if Options.valgrind:
            self.parseOneLibValgrind(f, output, yjson_file)
        else:
            self.parseOneLib(f, output)

        f.close()

        output = [ to_unicode(s) for s in output ]
        with codecs.open(log_file, 'w', encoding='utf-8', errors='backslashreplace') as f: f.write(''.join(output))

    def parseOneLib(self, f, output):
        '''Parse output for OneLib run for regular unit tests.'''
        for line in f:
            l = line.decode('utf-8', errors="backslashreplace")
            print(l, end='')
            output.append(l)
            sys.stdout.flush()

    def parseOneLibValgrind(self, f, output, yjson_file):
        valgrind_errors = 0
        for line in f:
            if line.startswith("=="):
                print(line.strip())
                if line.find("ERROR SUMMARY") != -1:
                    valgrind_errors += int( line.split()[3] )
                output.append(line)
                sys.stdout.flush()

        #The yjson file might not be the same name - e.g. if we're running a single suite
        if valgrind_errors > 0 and os.path.isfile(yjson_file):
            data = json.loads(open(yjson_file).read())
            failures = []
            for test in data["ALL_TESTS"]:
                if any( [ test.startswith( testname ) for testname in Options.test_list ] ):
                    failures.append( test )
            data["FAILED_TESTS"] = failures
            yjson = open( yjson_file, "w" )
            yjson.write( json.dumps(data) )
            yjson.close()

    def outfileOneSuite(self, lib, suite, extension):
        filename = self.testpath + '/' + lib + '.' + suite + '.' + extension
        filename = filename.replace(':','__')
        return filename

    def runOneSuite(self, lib, suite):
        assert ':' not in suite
        path = "cd " + self.testpath + " && "
        exe = self.buildCommandLine(lib, suite)

        log_file = self.outfileOneSuite( lib, suite, 'log' )
        yjson_file = self.outfileOneSuite( lib, suite, 'json' )

        if os.path.isfile(yjson_file): os.remove(yjson_file)
        if os.path.isfile(log_file): os.remove(log_file)

        #output = [ ]
        output = [ "Running %s:%s unit tests...\n" % (lib, suite) ]
        #print output

        if Options.timeout > 0:
            timelimit = sys.executable + ' ' + os.path.abspath('test/timelimit.py') + ' ' + str( Options.timeout * SUITE_FACTOR * _timeout_factors_.get(lib[:-len('.test')] + ':' + suite, 1) ) + ' '
        else:
            timelimit = ""
        #print "timelimit:", timelimit

        command_line = path + ' ' + timelimit + exe + " 1>&2"
        #print 'Executing:', command_line
        f = subprocess.Popen(command_line, bufsize=0, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).stderr
        if Options.valgrind:
            self.parseOneSuiteValgrind(f, output, lib, suite, yjson_file)
        else:
            self.parseOneSuite(f, output)
        f.close()

        output = [ to_unicode(s) for s in output ]
        with codecs.open(log_file, 'w', encoding='utf-8', errors='backslashreplace') as f: f.write(''.join(output))




    def parseOneSuite(self, f, output):
        for line in f:
            l = line.decode('utf-8', errors="backslashreplace")
            if l != 'All tests passed!\n':
                print(l, end='')
                output.append(l)
                sys.stdout.flush()

    def parseOneSuiteValgrind(self, f, output, lib, suite, yjson_file):
        valgrind_errors = 0
        for line in f:
            if line != 'All tests passed!\n':
                if line.find("ERROR SUMMARY") != -1:
                    valgrind_errors += int( line.split()[3] )
                    if valgrind_errors != 0:
                        print(lib, suite, line.strip())
                output.append(line)
                sys.stdout.flush()

        # Running suites, we really can't determine which of the sub-tests caused the Valgrind error - so mark them all failed.
        if valgrind_errors > 0:
            data = json.loads(open(yjson_file).read())
            failures = []
            for test in data["ALL_TESTS"]:
                if test.startswith( suite ):
                    failures.append( test )
            data["FAILED_TESTS"] = failures
            yjson = open( yjson_file, "w" )
            yjson.write( json.dumps(data) )
            yjson.close()

    # Obtain the list of valid tests.
    # Will initialize self.libs_for_test; self.all_tests and self.all_suites
    def obtainTestList(self):
        # Contains both tests (with colon) and suite names
        # Is a one->many mapping, just in case the same suite is named in multiple libraries
        self.libs_for_test = {} # dictionary of strings -> sets()

        self.all_tests = []
        self.all_suites = set()

        # To speed up running, launch library finding across multiple jobs
        for lib in UnitTestExecutable:
            testlist_output = self.testpath + '/' + lib + '.testlist_output'
            if os.path.isfile(testlist_output): os.remove(testlist_output)
            pid = self.mfork()
            if not pid:  # we are child process
                exit_code, command_output = execute('Getting list of test for {lib}...'.format(**locals()), self.testpath + '/' + lib + ' _ListAllTests_', return_='tuple', silent=True)

                if exit_code:
                    # The output doesn't contain a test list, only an error message
                    print('Could not acquire list of available tests from library {lib}, request terminated with error:\n{command_output}\nTerminating...'.format(**locals()))
                    sys.exit(1)
                else:
                    with open(testlist_output, 'w') as f:
                        f.write(command_output)
                    sys.exit(0)

        self.wait_all() # waiting for all child process to terminate...

        # Now parse
        for lib in UnitTestExecutable:
            testlist_output = self.testpath + '/' + lib + '.testlist_output'
            if not os.path.exists( testlist_output ): continue
            with open(testlist_output) as f:
                command_output = f.read();
            os.remove(testlist_output)

            for test in command_output.split():
                suite = test.split(':')[0]

                self.all_tests.append(test)
                self.all_suites.add(suite)
                self.libs_for_test.setdefault(test, set()).add( lib )
                self.libs_for_test.setdefault(suite, set()).add( lib )

        self.all_tests.sort( key = lambda x: x.lower() )
        self.all_suites = sorted( self.all_suites );

    # Run unit test.
    def runUnitTests(self):
        self.getPlatformID()
        #self.log( "Run unit tests...\n")
        self.unitTestLog = "================================ UnitTest Results ================================\n"

        logs_yjsons = {}

        # Fill out the test info self variables.
        self.obtainTestList()

        # Check tests, and default to running all suites if none specified
        for test_name in Options.test_list:
            if test_name not in self.libs_for_test:
                print('\nTest %s not found!' % test_name)
                print( "Available tests are: {tests}".format(tests = ' '.join(self.all_tests) ) )
                sys.exit(1)

        if len( Options.test_list ) == 0:
            Options.test_list = self.all_suites

        # Actually launch the tests
        if Options.one:
            for test_name in Options.test_list:
                for lib in self.libs_for_test.get(test_name,set()):

                    log_file = self.testpath + '/' + lib + '.log'
                    yjson_file = self.testpath + '/' + lib + '.json'

                    logs_yjsons[lib] = (log_file, yjson_file)

                    pid = self.mfork()
                    if not pid:
                        self.runOneLibUnitTests(lib, test_name, yjson_file, log_file)
                        sys.exit(0)

            self.wait_all()

        else: # running Unit test on multiple CPU's, new style, fully parallel
            for suite in Options.test_list:
                for lib in self.libs_for_test.get( suite, set() ):
                    pid = self.mfork()
                    if not pid:  # we are child process
                        if Options.debug: print('DEBUG MODE ENABLED: Skipping tests run: ' + lib + ' ' + suite)
                        else: self.runOneSuite(lib, suite)
                        sys.exit(0)

            self.wait_all() # waiting for all child process to termintate...

            # Now all tests should be finished, all we have to do is to create a log file and aggegated JSON file to emulate single CPU out run
            all_yjson = {}
            all_json = dict(tests={}, summary=dict(total=0, failed=0), config=dict(test_path=self.testpath), runtime={})
            _failed_ = 'failed'
            _finished_ = 'passed'

            for lib in UnitTestExecutable:
                ntests_for_lib = 0
                log_file = self.testpath + '/' + lib + '.log'
                yjson_file = self.testpath + '/' + lib + '.json'

                logs_yjsons[lib] = (log_file, yjson_file)

                suites_for_lib = [ suite for suite in Options.test_list if lib in self.libs_for_test.get( suite, set() ) ]

                yjson_data = {}
                with open(log_file, 'w') as log_file_h:
                    for suite in suites_for_lib:
                        suite_log = open( self.outfileOneSuite( lib, suite, 'log') ).read()
                        log_file_h.write( suite_log )

                        yjson_suite_file = self.outfileOneSuite( lib, suite, 'json')
                        print('trying: ', yjson_suite_file)
                        with open(yjson_suite_file) as f:
                            data = json.loads( f.read() )

                        for k in data:
                            if k in yjson_data:
                                if type(data[k]) == list:
                                    yjson_data[k] = list( set(yjson_data[k] + data[k]) )
                                else: yjson_data[k].update(data[k])

                            else: yjson_data[k] = data[k]

                        def json_key_from_test(test): return lib[:-len('.test')] + ':' + test.replace(':', ':')  # we might want different separator then ':' in the future

                        for test in data['ALL_TESTS']:
                            if not test.startswith(suite+':'): continue # Only count tests for this suite
                            ntests_for_lib += 1
                            json_key = json_key_from_test(test)
                            if json_key not in all_json['tests']: all_json['tests'][json_key] = dict(state=_finished_, log='')

                        for test in data['FAILED_TESTS']:
                            json_key = json_key_from_test(test)
                            all_json['tests'][json_key] = dict(state=_failed_, log=suite_log)

                if 'runtime' in yjson_data: all_json['runtime'].update( yjson_data['runtime'] )

                with open(yjson_file, 'w') as yjson_file_h:
                    yjson_file_h.write( json.dumps(yjson_data) )

                #if 'ALL_TESTS' in yjson_data:
                if yjson_data:
                    self.results[lib] = OI(testCount=ntests_for_lib, testFailed=len(yjson_data['FAILED_TESTS']), failedTestsList=yjson_data['FAILED_TESTS'], name=lib)
                    all_yjson[lib] = yjson_data

                logs_yjsons[lib] = (log_file, yjson_file)

            #print 'All_yjson:', all_yjson
            # with open('.unit_test_results.json', 'w') as f:
            #     f.write( json.dumps(all_yjson) )

            for k in list( all_json['runtime'].keys() ): all_json['runtime'][ k.replace('./', '').replace('.test', '') ] = all_json['runtime'].pop(k)

            for t in self.results:
                all_json['summary']['total']   += self.results[t].testCount
                all_json['summary']['failed'] += self.results[t].testFailed

            all_json['summary']['failed_tests'] = [ t for t in all_json['tests'] if all_json['tests'][t]['state'] == _failed_ ]

            with open('.unit_test_results.json', 'w') as f:
                json.dump(all_json, f, sort_keys=True, indent=2)


    def printSummary(self):
        total = 0
        failed = 0
        failedTestsList = []

        for t in self.results:
            total += self.results[t].testCount
            failed += self.results[t].testFailed
            failedTestsList.extend( [ self.results[t].name + ": " + r  for r in self.results[t].failedTestsList] )
            #print "--> ", self.results[t].failedTestsList

        print("-------- Unit test summary --------")
        print("Total number of tests:", total)
        print("  number tests passed:", total-failed)
        print("  number tests failed:", failed)
        if failedTestsList:
            print("  failed tests:")
            for t in failedTestsList:
                print("   ", t)
        if total == 0:
            print("Success rate: 0%")
        else:
            print("Success rate: %s%%" % ((total-failed)*100/total))
        print("---------- End of Unit test summary")


    def buildCommandLine(self, lib, test):
        if Options.valgrind:
            preamble = [ Options.valgrind_path ]
            if Options.trackorigins:
                preamble.append( "--track-origins=yes" )
            if Options.leakcheck:
                preamble.append( "--leak-check=full" )
        else:
            preamble = []

        return " ".join( preamble + ["./" + lib, test] +  self.genericTestFlags() )

    def genericTestFlags(self):
        """Generate generic command line flags (database/tracer/etc...) for tests."""

        flags = ["--database", Options.database]

        if Options.mute:
            flags.extend(["-mute"] +  Options.mute)

        if Options.unmute:
            flags.extend(["-unmute"] +  Options.unmute)

        if Options.levels:
            flags.extend(["-out:levels"] + Options.levels)

        if Options.level:
            flags.extend(["-out:level"] + Options.level)

        #No personal flag configs
        flags.append("-no_fconfig")

        return flags

def parse_valgrind_options(options, option_parser):
    if( options.timeout == option_parser.get_default_values().timeout ):
        print("Default value for timeout used - turning off timeout for valgrind.") # Valgrind runs take a long time.
        options.timeout = 0
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

def main(args):
    ''' Script to run Unit tests in Rosetta.
    By default the script will run all unit tests. You can run specific tests by listing them on the command line.
    It is HIGHLY RECOMMENDED to run a unit test suite (no colon) instead of one unit test!
    '''
    parser = OptionParser(usage="usage: %prog [OPTIONS] [TESTS]")
    parser.set_description(main.__doc__)

    parser.add_option("-d", '--database',
      default="", #processed below
      action="store",
      help="Path to Rosetta database. (default: $ROSETTA3_DB, ../database, ~/database, ../rosetta_database, ~/rosetta_database, )",
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
      default=False, action="store_true",
      help="Run using the older 'one-by-one' method. \
            This will automatically be enabled if you specified a subtest (with a colon) instead of a full suite."
    )

    parser.add_option("-j", "--jobs",
      default=1,
      type="int",
      help="Number of processors to use when running tests (default: 1)",
    )

    parser.add_option('--mute',
      default=[],
      action="append",
      help="Mute specified tracer channels. (Note: all channels will be muted by default when more than one test is run.)",
    )

    parser.add_option('--unmute',
      default=[],
      action="append",
      help="UnMute specified tracer channels. (Useful when debugging and running more than one test.)",
    )

    parser.add_option('--levels',
      default=[],
      action="append",
      help="Tracer out:levels configuration.",
    )

    parser.add_option('--level',
      default=[],
      action="append",
      help="Tracer out:level configuration (For example, 400/500)",
    )

    parser.add_option("-c", '--compiler',
      default=None,
      action="store",
      help="Name of the compiler used.",
    )

    parser.add_option('--cxx_ver',
      default=None,
      action="store",
      help="Version of the compiler used.",
    )

    parser.add_option("-C", '--CMake',
      default=False,
      action="store_true",
      help="Was CMake used to build the unit tests? (i.e. obey the --testpath option).",
    )

    parser.add_option("-T", '--testpath',
      default="cmake/build_debug/",
      action="store",
      help="The relative directory where the unit tests were built into. (default: 'cmake/build_debug/')",
    )

    parser.add_option("--debug",
      action="store_true", default=False,
      help="Enable debug mode: skip test running, only generate summary from previously created JSON and log files.",
    )

    parser.add_option("--timeout",
      action="store",type="float", default=6,
      help="Automatically cancel test runs if they are taking longer than the given number of minutes. A value of zero turns off the timeout. Note that some tests might have custom (longer) timeout value for historical reasons."
    )

    parser.add_option("--valgrind",
      default=False, action="store_true",
      help="Enable valgrind checking mode, instead of regular unit testing.",
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

    (options, args) = parser.parse_args(args=args[1:])

    options.test_list = args
    # Need old style to run colon-containing tests, due to lack of json output
    if not options.one and any( ':' in t for t in options.test_list ):
        print("A subtest (with a colon) was specified: automatically turning on one-by-one running.")
        options.one = True

    #When specifiying more than one test (or the default of all tests), mute tests unless overridden
    if len(options.test_list) != 1 and options.mute == [] and options.unmute == [] and options.levels == [] and options.level == []:
        print("More than one test to run - enabling mute.")
        options.mute.append("all")

    if options.database == parser.get_default_values().database:
        if path.isdir( path.join( "..", "database") ):
            options.database = path.abspath(path.join( "..", "database"))

        elif os.environ.get('ROSETTA3_DB') is not None and \
                path.isdir(os.environ.get('ROSETTA3_DB')):
            options.database = os.environ.get('ROSETTA3_DB')
        elif path.isdir( path.join( "..", "database") ):
            options.database = path.abspath(path.join( "..", "database"))
        elif path.isdir( path.join( path.expanduser("~"), "database") ):
            options.database = path.join( path.expanduser("~"), "database")
        elif path.isdir( path.join( "..", "rosetta_database") ):
            options.database = path.abspath(path.join( "..", "rosetta_database"))
        else:
            raise ValueError("Can't find database at %s; please set $ROSETTA3_DB or use -d" % options.database)

    else: options.database = path.abspath( options.database )

    for fl in 'bbdep02.May.sortlib.Dunbrack02.lib.bin bbdep02.May.sortlib-correct.12.2010.Dunbrack02.lib.bin ExtendedOpt1-5/Dunbrack10.lib.bin beta_nov2016/Dunbrack10.lib.bin shapovalov/StpDwn_0-0-0/Dunbrack10.lib.bin '.split():
        if os.path.isfile(options.database + '/rotamer/' + fl): break

    else:
        print('WARNING: Clean checkout detected! (no database binary files found) Adding extra time to timeout time...')
        if options.timeout: options.timeout += 16

    if not options.valgrind and ( options.trackorigins or options.leakcheck or options.valgrind_path is not None ):
        print("Valgrind specific options set, but Valgrind mode not enabled. Enabling it for you.")
        options.valgrind = True
    if options.valgrind:
        parse_valgrind_options(options,parser)

    # Fiddle with timeout defaults
    if len(options.test_list) > 0 and options.timeout == parser.get_default_values().timeout:
        options.timeout = 0 # Don't do timeout (by default) if we've gone to the hassle of manually specifying tests

    global Options;  Options = options

    T = Tester()
    T.runUnitTests()
    if not Options.one : T.printSummary()

    print("Done!")


if __name__ == "__main__": main(sys.argv)
