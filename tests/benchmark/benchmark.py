#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   benchmark.py
## @brief  Run arbitrary Rosetta testing script
## @author Sergey Lyskov


import os, sys, imp, shutil, json, platform, commands

import argparse

from tests import *  # Tests states and key names

# Calculating value of Platform dict
Platform = {}
if sys.platform.startswith("linux"): Platform['os'] = 'Linux'  # can be linux1, linux2, etc
elif sys.platform == "darwin" :      Platform['os'] = 'Mac'
elif sys.platform == "cygwin" :      Platform['os'] = 'Cygwin'
elif sys.platform == "win32" :       Platform['os'] = 'Windows'
else:                                Platform['os'] = '_unknown_'

#Platform['arch'] = platform.architecture()[0][:2]  # PlatformBits
Platform['compiler'] = 'gcc'



def main(args):
    ''' Script to Run arbitrary Rosetta test
    '''
    parser = argparse.ArgumentParser(usage="Generate pyrosetta distribution.")

    parser.add_argument('-j', '--jobs',
      default=1, type=int,
      help="Number of processors to use on when building. (default: use all avaliable memory)",
    )

    parser.add_argument('-m', '--memory',
      default=0, type=float,
      help="Amount of memory to use (default: use 2Gb per job",
    )


    parser.add_argument("--extras", default='', help="Specify scons extras separated by ',': like --extras=mpi,static" )

    parser.add_argument("--debug", action="store_true", dest="debug", default=False, help="Run specified test in debug mode (not with debug build!) this mean different things and depend on the test. Could be: skip the build phase, skip some of the test phases and so on. [off by default]" )

    parser.add_argument("--suffix", default='', help="Specify ending suffix for test output dir. This is useful when you want to save test results in different dir for later comparison." )

    parser.add_argument("--compare", nargs=2, help="Do not run the tests but instead compare precious results. Use --compare suffix1 suffix2" )

    parser.add_argument('args', nargs=argparse.REMAINDER)

    global Options;
    Options = parser.parse_args(args=args[1:])

    if Options.suffix: Options.suffix = '.' + Options.suffix

    Platform['extras'] = Options.extras.split(',') if Options.extras else []

    if Options.memory: memory = Options.memory
    elif Platform['os'] == 'Linux': memory = float(commands.getoutput('free -m').split('\n')[1].split()[1]) / 1024
    elif Platform['os'] == 'Mac':   memory = float(commands.getoutput('sysctl -a | grep hw.memsize').split()[2]) / 1024 / 1024 / 1024

    config = dict(cpu_count=Options.jobs, memory=memory)
    print('Config:{}, Platform:{}'.format(json.dumps(config, sort_keys=True), Platform))

    if Options.compare: print('Comparing tests {} with suffixes: {}'.format(Options.args, Options.compare) )
    else: print('Running tests: {}'.format(Options.args) )


    # def compare(fun):
    #     working_dir_1 = os.path.abspath('./results/' + test + '.' + Options.compare[0])
    #     working_dir_2 = os.path.abspath('./results/' + test + '.' + Options.compare[1])
    #     res_1 = json.load( file(working_dir_1 + '/.results.json') )
    #     res_2 = json.load( file(working_dir_1 + '/.results.json') )
    #     res = fun(res_1, working_dir_1, res_2, working_dir_2)
    #     print json.dumps(res, sort_keys=True, indent=2)
    #     sys.exit(0)


    for test in Options.args:
        if test.startswith('tests/'): test = test.partition('tests/')[2][:-3]  # removing dir prefix and .py suffix

        if test.count('.'): suite_name, test_name = test.split('.')
        else: suite_name, test_name = test, ''

        file_name = 'tests/' + suite_name + '.py'
        test_suite = imp.load_source('test_suite', file_name)

        if Options.compare:
            working_dir_1 = os.path.abspath('./results/' + test + '.' + Options.compare[0])
            working_dir_2 = os.path.abspath('./results/' + test + '.' + Options.compare[1])
            res_1 = json.load( file(working_dir_1 + '/.results.json') )
            res_2 = json.load( file(working_dir_1 + '/.results.json') )
            res = test_suite.compare(test, res_1, working_dir_1, res_2, working_dir_2)

            with file(working_dir_1+'/.compare.json', 'w') as f: json.dump(res, f, sort_keys=True, indent=2)

            print('Comparison finished with results:\n{}'.format( json.dumps(res, sort_keys=True, indent=2) ) )

            if 'summary' in res: print('Summary section:\n{}'.format( json.dumps(res['summary'], sort_keys=True, indent=2) ) )
            print 'Output this comparison saved to {0}/.compare.json'.format(working_dir_1)


        else:
            working_dir = os.path.abspath('./results/' + test + Options.suffix)
            if os.path.isdir(working_dir): shutil.rmtree(working_dir);  #print('Removing old job dir %s...' % working_dir)  # remove old dir if any
            os.makedirs(working_dir)

            try: api_version = test_suite._api_version_
            except AttributeError as e: api_version = ''

            if api_version < '1.0':
                res = test_suite.run(test=test_name, rosetta_dir=os.path.abspath('../..'), working_dir=working_dir, platform=dict(Platform), jobs=Options.jobs, verbose=True, debug=Options.debug)
            else:
                res = test_suite.run(test=test_name, rosetta_dir=os.path.abspath('../..'), working_dir=working_dir, platform=dict(Platform), config=config, verbose=True, debug=Options.debug)

            if res[_StateKey_] not in _S_Values_: print 'Warning!!! Test {} failed with unknow result code: {}'.format(t, res[_StateKey_])
            else: print 'Test {} finished with output:\n{}'.format(test, json.dumps(res, sort_keys=True, indent=2))

            print 'Output and log of this test saved to {0}/.results.json and {0}/.output.log'.format(working_dir)

            with file(working_dir+'/.output.log', 'w') as f: f.write(res[_LogKey_])
            with file(working_dir+'/.results.json', 'w') as f: json.dump(res, f, sort_keys=True, indent=2)


        # if not test.startswith('tests/')  and  test.count('.'):  # regular Test or a TestSuite?
        #     suite_name, test_name = test.split('.')
        #     file_name = 'tests/' + suite_name + '.py'
        #     test_suite = imp.load_source('test_suite', file_name)

        #     if Options.compare: compare(test_suite.compare_tests)
        #     else:
        #         working_dir = os.path.abspath('./results/' + test + Options.suffix)
        #         if os.path.isdir(working_dir): shutil.rmtree(working_dir);  #print('Removing old job dir %s...' % working_dir)  # remove old dir if any
        #         os.makedirs(working_dir)

        #         res = test_suite.run_test(test=test_name, rosetta_dir=os.path.abspath('../..'), working_dir=working_dir, platform=dict(Platform), jobs=Options.jobs, verbose=True, debug=Options.debug)

        #         if res[_StateKey_] not in _S_Values_: print 'Warning!!! Test {} failed with unknow result code: {}'.format(t, res[_StateKey_])
        #         else: print 'Test {} finished with output:\n{}'.format(test, json.dumps(res, sort_keys=True, indent=2))


        # else:  # TestSuite
        #     if not test.startswith('tests/'): file_name = 'tests/' + test + '.py'
        #     else: file_name = test;  test = test.partition('tests/')[2][:-3]  # removing dir prefix

        #     test_suite = imp.load_source('test_suite', file_name)

        #     if Options.compare: compare(test_suite.compare_test_suites)
        #     else:
        #         working_dir = os.path.abspath('./results/' + test + Options.suffix)
        #         if os.path.isdir(working_dir): shutil.rmtree(working_dir);  #print('Removing old job dir %s...' % working_dir)  # remove old dir if any
        #         os.makedirs(working_dir)

        #         res = test_suite.run_test_suite(rosetta_dir=os.path.abspath('../..'), working_dir=working_dir, platform=dict(Platform), jobs=Options.jobs, verbose=True, debug=Options.debug)

        #         if res[_StateKey_] not in _S_Values_: print 'Warning!!! TestSuite {} failed with unknow result code: {}'.format(t, res[_StateKey_])
        #         else:
        #             print 'Test {} finished with output:\n{}\n'.format(test, json.dumps(res, sort_keys=True, indent=2))


if __name__ == "__main__": main(sys.argv)
