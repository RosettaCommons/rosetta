#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   benchmark.py
## @brief  Run arbitrary Rosetta testing script
## @author Sergey Lyskov


import os, os.path, sys, imp, shutil, json, platform, commands
import codecs

from ConfigParser import ConfigParser
import argparse

from tests import *  # Tests states and key names
from hpc_drivers import *


# Calculating value of Platform dict
Platform = {}
if sys.platform.startswith("linux"):
    Platform['os'] = 'ubuntu' if os.path.isfile('/etc/lsb-release') and 'Ubuntu' in open('/etc/lsb-release').read() else 'linux'  # can be linux1, linux2, etc
elif sys.platform == "darwin" :      Platform['os'] = 'mac'
elif sys.platform == "cygwin" :      Platform['os'] = 'cygwin'
elif sys.platform == "win32" :       Platform['os'] = 'windows'
else:                                Platform['os'] = 'unknown'

#Platform['arch'] = platform.architecture()[0][:2]  # PlatformBits
Platform['compiler'] = 'gcc' if Platform['os'] == 'linux' else 'clang'

Platform['python'] = sys.executable

def print_(x): print x


def main(args):
    ''' Script to Run arbitrary Rosetta test
    '''
    parser = argparse.ArgumentParser(usage="Main testing script to run tests in the tests directory. "
                                           "Use the --skip-compile to skip the build phase when testing locally. "
                                           "Example Command: /benchmark.py -j2 integration.valgrind")

    parser.add_argument('-j', '--jobs',
      default=1, type=int,
      help="Number of processors to use on when building. (default: use all avaliable memory)",
    )

    parser.add_argument('-m', '--memory',
      default=0, type=int,
      help="Amount of memory to use (default: use 2Gb per job",
    )

    parser.add_argument('--compiler', default=Platform['compiler'], help="Compiler to use")

    parser.add_argument('--python', default=('{}.{}'.format(*sys.version_info) if Platform['os'] == 'mac' else '3.6'), help="Python interpreter to use")

    parser.add_argument("--extras", default='', help="Specify scons extras separated by ',': like --extras=mpi,static" )

    #parser.add_argument("--options", default='', help="""Specify JSON string of for platform options: like --options='{"py":"monolith"}'""" )

    parser.add_argument("--debug", action="store_true", dest="debug", default=False, help="Run specified test in debug mode (not with debug build!) this mean different things and depend on the test. Could be: skip the build phase, skip some of the test phases and so on. [off by default]" )

    parser.add_argument("--suffix", default='', help="Specify ending suffix for test output dir. This is useful when you want to save test results in different dir for later comparison." )

    parser.add_argument("--compare", nargs=2, help="Do not run the tests but instead compare previous results. Use --compare suffix1 suffix2" )

    parser.add_argument("--config", default='benchmark.{os}.ini'.format(os=Platform['os']), action="store", help="Location of .ini file with additional options configuration. Optional.")

    parser.add_argument("--skip-compile", dest='skip_compile', default=None, action="store_true", help="Skip the compilation phase. Assumes the binaries are already compiled locally.")

    parser.add_argument('args', nargs=argparse.REMAINDER)

    global Options;
    Options = parser.parse_args(args=args[1:])

    if any( [a.startswith('-') for a in Options.args] ) :
        print '\nWARNING WARNING WARNING WARNING\n'
        print '\tInterpreting', ' '.join(["'"+a+"'" for a in Options.args if a.startswith('-')]), 'as test name(s), rather than as option(s).'
        print "\tTry moving it before any test name, if that's not what you want."
        print '\nWARNING WARNING WARNING WARNING\n'

    if Options.suffix: Options.suffix = '.' + Options.suffix

    Platform['extras'] = Options.extras.split(',') if Options.extras else []
    Platform['python'] = Options.python
    #Platform['options'] = json.loads( Options.options ) if Options.options else {}

    if Options.memory: memory = Options.memory
    elif Platform['os'] in ['linux', 'ubuntu']: memory = int(commands.getoutput('free -m').split('\n')[1].split()[1]) / 1024
    elif Platform['os'] == 'mac':   memory = int(commands.getoutput('sysctl -a | grep hw.memsize').split()[1]) / 1024 / 1024 / 1024

    Platform['compiler'] = Options.compiler

    if os.path.isfile(Options.config):
        Config = ConfigParser( dict(here=os.path.abspath('./') ) )
        Config.readfp( file(Options.config) )

    else:
        Config = ConfigParser()
        Config.set('DEFAULT', 'hpc_driver', 'MultiCore')
        Config.set('DEFAULT', 'branch',     'unknown')
        Config.set('DEFAULT', 'revision',   '42')
        Config.set('DEFAULT', 'user_name',  'Jane Roe')
        Config.set('DEFAULT', 'user_email', 'jane.roe@university.edu')
        Config.add_section('config')

    Config.set('DEFAULT', 'cpu_count', str(Options.jobs) )
    Config.set('DEFAULT', 'memory',    str(memory) )

    config = Config.items('config')
    config = dict(config, cpu_count=Options.jobs, memory=memory)
    if 'prefix' not in config: config['prefix'] = os.path.abspath('./results/prefix')
    if Options.skip_compile is not None:
        config['skip_compile'] = Options.skip_compile

    print('Config:{}, Platform:{}'.format(json.dumps(config, sort_keys=True), Platform))

    if Options.compare: print('Comparing tests {} with suffixes: {}'.format(Options.args, Options.compare) )
    else: print('Running tests: {}'.format(Options.args) )


    # def compare(fun):
    #     working_dir_1 = os.path.abspath('./results/' + test + '.' + Options.compare[0])
    #     working_dir_2 = os.path.abspath('./results/' + test + '.' + Options.compare[1])
    #     res_1 = json.load( file(working_dir_1 + '/output.results.json') )
    #     res_2 = json.load( file(working_dir_1 + '/output.results.json') )
    #     res = fun(res_1, working_dir_1, res_2, working_dir_2)
    #     print json.dumps(res, sort_keys=True, indent=2)
    #     sys.exit(0)


    for test in Options.args:
        if test.startswith('tests/'): test = test.partition('tests/')[2][:-3]  # removing dir prefix and .py suffix

        suite_name = ''
        while test:
            suite_name_add_on, _, test_name = test.partition('.')
            #print 'suite_name: {suite_name}, suite_name_add_on: {suite_name_add_on}'.format(**vars())
            suite_name += '/' + suite_name_add_on

            file_name = 'tests' + suite_name + '.py'
            if os.path.isfile(file_name): break

            file_name = 'tests' + suite_name + '/command.py'
            if os.path.isfile(file_name): break

            test = test_name

        print 'Loading test from:', file_name
        test_suite = imp.load_source('test_suite', file_name)

        if Options.compare:
            working_dir_1 = os.path.abspath('./results/' + test + '.' + Options.compare[0])
            res_1 = json.load( file(working_dir_1 + '/output.results.json') )["results"]

            if Options.compare[1] == 'None': res_2, working_dir_2 = None, None
            else:
                working_dir_2 = os.path.abspath('./results/' + test + '.' + Options.compare[1])
                res_2 = json.load( file(working_dir_2 + '/output.results.json') )["results"]

            res = test_suite.compare(test, res_1, working_dir_1, res_2, working_dir_2)

            with file(working_dir_1+'/output.compare.json', 'w') as f: json.dump(res, f, sort_keys=True, indent=2)

            print('Comparison finished with results:\n{}'.format( json.dumps(res, sort_keys=True, indent=2) ) )

            if 'summary' in res: print('Summary section:\n{}'.format( json.dumps(res['summary'], sort_keys=True, indent=2) ) )
            print 'Output this comparison saved to {0}/output.compare.json'.format(working_dir_1)


        else:
            working_dir = os.path.abspath('./results/' + test + Options.suffix)
            if os.path.isdir(working_dir): shutil.rmtree(working_dir);  #print('Removing old job dir %s...' % working_dir)  # remove old dir if any
            os.makedirs(working_dir)

            hpc_driver = eval(config['hpc_driver'] + '_HPC_Driver')(working_dir, Config, tracer=lambda x: print_(x), set_daemon_message=lambda x:None)

            api_version = test_suite._api_version_ if hasattr(test_suite, '_api_version_') else ''

            if api_version < '1.0':
                res = test_suite.run(test=test_name, rosetta_dir=os.path.abspath('../..'), working_dir=working_dir, platform=dict(Platform), jobs=Options.jobs, verbose=True, debug=Options.debug)
            else:
                res = test_suite.run(test=test_name, rosetta_dir=os.path.abspath('../..'), working_dir=working_dir, platform=dict(Platform), config=config, hpc_driver=hpc_driver, verbose=True, debug=Options.debug)

            if res[_StateKey_] not in _S_Values_: print 'Warning!!! Test {} failed with unknow result code: {}'.format(t, res[_StateKey_])
            else: print 'Test {} finished with output:\n{}'.format(test, json.dumps(res, sort_keys=True, indent=2))

            print 'Output and log of this test saved to {0}/output.results.json and {0}/output.log'.format(working_dir)

            # Caution! Some of the strings in the result object may be unicode.
            with codecs.open(working_dir+'/output.log'.format(test), 'w', encoding='utf-8', errors='replace') as f: # Be robust to unicode in the log messages
                    f.write(res[_LogKey_])
            # Json by default serializes to an ascii-encoded format
            with file(working_dir+'/output.results.json', 'w') as f: json.dump(res, f, sort_keys=True, encoding="utf-8", indent=2)


if __name__ == "__main__": main(sys.argv)
