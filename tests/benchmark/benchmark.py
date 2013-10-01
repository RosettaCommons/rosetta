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


import os, sys, imp, shutil, json, platform

import argparse

# ⚔ do not change wording below, it have to stay in sync with upstream (up to benchmark-model).
# Copied from benchmark-model, standard state code's for tests results.
_S_Draft_         = 'draft'
_S_Queued_        = 'queued'
_S_Running_       = 'running'
_S_Finished_      = 'finished'
_S_Failed_        = 'failed'
_S_BuildFailed_   = 'build failed'
_S_ScriptFailed_  = 'script failed'

_S_Values_ = [_S_Draft_, _S_Queued_, _S_Running_, _S_Finished_, _S_Failed_, _S_BuildFailed_, _S_ScriptFailed_]

_StateKey_    = 'state'
_ResultsKey_  = 'results'
_LogKey_      = 'log'


class NT:  # named tuple
    def __init__(self, **entries): self.__dict__.update(entries)
    def __repr__(self):
        r = 'NT: |'
        for i in dir(self):
            if not i.startswith('__') and not isinstance(getattr(self, i), types.MethodType): r += '%s --> %s, ' % (i, getattr(self, i))
        return r[:-2]+'|'


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
      help="Number of processors to use on when building. (default: %(default)s)",
    )

    parser.add_argument('args', nargs=argparse.REMAINDER)

    global Options;
    Options = parser.parse_args(args=args[1:])

    print('Platform: {}'.format(Platform))
    print('Running tests: {}'.format(Options.args) )

    for test in Options.args:
        if not test.startswith('tests/')  and  test.count('.'):  # regular Test or a TestSuite?
            suite_name, test_name = test.split('.')
            file_name = 'tests/' + suite_name + '.py'
            test_suite = imp.load_source('test_suite', file_name)

            working_dir = os.path.abspath('./results/' + test)
            if os.path.isdir(working_dir): shutil.rmtree(working_dir);  #print('Removing old job dir %s...' % working_dir)  # remove old dir if any
            os.makedirs(working_dir)

            res = test_suite.run_test(test=test_name, rosetta_dir=os.path.abspath('../..'), working_dir=working_dir, platform=dict(Platform), jobs=Options.jobs, verbose=True)

            if res[_StateKey_] not in _S_Values_: print 'Warning!!! Test {} failed with unknow result code: {}'.format(t, res[_StateKey_])
            else: print 'Test {} finished with output:\n{}'.format(test, json.dumps(res, sort_keys=True, indent=2))


        else:  # TestSuite
            if not test.startswith('tests/'): file_name = 'tests/' + test + '.py'
            else: file_name = test;  test = test.partition('tests/')[2][:-3]  # removing dir prefix

            test_suite = imp.load_source('test_suite', file_name)

            working_dir = os.path.abspath('./results/' + test)
            if os.path.isdir(working_dir): shutil.rmtree(working_dir);  #print('Removing old job dir %s...' % working_dir)  # remove old dir if any
            os.makedirs(working_dir)

            res = test_suite.run_test_suite(rosetta_dir=os.path.abspath('../..'), working_dir=working_dir, platform=dict(Platform), jobs=Options.jobs, verbose=True)

            if res[_StateKey_] not in _ResultCodes_.values(): print 'Warning!!! TestSuite {} failed with unknow result code: {}'.format(t, res[_StateKey_])
            else:
                js = json.dumps(res, sort_keys=True, indent=2)
                file( working_dir + '/.results.json', 'w').write(js)
                print 'Test {} finished with output:\n{}\n[This output is saved in to {}/.results.json]'.format(test, js, working_dir)

        '''
        # Running as individual test... may be later...
        for t in test_suite.get_tests():
            working_dir = os.path.abspath('./results/' + ts.split('/')[-1][:-3]) + '/' + t
            if os.path.isdir(working_dir): shutil.rmtree(working_dir);  #print('Removing old job dir %s...' % working_dir)  # remove old dir if any
            os.makedirs(working_dir)

            res, output = test_suite.run(t, rosetta_dir=os.path.abspath('../..'), working_dir=working_dir, jobs=Options.jobs, verbose=True)

            if res not in _ResultCodes_.values(): print 'Warning!!! Test {} failed with unknow result code: {}'.format(t, res)
            else: print 'Test {} finished with output:\n{}\nAnd result code: {}'.format(t, output, res)
        '''


if __name__ == "__main__": main(sys.argv)
