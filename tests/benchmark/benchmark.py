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


import os, sys, imp, shutil, json

import argparse

# do not change this wording, they have to stay in sync with upstream (up to benchmark-model).
_ResultCodes_ = dict(
    _Finished_     = '_Finished_',
    _Failed_       = '_Failed_',
    _BuildFailed_  = '_BuildFailed_',
    _ScriptFailed_ = '_ScriptFailed_'
    )


class NT:  # named tuple
    def __init__(self, **entries): self.__dict__.update(entries)
    def __repr__(self):
        r = 'NT: |'
        for i in dir(self):
            if not i.startswith('__') and not isinstance(getattr(self, i), types.MethodType): r += '%s --> %s, ' % (i, getattr(self, i))
        return r[:-2]+'|'


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

    print( 'Running test suites: {}'.format(Options.args) )

    for ts in Options.args:
        if not ts.startswith('tests/'): file_name = 'tests/' + ts + '.py'
        else: file_name = ts;  ts = ts.partition('tests/')[2][:-3]  # removing dir prefix

        #print file_name, ts

        test_suite = imp.load_source('test_suite', file_name)


        working_dir = os.path.abspath('./results/' + ts)
        if os.path.isdir(working_dir): shutil.rmtree(working_dir);  #print('Removing old job dir %s...' % working_dir)  # remove old dir if any
        os.makedirs(working_dir)

        res = test_suite.run_test_suite(rosetta_dir=os.path.abspath('../..'), working_dir=working_dir, jobs=Options.jobs, verbose=True)


        if res['suite_result'] not in _ResultCodes_.values(): print 'Warning!!! TestSuite {} failed with unknow result code: {}'.format(t, res['suite_result'])
        else: print 'Test {} finished with output:\n{}'.format(ts, json.dumps(res, sort_keys=True, indent=2))

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
