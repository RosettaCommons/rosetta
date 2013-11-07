#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   tests/__init__.py
## @brief  Common constats and types for all test types
## @author Sergey Lyskov

# ⚔ do not change wording below, it have to stay in sync with upstream (up to benchmark-model).
# Copied from benchmark-model, standard state code's for tests results.

__all__ = ['_S_Values_', '_S_draft_', '_S_queued_', '_S_running_', '_S_finished_', '_S_failed_', '_S_build_failed_', '_S_script_failed_',
           '_StateKey_', '_ResultsKey_', '_LogKey_'
]


_S_draft_                 = 'draft'
_S_queued_                = 'queued'
_S_running_               = 'running'
_S_finished_              = 'finished'
_S_failed_                = 'failed'
_S_build_failed_          = 'build failed'
_S_script_failed_         = 'script failed'
_S_queued_for_comparison_ = 'queued for comparison'

_S_Values_ = [_S_draft_, _S_queued_, _S_running_, _S_finished_, _S_failed_, _S_build_failed_, _S_script_failed_, _S_queued_for_comparison_]

_StateKey_    = 'state'
_ResultsKey_  = 'results'
_LogKey_      = 'log'


# Standard funtions and classes below ---------------------------------------------------------------------------------

class BenchmarkError(Exception): pass


class NT:  # named tuple
    def __init__(self, **entries): self.__dict__.update(entries)
    def __repr__(self):
        r = 'NT: |'
        for i in dir(self):
            if not i.startswith('__') and not isinstance(getattr(self, i), types.MethodType): r += '%s --> %s, ' % (i, getattr(self, i))
        return r[:-2]+'|'


def Tracer(verbose=False):
    def print_(x): print x
    return print_ if verbose else lambda x: None


def execute(message, commandline, return_=False, untilSuccesses=False):
    TR = Tracer()
    TR(message);  TR(commandline)
    while True:
        (res, output) = commands.getstatusoutput(commandline)
        TR(output)

        if res and untilSuccesses: pass  # Thats right - redability COUNT!
        else: break

        print "Error while executing %s: %s\n" % (message, output)
        print "Sleeping 60s... then I will retry..."
        time.sleep(60)

    if return_ == 'tuple': return(res, output)

    if res:
        TR("\nEncounter error while executing: " + commandline )
        if return_==True: return True
        else: raise BenchmarkBuildError()

    if return_ == 'output': return output
    else: return False
