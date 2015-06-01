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


import sys, signal

import os, commands, re, subprocess, time
from os import path
from optparse import OptionParser


class OI:
    def __init__(self, **entries): self.__dict__.update(entries)


class Runner:
    def __init__(self):
        self.jobs = []    # list of spawned process pid's
        self.output = ''  # output of current job

    def log(self, message):
        self.output += message
        if not Options.quiet: print message

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


    def signal_handler(self, signal_, f):
        print 'Ctrl-C pressed... killing child jobs...'
        for pid in self.jobs:
            os.killpg(os.getpgid(pid), signal.SIGKILL)


    def runCommandLine(self, file_, line_, command_line):
        self.log("Running %s:%s %s" % (file_, line_, command_line) )

        f = subprocess.Popen(command_line, bufsize=0, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).stderr
        for line in f:
            self.log(line)
            sys.stdout.flush()
        f.close()


    def run(self, files):
        signal.signal(signal.SIGINT, self.signal_handler)
        for f in files:
            prefix = Options.prefix
            if len(files) > 1  and  prefix:
                prefix += '_' + f

            for i, l in enumerate(file(f)):
                pid = self.mfork()
                if not pid:  # we are child process
                    self.runCommandLine(f, i, l)
                    if prefix: file(prefix + '_%.4d' % i, 'w').write(self.output)
                    sys.exit(0)

        for p in self.jobs: os.waitpid(p, 0)  # waiting for all child process to termintate...



def main(args):
    ''' Script to run Jobs in parallel.
    '''
    parser = OptionParser(usage="usage: %prog [OPTIONS] file_with_command_lines [file2] [file3] ...")
    parser.set_description(main.__doc__)


    parser.add_option("-j", "--jobs",
      default=1,
      type="int",
      help="Number of processors to use when running testss (default: 1)",
    )

    parser.add_option("-p", '--prefix',
      default='',
      action="store",
      help="Specify the prefix for files name where output is saved. Default is '' - which mean no output is saved.",
    )

    parser.add_option('-q', "--quiet", action="store_true", dest="quiet", default=False,
      help="Suppress (mute) output to std::out."
    )


    (options, args) = parser.parse_args(args=args[1:])

    global Options;  Options = options

    if len(args) > 0:
        tests = args
    else:
        print "Must specify file name with command lines to run!"
        sys.exit(1)

    R = Runner()
    R.run(args)



if __name__ == "__main__": main(sys.argv)
