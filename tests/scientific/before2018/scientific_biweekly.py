#!/usr/bin/env python
# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

# For help, run with -h.  Requires Python 2.4+, like the unit tests.

import sys
if not hasattr(sys, "version_info") or sys.version_info < (2,4):
    raise ValueError("Script requires Python 2.4 or higher!")

import os, shutil, threading, subprocess, signal, time, re, random
from os import path
from optparse import OptionParser, IndentedHelpFormatter


def main(argv):
    '''
A simple system for running scientific tests on Rosetta.

Each test has its own subdirectory in tests/, which serves as its name.
Each test has a file named "command", which should have the command to be executed.
Variable substitution is possible using Python printf format; see examples.
To create a new test, just make a new directory with a "command" file and other needed files.
See tests/HOW_TO_MAKE_TESTS for more information.

When the tests run, a statistics directory will be created which copys the tests directory.
When the tests finishes, there is a staResult in each test subdirectory of statistics. The staResult is
a plain text format for now and will be modified to be yaml format in the future.

Intended usage is to run the tests once with a clean checkout, and again after changes.
The exact results of many runs are hardware-dependent, so "expected" results cannot be kept in SVN.

EXAMPLES:

rm -rf statistics/; ./scientific.py    # create reference results using only default settings

./scientific_biweekly.py    # test reference results against new results

./scientific_biweekly.py -j2    # again, but using 2 processors and custom database location

./scientific_biweekly.py  monomer_ddg   # only run the "monomer_ddg" test
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
      help="Path to Rosetta database. (default: $ROSETTA3_DB, ../../database)",
    )

    parser.add_option("-m", "--mini_home",
      #default=path.join( path.expanduser("~"), "mini"),
      default= path.join( path.dirname( path.dirname( path.dirname(path.abspath(sys.argv[0])) ) ), 'source'),
      help="Directory where Rosetta is found (default: ../../source/)",
    )
    parser.add_option("-j", "--num_procs",
      default=1,
      type="int",
      help="number of processors to use on local machine (default: 1)",
    )
    parser.add_option("--host",
      default=[],
      action="append",
      help="ssh to HOST to run some tests. Assumes shared filesystem. May specify multiple times. May specify HOST or USER@HOST.",
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
      help="Maximum runtime for each test, in minutes (default: no limit)",
      metavar="MINUTES",
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
      default="default",
      dest="extras",
      help="in selecting binaries, which options were specified? (default: default)",
    )
    parser.add_option("--score-function",
      default="score12prime",
      dest="scfxn",
      help="Apply a score function parametrization in score_functions/ to the benchmark (default: score12prime)",
    )
    parser.add_option("--daemon", action="store_true", dest="daemon", default=False,
      help="generate daemon friendly output (off by default)"
    )
    parser.add_option("--yaml",
      default=None,
      help="Save results to specified file in YAML format. (default: None)",
    )
    parser.add_option("--output-dir",
      default="biweekly_statistics",
      #dest="outdir",
      help="Where the output results should be written. (default: biweekly_statistics)"
    )
    parser.add_option("--clean-output-dir",
      default=False,
      action="store_true",
      dest="clean_outdir",
      help="If the output directory exists should it be cleaned before running the test? (default: False)"
    )

    (options, args) = parser.parse_args(args=argv)

    # this is for Baker lab shortcur
    if options.digs > 0:
        options.num_procs = 0 # don't use local processors too, b/c local *is* a dig
        # From the "freewhip" script on the whips themselves:
        # digs = ["dig01","dig02","dig03","dig04","dig05","dig06","dig07","dig08", "dig09","dig11","dig12","dig13","dig14","dig15","dig16", "dig17","dig19","dig20","dig21","dig22","dig23", "dig25","dig26"] * 2
        digs = [ "dig"+str(x) for x in range(1,33) ] * 4 # we use up to 4 processors per box
        assert options.digs <= len(digs)
        random.shuffle(digs)
        options.host.extend( digs[:options.digs] )

    if options.database == parser.get_default_values().database:
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

    # Make sure the current directory is the script directory:
    # Using argv[] here causes problems when people try to run the script as "python scientific.py ..."
    #os.chdir( path.dirname(sys.argv[0]) ) # argv[0] is the script name
    if not path.isdir("biweekly"):
        print "You must run this script from main/tests/scientific/"
        return 2

    #If the output directory directory doesn't exist, create one;
    #else remove the old one and create a new one
    outdir = options.output_dir
    if not path.isdir(outdir): os.mkdir(outdir)
    else:
        if options.clean_outdir:
            print "Cleaning output_dir: %s " % outdir
            shutil.rmtree(outdir)
            os.mkdir(outdir)

    # Each test consists of a directory with a "command" file in it.
    if len(args) > 0:
        tests = args
    else:
        tests = [ d for d in os.listdir("biweekly") if not d.startswith(".") and path.isdir(path.join("biweekly", d)) ]

    # start to run the tests
    queue = Queue()

    for test in tests:
        if not path.isdir( path.join("biweekly", test)):
            print "ERROR: The specified scientific benchmark '%s' must correspond to the directory:" % (test)
            print "ERROR:    biweekly/%s" % (test)
            print "ERROR: Available scientific bencharks are:"
            for test_dir in sorted(os.listdir( "biweekly" )):
                print "ERROR:    %s" % test_dir
            return 1
    

    for test in tests:
        queue.put(test)
        if path.isdir(path.join(outdir, test)): shutil.rmtree(path.join(outdir, test))
        copytree( path.join("biweekly", test), path.join(outdir, test),
                  accept=lambda src, dst: path.basename(src) != '.svn' )
    
    if not path.isdir(path.join("score_functions", options.scfxn)) or \
            not path.isfile(path.join("score_functions", options.scfxn, "flags")) or \
            not path.isfile(path.join("score_functions", options.scfxn, "weights.wts")):
        print "ERROR: The specified score function '%s' must correspond to the directory:" % (options.scfxn)
        print "ERROR:   'score_functions/%s'" % (options.scfxn)
        print "ERROR: containing the following files:"
        print "ERROR:   'flags' : a flags file specifying command line parameters"
        print "ERROR:   'weights.wts' : a score function weights file"
        print "ERROR: Available score functions are:"
        for scfxn_dir in sorted(os.listdir("score_functions")):
            print "ERROR:   %s" % scfxn_dir
        return 1
    else:
        if path.isdir( path.join(outdir, test, options.scfxn) ):
            shutil.rmtree( path.join(outdir, test, options.scfxn) )
        copytree( path.join("score_functions", options.scfxn), path.join(outdir, test, options.scfxn),
                  accept=lambda src, dst: path.basename(src) != '.svn' )

    # Start worker thread(s)
    for i in range(options.num_procs):
        worker = Worker(queue, outdir, options, timeout_minutes=options.timeout)
        thread = threading.Thread(target=worker.work)
        #thread.setDaemon(True) # shouldn't be necessary here
        thread.start()
    for host in options.host:
        worker = Worker(queue, outdir, options, host=host, timeout_minutes=options.timeout)
        thread = threading.Thread(target=worker.work)
        #thread.setDaemon(True) # shouldn't be necessary here
        thread.start()
        # Wait for them to finish
    queue.join()

    # Read and Analyze results
    '''
        if outdir == "biweekly_statistics":
            print "Just generated 'biweekly_statistic' results."
            if options.daemon:
                print "SUMMARY: TOTAL:%i PASSED:%i FAILED:%i." % (len(tests), len(tests), 0)

            if options.yaml:
                f = file(options.yaml, 'w');  f.write("{total : %s, failed : 0, details : {}}" % len(tests));  f.close()

            results = {}
            for test in tests:
                dir_statistics = path.join("biweekly_statistics", test)
                ## diff returns 0 on no differences, 1 if there are differences
                sta_file = open(path.join( path.join(outdir, test), "staResult"), "r")
                if options.daemon:
                    print "-------"+test+"-------"
                    if( path.getsize(path.join( path.join(outdir, test), "staResult"))):
                        for sta in sta_file.readlines():
                            print sta
                    else:
                        print "sta file was not generated!"

                else:
                    print "-----"+test+"-----"
                    if(path.getsize(path.join( path.join(outdir, test), "staResult"))>0):
                        for sta in sta_file.readlines():
                            print sta
                    else:
                        print "sta file was not generated!"


            if options.daemon:
                print "SUMMARY: TOTAL:%i" % (len(tests))
            else:
                print "All tests passed."

            if options.yaml:
                f = file(options.yaml, 'w')
                f.write("{total : %s}" % (len(tests)))
                f.close()
        else: print "wrong out dir\n"
    '''
    return 0


class Worker:
    def __init__(self, queue, outdir, opts, host=None, timeout_minutes=0):
        self.queue = queue
        self.outdir = outdir
        self.opts = opts
        self.host = host
        self.timeout = timeout_minutes * 60

    def work(self):
        try:
            while True:
                test = self.queue.get_nowait()
                try: # Actually catch exception and ignore it.  Python 2.4 can't use "except" and "finally" together.
                    try: # Make sure job is marked done even if we throw an exception
                        # Variables that may be referrenced in the cmd string:
                        workdir = path.abspath( path.join(self.outdir, test) )
                        minidir = self.opts.mini_home
                        database = self.opts.database
                        bin = path.join(minidir, "bin")
                        pyapps = path.join(minidir, "src", "python", "apps")
                        if sys.platform.startswith("linux"): platform = "linux" # can be linux1, linux2, etc
                        elif sys.platform == "darwin": platform = "macos"
                        else: platform = "_unknown_"
                        compiler = self.opts.compiler
                        mode = self.opts.mode
                        extras = self.opts.extras
                        binext = extras+"."+platform+compiler+mode
                        scfxn = self.opts.scfxn
                        # Read the command from the file "command"
                        cmd = file(path.join(workdir, "command")).read().strip()
                        cmd = cmd % vars() # variable substitution using Python printf style
                        cmd_line_sh = path.join(workdir, "command.sh")
                        f = file(cmd_line_sh, 'w');  f.write(cmd);  f.close() # writing back so test can be easily re-run by user lately...
                        #if "'" in cmd: raise ValueError("Can't use single quotes in command strings!")
                        #print cmd; print
                        if self.host is None:
                            print "Running %s on localhost ..." % test
                            proc = subprocess.Popen(["bash", "-c", cmd])#, cwd=workdir)
                        # Can't use cwd=workdir b/c it modifies *local* dir, not remote dir.
                        else:
                            print "Running %s on %s ..." % (test, self.host)
                            # A horrible hack b/c SSH doesn't honor login scripts like .bash_profile
                            # when executing specific remote commands.
                            # This causes problems with e.g. the custom Python install on the Whips.
                            #we replace the default remote PATH with the current local one.
                            cmd = 'PATH="%s"\n%s' % (os.environ["PATH"], cmd)
                            proc = subprocess.Popen(["ssh", self.host, cmd])#, cwd=workdir)
                        if self.timeout == 0:
                            retcode = proc.wait() # does this block all threads?
                        else:
                            start = time.time()
                            while time.time() - start <= self.timeout:
                                retcode = proc.poll()
                                if retcode is not None: break
                                time.sleep(1)
                            if retcode is None:
                                print "*** Test %s exceeded the timeout and will be killed!" % test
                                os.kill(proc.pid, signal.SIGTERM)
                    finally: # inner try
                        print "Finished %s" % test
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
    except WindowsError:
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
#################################################################################
#

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
