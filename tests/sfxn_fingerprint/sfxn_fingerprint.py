#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

# For help, run with -h.  Requires Python 2.4+, like the unit tests.

import sys
if not hasattr(sys, "version_info") or sys.version_info < (2,4):
    raise ValueError("Script requires Python 2.4 or higher!")

import os, shutil, threading, subprocess, signal, time, re, random, datetime
from os import path
from optparse import OptionParser, IndentedHelpFormatter


def main(argv):
    """
System for running ScoreVestion's tests on Mini.

Each test has its own subdirectory in scores/, which serves as its name.
Each test has a file named "command", which should have the command to be executed.
Variable substitution is possible using Python printf format; see examples.
To create a new test, just make a new directory with a "command" file and other needed files.
See score/EXAMPLE for more information.

When the tests are run, one of two things happens.
(1) If the "ref" directory does not exist, it is created.
Each subdirectory of scores/ is copied to ref/, and the test is executed there.
(2) Otherwise, the "new" directory is wiped out and re-created, and tests are run in new/ instead.
Afterwards, "diff" is used to compare each subdirectory of ref/ with the corresponding subdirectory in new/.
A test is considered passed if there are no differences, and failed otherwise.

Intended usage is to run the tests once with a clean checkout, and again after changes.
The exact results of many runs are hardware-dependent, so "expected" results cannot be kept in SVN.
If tests fail in expected ways, the appropriate subdirectories can be copied from new/ to ref/ to update the expected results.

EXAMPLES:

rm -r ref/; ./ScoreVersion.py    # create reference results using only default settings

./ScoreVersion.py    # test reference results against new results

./ScoreVersion.py -d ~/minidb -j2    # again, but using 2 processors and custom database location

./ScoreVersion.py score12    # only run the "score12" test
    """
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
      help="Path to Rosetta database. (default: $ROSETTA_DB, ~/rosetta_database)",
    )

    parser.add_option("-m", "--mini_home",
      #default=path.join( path.expanduser("~"), "mini"),
      #default=path.dirname( path.dirname( path.dirname(path.abspath(sys.argv[0])) ) ),
      default= path.join( path.dirname( path.dirname( path.dirname(path.abspath(sys.argv[0])) ) ), 'source'),
      help="Directory where Mini is found (default: ../../source/)",
    )
    parser.add_option("--data",
      default = path.dirname( path.dirname(path.abspath(sys.argv[0]) ) + '/data/'),
      help="Directory where data folder for ScoreVersion is found (default: ./data/)",
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
      default="",
      dest="extras",
      help="in selecting binaries, which options were specified? (default: default)",
    )
    #parser.add_option("--daemon", action="store_true", dest="daemon", default=False,
    #  help="generate daemon friendly output (off by default)"
    #)
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

    (options, args) = parser.parse_args(args=argv)

    options.mini_home = path.abspath( options.mini_home )
    print 'Using Rosetta source dir at:', options.mini_home


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
        elif path.isdir( path.join( path.expanduser("~"), "rosetta_database") ):
            options.database = path.join( path.expanduser("~"), "rosetta_database")
        else:
            print "Can't find database at %s; please set $ROSETTA3_DB or use -d" % options.database
            return 1
    # Normalize path before we change directories!
    options.database = path.abspath(options.database)

    # Make sure the current directory is the script directory:
    # Using argv[] here causes problems when people try to run the script as "python integration.py ..."
    #os.chdir( path.dirname(sys.argv[0]) ) # argv[0] is the script name
    if not path.isdir("scores"):
        print "You must run this script from mini/tests/ScoreVersion/"
        return 2

    # If the "ref" directory doesn't exist, compute that;
    # else compute the "new" directory.

    #if not path.isdir("ref"): outdir = "ref"
    #else: outdir = "new"
    # scratch that ^^^^^, output will always go to new dir for now and rename it later.
    if not path.isdir("ref"): rename_to_ref = True
    else: rename_to_ref = False
    outdir = "new"

    if not options.compareonly:
        # Remove everything in the current outdir, then re-create it empty
        if path.isdir(outdir): shutil.rmtree(outdir)
        os.mkdir(outdir)

    # Each test consists of a directory with a "command" file in it.
    if len(args) > 0:
        tests = args
    else:
        tests = [ d for d in os.listdir("scores") if not d.startswith(".") and path.isdir(path.join("scores", d)) ]

    if not options.compareonly:
        queue = Queue()
        queue.TotalNumberOfTasks = len(tests)
        for test in tests:
            queue.put(test)
            #shutil.copytree( path.join("tests", test), path.join(outdir, test) )
            copytree( path.join("scores", test), path.join(outdir, test),
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


    # removing absolute paths to root Rosetta checkout from tests results and replacing it with 'ROSETTA'
    rosetta_dir = os.path.abspath('../..')
    for test in tests:
        for dir_, _, files in os.walk( path.join("new", test) ):
            for f in files:
                if f == 'command.sh': continue
                fname = dir_ + '/' + f
                data = file(fname).read()
                if rosetta_dir in data:
                    with file(fname, 'w') as f: f.write( data.replace(rosetta_dir, 'ROSETTA_MAIN') )

    # Analyze results
    print

    #if outdir == "ref":
    if rename_to_ref:
        os.renames('new', 'ref')

        print "Just generated 'ref' results;  run again after making changes."
        #if options.daemon:
        #    print "SUMMARY: TOTAL:%i PASSED:%i FAILED:%i." % (len(tests), len(tests), 0)

        if options.yaml:
            f = file(options.yaml, 'w');  f.write("{total : %s, failed : 0, details : {}}" % len(tests));  f.close()

    else:
        if options.skip_comparison:
            print 'Skipping comparison phase because command line option "--skip-comparison" was specified...'

        else:
            diffs = 0
            results = {}
            full_log = ''
            for test in tests:
                dir_before = path.join("ref", test)
                dir_after = path.join(outdir, test)
                # diff returns 0 on no differences, 1 if there are differences

                flags = ["-rq"]
                if options.fulldiff: flags = ["-r"]
                flags += ['--exclude=command.sh', '--exclude=*.ignore']

                proc = subprocess.Popen(["diff"] + flags + [dir_before, dir_after], stdout=subprocess.PIPE)

                full_log_msg = "FAIL %s\n" % test


                msg = "FAIL %s\n" % test
                for diff in proc.stdout.readlines():
                    msg += "     %s\n" % diff.strip()
                    full_log_msg += "     %s\n" % diff.strip()

                full_log_msg += '\n'

                result = proc.wait()
                results[test] = result

                if result == 0:
                    print "ok   %s" % test
                    full_log += "ok   %s\n" % test
                else:
                    print msg
                    full_log += full_log_msg

                    diffs += 1

            if diffs:
                print "%i test(s) failed.  Use 'diff' to compare results." % diffs
            else:
                print "All tests passed."

            if options.yaml:
                f = file(options.yaml, 'w')
                brief = makeBriefResults(full_log)
                brief = brief.replace('"', '\\"')
                brief = '"' + brief.replace('\n', '\\n') + '"'
                f.write("{total : %s, failed : %s, details : %s, brief : %s}" % (len(tests), diffs, results, brief) )
                f.close()

    return 0


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
    return res


# Wrap new line simbols inside given strings by using '\' character
def wrapNewLine(s):
    r = ''



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
                    start = time.time() # initial guess at start time, in case of exception
                    try: # Make sure job is marked done even if we throw an exception
                        # Variables that may be referrenced in the cmd string:
                        python = sys.executable
                        workdir = path.abspath( path.join(self.outdir, test) )
                        minidir = self.opts.mini_home
                        database = self.opts.database
                        data = self.opts.data
                        bin = path.join(minidir, "bin")
                        pyapps = path.join(minidir, "src", "python", "apps")
                        if sys.platform.startswith("linux"): platform = "linux" # can be linux1, linux2, etc
                        elif sys.platform == "darwin": platform = "macos"
                        elif sys.platform == "cygwin": platform = "cygwin"
                        else: platform = "_unknown_"
                        compiler = self.opts.compiler
                        mode = self.opts.mode
                        extras = self.opts.extras if self.opts.extras else 'default'
                        binext = extras+"."+platform+compiler+mode
                        # Read the command from the file "command"
                        cmd = file(path.join(workdir, "command")).read().strip()
                        cmd = cmd % vars() # variable substitution using Python printf style
                        cmd_line_sh = path.join(workdir, "command.sh")
                        f = file(cmd_line_sh, 'w');  f.write(cmd);  f.close() # writing back so test can be easily re-run by user lately...
                        #if "'" in cmd: raise ValueError("Can't use single quotes in command strings!")
                        #print cmd; print
                        if self.host is None:
                            print "Running %s on localhost ..." % test
                            proc = subprocess.Popen(["bash",  cmd_line_sh], preexec_fn=os.setpgrp)
                        # Can't use cwd=workdir b/c it modifies *local* dir, not remote dir.
                        else:
                            print "Running %s on %s ..." % (test, self.host)
                            # A horrible hack b/c SSH doesn't honor login scripts like .bash_profile
                            # when executing specific remote commands.
                            # This causes problems with e.g. the custom Python install on the Whips.
                            # So we replace the default remote PATH with the current local one.
                            cmd = 'PATH="%s"\n%s' % (os.environ["PATH"], cmd)
                            proc = subprocess.Popen(["ssh", self.host, cmd], preexec_fn=os.setpgrp)#, cwd=workdir)
                        start = time.time() # refined start time
                        if self.timeout == 0:
                            retcode = proc.wait() # does this block all threads?
                        else:
                            while time.time() - start <= self.timeout:
                                retcode = proc.poll()
                                if retcode is not None: break
                                time.sleep(1)
                            if retcode is None:
                                print "*** Test %s exceeded the timeout and will be killed! [%s]\n" % (test, datetime.datetime.now())
                                #os.kill(proc.pid, signal.SIGTERM)
                                os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
                        if retcode != 0 and retcode is not None:
                            error_string = "*** Test %s did not run!  Check your --mode flag and paths. [%s]\n" % (test, datetime.datetime.now())
                            print error_string,

                            # Writing error_string to a file, so integration test should fail for sure
                            file(path.join(workdir, ".test_did_not_run.log"), 'w').write(error_string)


                    finally: # inner try
                        percent = (100* (self.queue.TotalNumberOfTasks-self.queue.qsize())) / self.queue.TotalNumberOfTasks
                        print "Finished %s in %i seconds\t [~%s test (%s%%) started, %s in queue]" % (test, time.time() - start, self.queue.TotalNumberOfTasks-self.queue.qsize(), percent, self.queue.qsize())
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
    #except WindowsError:
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
################################################################################


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
