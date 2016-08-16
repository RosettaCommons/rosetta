#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

import sys
if not hasattr(sys, "version_info") or sys.version_info < (2,4):
    raise ValueError("Script requires Python 2.4 or higher!")

import os, threading, subprocess, signal, time, re, random
from os import path
from optparse import OptionParser, IndentedHelpFormatter


def main(argv):
    '''
A very simple system for executing multiple commands in parallel on one machine (or the digs).
Just supply a bunch of commands on standard input, and the number of cores to use.
For a (stupid) example:

    (for((i=0;i<2;i++)); do echo sleep 5; done) | time parallel.py

The script *can* be interrupted with Ctrl-C, but it's delicate --
so only hit Ctrl-C *once*, or you may leave some orphan jobs running.
    '''
    parser = OptionParser(usage="usage: %prog [OPTIONS] < cmds.txt", formatter=ParagraphHelpFormatter())
    parser.set_description(main.__doc__)
    # parser.add_option("-short", ["--long"],
    #   action="store|store_true|store_false|append",
    #   default=True|False|...
    #   type="string|int|float",
    #   dest="opt_name",
    #   help="store value in PLACE",
    #   metavar="PLACE",
    # )
    parser.add_option("-j", "--num-procs",
      default=2,
      type="int",
      help="number of processors to use on local machine (default: 2)",
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

    (options, args) = parser.parse_args(args=argv)

    if options.digs > 0:
        options.num_procs = 0 # don't use local processors too, b/c local *is* a dig
        digs = [ "dig"+str(x) for x in range(1,33) ] * 4 # we use up to 4 processors per box
        assert options.digs <= len(digs)
        random.shuffle(digs)
        options.host.extend( digs[:options.digs] )

    if len(args) > 0:
        parser.print_help()
        print "No command line args, please;  send cmds to stdin."
        return 1

    queue = Queue()
    for cmd in sys.stdin:
        queue.put( cmd.strip() )

    # Start worker thread(s)
    workers = []
    for i in range(options.num_procs):
        worker = Worker(queue, options, timeout_minutes=options.timeout)
        workers.append(worker)
        thread = threading.Thread(target=worker.work)
        #thread.setDaemon(True) # shouldn't be necessary here
        thread.start()
    for host in options.host:
        worker = Worker(queue, options, host=host, timeout_minutes=options.timeout)
        workers.append(worker)
        thread = threading.Thread(target=worker.work)
        #thread.setDaemon(True) # shouldn't be necessary here
        thread.start()

    # Wait for them to finish
    try:
        queue.join()
    except KeyboardInterrupt, e:
        # If interrupted by Ctrl-C, ask workers to kill their running jobs
        for worker in workers:
            worker.abort()
        raise e

    return 0


class Worker:
    def __init__(self, queue, opts, host=None, timeout_minutes=0):
        self.queue = queue
        self.opts = opts
        self.host = host
        self.timeout = timeout_minutes * 60
        self._abort = threading.Event()

    def abort(self):
        self._abort.set()

    def work(self):
        try:
            while True:
                cmd = self.queue.get_nowait()
                try: # Actually catch exception and ignore it.  Python 2.4 can't use "except" and "finally" together.
                    try: # Make sure job is marked done even if we throw an exception
                        workdir = path.abspath(".")
                        if self.host is None:
                            print "Running %s on localhost ..." % cmd
                            proc = subprocess.Popen(["bash", "-c", cmd], preexec_fn=os.setpgrp)
                        else:
                            print "Running %s on %s ..." % (cmd, self.host)
                            # A horrible hack b/c SSH doesn't honor login scripts like .bash_profile
                            # when executing specific remote commands.
                            # This causes problems with e.g. the custom Python install on the Whips.
                            # So we replace the default remote PATH with the current local one.
                            cmd = 'cd %s\nPATH="%s"\n%s' % (workdir, os.environ["PATH"], cmd)
                            proc = subprocess.Popen(["ssh", self.host, cmd], preexec_fn=os.setpgrp)
                        start = time.time()
                        while self.timeout == 0 or time.time() - start <= self.timeout:
                            retcode = proc.poll()
                            if retcode is not None: break
                            if self._abort.isSet(): break
                            time.sleep(1)
                        if retcode is None:
                            print "*** Cmd %s exceeded the timeout and will be killed!" % cmd
                            #os.kill(proc.pid, signal.SIGTERM)
                            os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
                    finally: # inner try
                        print "Finished %s" % cmd
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

################################################################################
# Python 2.4 lacks support for join() / task_done() in the Queue class,
# so I pasted the 2.5 implementation here.
# Also, I fixed a "bug" in join() that prevents Ctrl-C interrupts.

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
                # Using a timeout value here (any value) allows keyboard interruption with Ctrl-C
                #self.all_tasks_done.wait()
                self.all_tasks_done.wait(1000)
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
