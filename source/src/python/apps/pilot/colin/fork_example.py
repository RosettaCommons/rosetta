#!/usr/bin/env python

import rosetta_py.utility.ForkManager as ForkManager
import time

# This is an example of parallelizing the following for loop with ForkManager:
#
#for i in xrange(100):
#
#	jobname = "job" + str(i)
#	print jobname + " started"
#
#	for t in xrange(1, 6):
#		time.sleep(1)
#		print t, "seconds elapsed..."
#
# Below, 4 lines are added to enable parallelization in 4 slave processes:

# Setting max_slaves to 0 disables parallelization altogether
fm = ForkManager.ForkManager(max_slaves=4, incremental_output = True,
                             only_output_bad = False, time_limit = 600)

for i in xrange(100):

	jobname = "job" + str(i)
	print jobname + " started"

	if fm.fork(jobname):

		# We're in the child process now, do the processing
		for t in xrange(1, 6):
			time.sleep(1)
			print str(t) + " seconds elapsed..."

		fm.exit_slave()

print "Waiting for all slaves to finish"
fm.wait()
