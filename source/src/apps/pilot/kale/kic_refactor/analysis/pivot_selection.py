#!/usr/bin/env python

# I suspect that there is bias in the algorithm that picks the windows within a 
# loop where KIC gets applied.  This script attempts to recapitulate that 
# algorithm, but it doesn't do it quite right.  I had to change how kic_start 
# is calculated in the else-statement to things to work at all, but I think 
# something went wrong.  I've now probably got an off-by-one bug somewhere.
#
# More work needs to be done to determine whether or not a bias exists.  If it 
# does, I also need to devise an algorithm that gets rid of it.  Some thoughts:
#
# There are a finite number of possible start and end points.  I could 
# sampling each one and pick one using a random number.  This would definitely 
# be free of bias, but it could be hard to map random numbers to end points.
# 
# Actually, this is not completely unbiased.  The constraint that adjacent 
# indices cannot be picked seems to add some bias towards the edges.  I think 
# it's just a plus-one thing, though, nothing systematic.

import argparse
import random
import pylab

parser = argparse.ArgumentParser()
parser.add_argument('--iterations', '-i', type=int, default=10000)
parser.add_argument('--loop', '-l', type=int, nargs=2, default=(0, 10))
arguments = parser.parse_args()

loop_start, loop_end = arguments.loop
kic_starts, kic_ends = [], []
max_seglen = 12

for iteration in range(arguments.iterations):

    if iteration % 2:
        kic_start = random.randint(loop_start, loop_end - 2)
        kic_end = random.randint(
                kic_start + 2, min(kic_start + max_seglen - 1, loop_end))
    else:
        kic_end = random.randint(loop_start + 2, loop_end)
        kic_start = random.randint(
                max(kic_end - max_seglen + 1, loop_start), kic_end - 2)

    kic_starts.append(kic_start)
    kic_ends.append(kic_end)

pylab.hist((kic_starts, kic_ends), normed=True)
pylab.savefig('pivot-selection.pdf')
pylab.show()
