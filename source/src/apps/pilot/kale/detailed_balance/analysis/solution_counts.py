#!/usr/bin/env python

import os
import argparse
import numpy

from pylab import *

parser = argparse.ArgumentParser()
parser.add_argument('job')
parser.add_argument('--no-plot', dest='plot', action='store_false')
arguments = parser.parse_args()

path = os.path.join(arguments.job, 'solutions.dat')

with open(path) as file:
    lines = file.readlines()
    data = array([int(x) for x in lines])

print data[data > 16]
print max(data)

hist(data, bins=arange(max(data) + 2) - 0.5)
title(arguments.job)
xlim(-0.5, max(data) + 2)
xticks(arange(0, max(data) + 2, 2))

if arguments.plot and not os.fork():
    show()


