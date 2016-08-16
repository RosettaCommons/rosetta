#!/usr/bin/env python
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

from popen2 import popen2
import random
import math
from sys import argv

####
##
## Usage: ./compate_boot_differences.py <test.boot1> <control.boot1>
##
## Outputs the counts for better, worse, and no change
##
####



input1 = open(argv[1], 'r', -1).readlines()
input2 = open(argv[2], 'r', -1).readlines()


dict1 = {}
for i in input1:
    line = i.split()

    if len(line) == 3:

        rms = float(line[2])
        stdev = float(line[1])

        tag = line[0]

        dict1[tag] = [rms,stdev]

dict2 = {}
for i in input2:
    line = i.split()

    if len(line) == 3:

        rms = float(line[2])
        stdev = float(line[1])

        tag = line[0]

        dict2[tag] = [rms,stdev]

better = 0
worse = 0
no_change = 0
for i in dict1.keys():
    if dict2.has_key(i):
        diff = dict1[i][0]-dict2[i][0]

        error = math.sqrt(pow(dict1[i][1],2)+pow(dict2[i][1],2))

        if abs(diff) > error*2:
            if diff < 0:
                better += 1
            else:
                worse += 1

        else:
            no_change += 1


print "Better: ",better
print "Worse: ",worse
print "No Change: ",no_change
