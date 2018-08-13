#!/usr/bin/env python
# :noTabs=true:
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
## Usage: ./analysis_script.py <file containing list of score files>
##
## Outputs a file containing best 1% rms and standard deviation values for each file in the input list
##
###

input = open(argv[1], 'r', -1)
list = input.readlines()
input.close()

proteins = {}

names = []

for i in list:

    na = i[:-1]

    bench = open(na,'r',-1).readlines()

    rms = 0
    if bench[0].split()[0] == "SEQUENCE:":
        for i in range(len(bench[1].split())):
            if (bench[1].split()[i] == "rms"):
                rms = i
                l = len(bench[0].split()) + 1
    else:
        for i in range(len(bench[0].split())):
            if (bench[0].split()[i] == "rms"):
                rms = i
                l = len(bench[0].split()) + 1

    for i in bench:
        line = i.split()

        if (line[0] == 'SCORE:') and (line[rms] != 'NULL') and (line[1] != 'SCORE') and (line[1] != 'score'):

            if proteins.has_key(na):
                proteins[na].append(float(line[rms]))
            else:
                names.append(na)
                proteins[na] = [float(line[rms])]

names.sort()

#############BOOTSTRAPPING

random.seed()

bn = 100 # number of boots

percentile = float(10)/float(1000) # 1% leading edge

boot = {}
bresults = {}
results = {}
n = 0

for i in names:
    boot[i] = []
    bresults[i] = []
    results[i] = -1
    n += 1

    for j in range(bn):
        boot[i].append([])

    for k in proteins[i]:
        for j in range(bn):
            boot[i][j].append(random.choice(proteins[i]))

for i in names:
    for b in boot[i]:
        b.sort()
        total = 0
        num = int(float(len(b))*float(percentile))
        if percentile == 0:
            num = 1
        for k in range(num):
            total += float(b[k])
        bresults[i].append(float(total)/(float(num) + 1e-200 ) )

for i in names:
    b = proteins[i]
    b.sort()
    total = 0
    num = int(float(len(b))*float(percentile))
    if percentile == 0:
            num = 1
    for k in range(num):
        total += float(b[k])
    results[i] = float(total) / ( float(num) + 1e-200 )

means = {}
stdev = {}

output = open(argv[1]+".boot1", 'w', -1)
output.write("#Average of bottom 1% by RMS\n#Protein |  Standard_Deviation | Real Mean\n")

for i in names:
    total = 0
    for b in bresults[i]:
        total += float(b)

    means[i] = float(total)/(float(len(bresults[i])))

    total = 0
    for b in bresults[i]:
        total += float(pow((float(b) - float(means[i])), 2))

    stdev[i] = math.sqrt(float(total)/float(len(bresults[i])))

    output.write(i+" "+str(stdev[i])+" "+str(results[i])+"\n")
