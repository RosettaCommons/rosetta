#!/usr/bin/env python
# encoding: utf-8

from __future__ import division

import os, sys, glob

try:
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('job')
    job_path = parser.parse_args().job
except ImportError:
    assert len(sys.argv) == 2
    job_path = sys.argv[1]

rebuild_rmsds = []
centroid_rmsds = []
fullatom_rmsds = []

subjob_pattern = os.path.join(job_path, '???')
subjobs = glob.glob(subjob_pattern)

for subjob_path in subjobs:
    log_path = os.path.join(subjob_path, 'output.txt')

    with open(log_path) as file:
        for line in file:
            if 'loop_rebuildrms' in line:
                rmsd = float(line.split()[2])
                rebuild_rmsds.append(rmsd)
            if 'loop_cenrms' in line:
                rmsd = float(line.split()[2])
                centroid_rmsds.append(rmsd)
            if 'loop_rms' in line:
                rmsd = float(line.split()[2])
                fullatom_rmsds.append(rmsd)

def mean(x):
    try: return sum(x) / len(x)
    except ZeroDivisionError: return float('NaN')

def subangstrom(x):
    count = 0

    for i in x:
        if i < 1:
            count += 1

    try: return 100 * count / len(x)
    except ZeroDivisionError: return float('NaN')


print job_path
print "Completed Jobs:        ", len(fullatom_rmsds)
print "Average Rebuild RMSD:  ", mean(rebuild_rmsds)
print "Average Centroid RMSD: ", mean(centroid_rmsds)
print "Average Full-atom RMSD:", mean(fullatom_rmsds)
print "Percent Subangstrom:    %.2f%%" % subangstrom(fullatom_rmsds)
