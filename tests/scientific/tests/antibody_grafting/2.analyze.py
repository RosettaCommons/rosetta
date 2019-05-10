#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  cartesian_relax/2.analyze.py
## @brief this script is part of cartesian_relax scientific test
## @author Sergey Lyskov

import collections, csv, os, sys, subprocess, math
import benchmark
from benchmark.util import quality_measures as qm

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into
config = benchmark.config()

results = {}
scorefiles = []
cutoffs_fraction_dict = {}
cutoffs_rms_dict = {}
target_results = collections.defaultdict(list)
failures = {}
off_by_more_than_X = collections.defaultdict(list)

# inputs are header labels from the scorefile, for instance "total_score" and "rmsd"
# => it figures out the column numbers from there
x_label = "rmsd"
y_label = "total_score"
outfile = "result.txt"
cutoffs = "cutoffs"

# scorefiles actually are single csv rmsd files for each target
scorefiles.extend( [ f'{working_dir}/output/{t}/native_comparison.txt' for t in targets ] )

# read cutoffs
regions = subprocess.getoutput( "grep -v '#' " + cutoffs + " | awk '{print $1}'" ).splitlines()
cutoffs_fraction = subprocess.getoutput( "grep -v '#' " + cutoffs + " | awk '{print $2}'" ).splitlines()
cutoffs_rmsd = subprocess.getoutput( "grep -v '#' " + cutoffs + " | awk '{print $3}'" ).splitlines()
cutoffs_fraction = map( float, cutoffs_fraction )
cutoffs_rmsd = map( float, cutoffs_rmsd )
cutoffs_fraction_dict.update( dict( zip ( regions, cutoffs_fraction )))
cutoffs_rms_dict.update( dict( zip ( regions, cutoffs_rmsd )))

# open results output file
f = open( outfile, "w" )

# go through scorefiles of targets
for i in range(0, len(scorefiles)):

    # single line csv files, read in as dict in one go
    with open(scorefiles[i], "r") as infile:
        t = list(csv.DictReader(infile))[0]

    # track results
    for k,v in t.items():
        target_results[k.strip()].append(float(v))

# drop OCD results as these are not yet reliable
# (also are not evaluated on the hard-coded sub-A cutoff)
del target_results['ocd']

# for each region evaluate # sub angstrom
for k in target_results.keys():
    frac = sum([ x < cutoffs_rms_dict[k] for x in target_results[k] ]) / float(len(target_results[k]))
    f.write("{} sub-{} fraction: {}, failures: ".format(k, cutoffs_rms_dict[k], frac))
    if frac < cutoffs_fraction_dict[k]: # region has failed test
        failures[k] = frac
    for i in range( 0, len(targets) ):
        if target_results[k][i] >= 1.0:
            f.write("{},".format(targets[i]))
            off_by_more_than_X[k].append(targets[i])
    f.write("\n")

# convert defaultdict to dict so it can be saved in JSON format
target_results = dict(target_results)
off_by_more_than_X = dict(off_by_more_than_X)
benchmark.save_variables('targets working_dir testname results scorefiles target_results cutoffs_rms_dict cutoffs_fraction_dict failures off_by_more_than_X')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
