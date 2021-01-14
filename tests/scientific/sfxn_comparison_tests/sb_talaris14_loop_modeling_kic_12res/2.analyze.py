#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

# @file  cartesian_relax/2.analyze.py
# @brief this script is part of cartesian_relax scientific test
# @author Sergey Lyskov
# @author Phuong T. Nguyen (tranphuonguns@gmail.com)

import os
import sys
import subprocess
import math
import numpy as np
import benchmark
from benchmark.util import quality_measures as qm

# Python black magic: load all variables saved by previous script into  s
benchmark.load_variables()
config = benchmark.config()

results = {}
scorefiles = []
#logfiles = []
cutoffs_rmsd_dict = {}
failures = []


def check_any_topNpercentscores_below_rmsdcutoff(score_rmsd_dict, cutoff, tag, filehandle, percentage):
    out = "Any " + tag + " < cutoff"
    filehandle.write(out + " " + str(cutoff) + "\t")
    # sort
    col = [score_rmsd_dict[key] for key in sorted(score_rmsd_dict.keys())]
    partial_col = col[:int(math.ceil(percentage * len(col) / 100.0))]
    if any(i <= cutoff for i in partial_col):
        value = True
    else:
        value = False
    filehandle.write(str(value) + "\n")
    return {out: value}


# inputs are header labels from the scorefile, for instance "total_score" and "rmsd"
# => it figures out the column numbers from there
x_label = "looprms"
y_label = "total_score"
outfile = "result.txt"
cutoffs = "cutoffs"

# scorefiles and logfiles
scorefiles.extend([f'{working_dir}/output/{t}/{t}.score' for t in targets])
#logfiles.extend( [ f'{working_dir}/hpc-logs/hpc.{testname}-{t}.*.log' for t in targets ] )

# get column numbers from labels, 1-indexed
x_index = str(subprocess.getoutput("grep " + x_label + " " + scorefiles[0]).split().index(x_label) + 1)
y_index = str(subprocess.getoutput("grep " + y_label + " " + scorefiles[0]).split().index(y_label) + 1)

# read cutoffs
protein = subprocess.getoutput("grep -v '#' " + cutoffs + " | awk '{print $1}'").splitlines()
cutoffs_rmsd = subprocess.getoutput("grep -v '#' " + cutoffs + " | awk '{print $2}'").splitlines()
cutoffs_rmsd = map(float, cutoffs_rmsd)
cutoffs_rmsd_dict.update(dict(zip(protein, cutoffs_rmsd)))

# open results output file
f = open(outfile, "w")

# go through scorefiles of targets
for i in range(0, len(scorefiles)):
    target_results = {}

    # read in score file, scores are sorted, first one is lowest
    x = subprocess.getoutput("grep -v SEQUENCE " + scorefiles[i] + " | grep -v " + y_label + " | sort -nk2 | awk '{print $" + x_index + "}'").splitlines()
    y = subprocess.getoutput("grep -v SEQUENCE " + scorefiles[i] + " | grep -v " + y_label + " | sort -nk2 | awk '{print $" + y_index + "}'").splitlines()
    x_nonsorted = subprocess.getoutput("grep -v SEQUENCE " + scorefiles[i] + " | grep -v " + y_label + " | awk '{print $" + x_index + "}'").splitlines()
    y_nonsorted = subprocess.getoutput("grep -v SEQUENCE " + scorefiles[i] + " | grep -v " + y_label + " | awk '{print $" + y_index + "}'").splitlines()

    # map values to floats (were strings)
    x = list(map(float, x))
    y = list(map(float, y))
    x_nonsorted = list(map(float, x_nonsorted))
    y_nonsorted = list(map(float, y_nonsorted))

    # make a dictionary of score vs rmsd, scores as keys
    score_rmsd = dict(zip(y_nonsorted, x_nonsorted))

    # check if any top 10% scores below rmsd cutoff
    f.write(targets[i] + "\t")
    val_topscores_rmsdcutoff = check_any_topNpercentscores_below_rmsdcutoff(score_rmsd, cutoffs_rmsd_dict[targets[i]], "top10% scores", f, 10)

    # add to failues
    if val_topscores_rmsdcutoff['Any top10% scores < cutoff'] == False:
        failures.append(targets[i])

    # check for RMSD range
    f.write(targets[i] + "\t")
    val_rms = qm.check_range(x, "rmsd", f)
    target_results.update(val_rms)

    # check for score range
    f.write(targets[i] + "\t")
    val_score = qm.check_range(y, "score", f)
    target_results.update(val_score)

    # check runtime
    # runtime = subprocess.getoutput( "grep \"reported success\" " + logfiles[i] + " | awk '{print $6}'" ).splitlines()
    # runtime = list( map( float, runtime ))
    # if len(runtime) > 0:
    #   print (targets[i], "\t"),
    #   val_runtime = check_range( runtime, "runtime" )
    #   target_results.update( val_runtime )

    # TODO: check discrimination score
    # discrimination score is likely better for abinitio sampling + refinement
    # few options:
    # 1) ff_metric script in this directory: easy + simple but need to test the metric
    # 2) discrimination score script in https://github.com/RosettaCommons/bakerlab_scripts/tree/master/boinc/scoring_methods
    # => gives me 0, a little too complicated
    # 3) https://github.com/RosettaCommons/bakerlab_scripts/blob/master/boinc/score_energy_landscape.py
    # => requires lots of dependent scripts, even more complicated

    results.update({targets[i]: target_results})
    f.write("\n")

f.close()

# Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
benchmark.save_variables(
    'targets nstruct working_dir testname results scorefiles cutoffs_rmsd_dict failures')
