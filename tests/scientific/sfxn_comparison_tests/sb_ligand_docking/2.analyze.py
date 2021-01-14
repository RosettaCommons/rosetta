#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  ligand_docking/2.analyze.py
## @brief this script is part of the ligand docking scientific test
## @author Sergey Lyskov

import os, sys, subprocess, math
import numpy as np
import benchmark
from benchmark.util import quality_measures as qm

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into 	s
config = benchmark.config()

#=======================================
def get_scoring_percentage( x, y, rmsd_cutoff, filehandle ):
    out = "score percent"
    filehandle.write(out + "\t")
    ## check if we sampled ANY near-native
    if min(x) >= rmsd_cutoff:
        score_percent = 100
        ## check at which scoring percentile we find a sub2A structure
    else:
        for j in sorted(y):
            if (x[y.index(j)]) < rmsd_cutoff:
                score_percent = round((y.index(j)) / len(y) * 100, 2)
                break
    filehandle.write(str(score_percent) + "\n")
    return {out : score_percent}
#=======================================

results = {}
scorefiles = []
# logfiles = []
cutoffs_rmsd_dict = {}
cutoffs_score_dict = {}
sampling_failures = {}
scoring_failures = {}

# inputs are header labels from the scorefile, for instance "interface_delta_X" and "ligand_rms_no_super_X"
# => it figures out the column numbers from there
x_label = "ligand_rms_no_super_X"
y_label = "interface_delta_X"
outfile = "result.txt"
cutoffs = "cutoffs"

for target in targets:
    for sfxn in sfxns:
        scorefiles.extend([f'{working_dir}/output/{sfxn}/{target}/{target}.sc'])
# logfiles.extend( [ f'{working_dir}/hpc-logs/hpc.{testname}-{sfxn}-{t}.*.log' for t in targets ] )

# read cutoffs
protein = subprocess.getoutput("grep -v '#' " + cutoffs + " | awk '{print $1}'").splitlines()
cutoffs_rmsd = subprocess.getoutput("grep -v '#' " + cutoffs + " | awk '{print $2}'").splitlines()
cutoffs_score = subprocess.getoutput("grep -v '#' " + cutoffs + " | awk '{print $3}'").splitlines()
cutoffs_rmsd = map(float, cutoffs_rmsd)
cutoffs_score = map(float, cutoffs_score)
cutoffs_rmsd_dict.update(dict(zip(protein, cutoffs_rmsd)))
cutoffs_score_dict.update(dict(zip(protein, cutoffs_score)))

# open results output file
f = open(outfile, "w")

# go through scorefiles of targets

sampling_failures[sfxn] = {}
scoring_failures[sfxn] = {}
for sfxn in sfxns:
    results[sfxn] = {}
    sfxn_sampling_failures = []
    sfxn_scoring_failures = []
    for target in targets:
        scorefile = f'{working_dir}/output/{sfxn}/{target}/{target}.sc'
        results[sfxn][target] = {}
        # get column numbers from labels, 1-indexed
        x_index = str(subprocess.getoutput("grep " + x_label + " " + scorefile).split().index(x_label) + 1)
        y_index = str(subprocess.getoutput("grep " + y_label + " " + scorefile).split().index(y_label) + 1)
        # read in score file, scores are sorted, first one is lowest
        x = subprocess.getoutput("grep -v SEQUENCE " + scorefile + " | grep -v " + y_label + " | sort -nk2 | awk '{print $" + x_index + "}'").splitlines()
        y = subprocess.getoutput("grep -v SEQUENCE " + scorefile + " | grep -v " + y_label + " | sort -nk2 | awk '{print $" + y_index + "}'").splitlines()
        # map values to floats (were strings)
        x = list(map(float, x))
        y = list(map(float, y))

        # when do we see sub2A output
        score_percent = get_scoring_percentage(x,y,int(2),f)
        results[sfxn][target]['score percent'] = score_percent

        # check for RMSD range
        f.write(target + "\t")
        val_rms = qm.check_range(x, "ligand_rms_no_super_X", f)
        results[sfxn][target]['rmsd data'] = val_rms

        # check for score range
        f.write(target + "\t")
        val_score = qm.check_range(y, "interface_delta_X", f)
        results[sfxn][target]['score data'] = val_score
        f.write("\n")

        if results[sfxn][target]['score percent']['score percent'] == 100:
            sfxn_sampling_failures.append(target)
        elif results[sfxn][target]['score percent']['score percent'] > 10 and results[sfxn][target]['score percent'][
            'score percent'] < 100:
            sfxn_scoring_failures.append(target)
        else:
            continue

    sampling_failures.update({sfxn : sfxn_sampling_failures})
    scoring_failures.update({sfxn : sfxn_scoring_failures})
f.close()
benchmark.save_variables(
    'targets nstruct working_dir testname results scorefiles cutoffs_rmsd_dict cutoffs_score_dict sampling_failures scoring_failures sfxns')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
