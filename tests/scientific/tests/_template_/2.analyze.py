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

import os, sys, subprocess, math
import numpy as np
import benchmark
from benchmark.util import quality_measures as qm

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into 	s
config = benchmark.config()

results = {}
scorefiles = []
#logfiles = []
cutoffs_rmsd_dict = {}
cutoffs_score_dict = {}
failures = []

# inputs are header labels from the scorefile, for instance "total_score" and "rmsd"
# => it figures out the column numbers from there
x_label = "rmsd"
y_label = "total_score"
outfile = "result.txt"
cutoffs = "cutoffs"

# scorefiles and logfiles
scorefiles.extend( [ f'{working_dir}/output/{t}/{t}.score' for t in targets ] )
#logfiles.extend( [ f'{working_dir}/hpc-logs/hpc.{testname}-{t}.*.log' for t in targets ] )

# get column numbers from labels, 1-indexed
x_index = str( subprocess.getoutput( "grep " + x_label + " " + scorefiles[0] ).split().index( x_label ) + 1 )
y_index = str( subprocess.getoutput( "grep " + y_label + " " + scorefiles[0] ).split().index( y_label ) + 1 )

# read cutoffs
protein = subprocess.getoutput( "grep -v '#' " + cutoffs + " | awk '{print $1}'" ).splitlines()
cutoffs_rmsd = subprocess.getoutput( "grep -v '#' " + cutoffs + " | awk '{print $2}'" ).splitlines()
cutoffs_score = subprocess.getoutput( "grep -v '#' " + cutoffs + " | awk '{print $3}'" ).splitlines()
cutoffs_rmsd = map( float, cutoffs_rmsd )
cutoffs_score = map( float, cutoffs_score )
cutoffs_rmsd_dict.update( dict( zip ( protein, cutoffs_rmsd )))
cutoffs_score_dict.update( dict( zip ( protein, cutoffs_score )))

# open results output file
f = open( outfile, "w" )

# go through scorefiles of targets
for i in range( 0, len( scorefiles ) ):

	target_results = {}

	# read in score file, scores are sorted, first one is lowest
	x = subprocess.getoutput( "grep -v SEQUENCE " + scorefiles[i] + " | grep -v " + y_label + " | sort -nk2 | awk '{print $" + x_index + "}'" ).splitlines()
	y = subprocess.getoutput( "grep -v SEQUENCE " + scorefiles[i] + " | grep -v " + y_label + " | sort -nk2 | awk '{print $" + y_index + "}'" ).splitlines()

	# map values to floats (were strings)
	x = list( map( float, x ))
	y = list( map( float, y ))

	# check for RMSDs below cutoff
	f.write( targets[i] + "\t" )
	val_cutoff = qm.check_all_values_below_cutoff( x, cutoffs_rmsd_dict[targets[i]], "rmsd", f )
	target_results.update( val_cutoff )

    # add to failues
	if val_cutoff['All rmsds < cutoff'] == False:
		failures.append( targets[i] )

	# check for scores below cutoff
	f.write( targets[i] + "\t" )
	val_cutoff = qm.check_all_values_below_cutoff( y, cutoffs_score_dict[targets[i]], "score", f )
	target_results.update( val_cutoff )

    # add to failures
	if val_cutoff['All scores < cutoff'] == False:
		failures.append( targets[i] )

	# check lowest scoring model has low RMSD
#		f.write( targets[i] + "\t" )
#		val_topscoring = check_rmsd_of_topscoring( x, cutoffs_rmsd_dict[targets[i]], f )
#		target_results.update( val_topscoring )

	# check for RMSD range
	f.write( targets[i] + "\t" )
	val_rms = qm.check_range( x, "rmsd", f )
	target_results.update( val_rms )

	# check for score range
	f.write( targets[i] + "\t" )
	val_score = qm.check_range( y, "score", f )
	target_results.update( val_score )

	# check runtime
	# runtime = subprocess.getoutput( "grep \"reported success\" " + logfiles[i] + " | awk '{print $6}'" ).splitlines()
	# runtime = list( map( float, runtime ))
	# if len(runtime) > 0:
	# 	print (targets[i], "\t"),
	# 	val_runtime = check_range( runtime, "runtime" )
	# 	target_results.update( val_runtime )

	# TODO: check discrimination score
	# discrimination score is likely better for abinitio sampling + refinement
	# few options:
	# 1) ff_metric script in this directory: easy + simple but need to test the metric
	# 2) discrimination score script in https://github.com/RosettaCommons/bakerlab_scripts/tree/master/boinc/scoring_methods
	# => gives me 0, a little too complicated
	# 3) https://github.com/RosettaCommons/bakerlab_scripts/blob/master/boinc/score_energy_landscape.py
	# => requires lots of dependent scripts, even more complicated

	results.update( {targets[i] : target_results} )
	f.write( "\n" )

f.close()

benchmark.save_variables('debug targets nstruct working_dir testname results scorefiles cutoffs_rmsd_dict cutoffs_score_dict failures')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
