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
print(config)

results = {}
scorefiles = []
cutoffs_rmsd_dict = {}
cutoffs_discrim_dict = {}
failures = []

# inputs are header labels from the scorefile, for instance "total_score" and "rmsd"
# => it figures out the column numbers from there
#x_label = "H3_new" # change to H3_RMS after testing
x_label = "Irms"
y_label = "I_sc"
outfile = "result.txt"
cutoffs = "cutoffs"

# scorefiles and logfiles
scorefiles.extend( [ f'{working_dir}/output/{t}/{t}_rmsd.sc' for t in targets ] )

# get column numbers from labels, 1-indexed
x_index = str( subprocess.getoutput( "grep " + x_label + " " + scorefiles[0] ).split().index( x_label ) + 1 )
y_index = str( subprocess.getoutput( "grep " + y_label + " " + scorefiles[0] ).split().index( y_label ) + 1 )

# read cutoffs
protein = subprocess.getoutput( "grep -v '#' " + cutoffs + " | awk '{print $1}'" ).splitlines()
cutoffs_rmsd = subprocess.getoutput( "grep -v '#' " + cutoffs + " | awk '{print $2}'" ).splitlines()
cutoffs_discrim = subprocess.getoutput( "grep -v '#' " + cutoffs + " | awk '{print $3}'" ).splitlines()
cutoffs_rmsd = map( float, cutoffs_rmsd )
cutoffs_discrim = map( float, cutoffs_discrim )
cutoffs_rmsd_dict.update( dict( zip ( protein, cutoffs_rmsd )))
cutoffs_discrim_dict.update( dict( zip ( protein, cutoffs_discrim )))

# open results output file
f = open( outfile, "w" )
f.write("target\t" + "min_rms\t" + "discrim\t\n") # min rms is of the interface

# go through scorefiles of targets
for i in range( 0, len( scorefiles ) ):
	target_results = {}

	# read in score file, scores are sorted, first one is lowest
	x = subprocess.getoutput( "grep -v SEQUENCE " + scorefiles[i] + " | grep -v " + y_label + " | sort -nk2 | awk '{print $" + x_index + "}'" ).splitlines()
	y = subprocess.getoutput( "grep -v SEQUENCE " + scorefiles[i] + " | grep -v " + y_label + " | sort -nk2 | awk '{print $" + y_index + "}'" ).splitlines()

	# map values to floats (were strings)
	x = list( map( float, x )) # x is Irms
	y = list( map( float, y )) # y is I_sc

	# check for minimum RMSD below cutoff
	target_results["min_rms"] = min(x)
	f.write( targets[i] + "\t" + str(target_results["min_rms"]) + "\t" )

    # add to failues
	if target_results["min_rms"] > cutoffs_rmsd_dict[targets[i]]:
		failures.append( targets[i] )

	# check for discrimination scores below cutoff -- can't do this in debug mode (really)
	if config['debug']:
		target_results["discrim"] = 0.0
		print("debuggin!")
	else: # [1, 2, 3, 4, 6] are the offsets from the min(rmsd) for discrim calculation
		target_results["discrim"] = qm.calc_Conway_discrim_score( x, y, [0, 1, 2, 3, 4, 6] )
	f.write( str(target_results["discrim"]) )

    # add to failures
	if target_results["discrim"] > cutoffs_discrim_dict[targets[i]]:
		failures.append( targets[i] )

	results.update( {targets[i] : target_results} )
	f.write( "\n" )

f.close()

benchmark.save_variables('debug targets nstruct working_dir testname results scorefiles cutoffs_discrim_dict cutoffs_rmsd_dict failures')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
