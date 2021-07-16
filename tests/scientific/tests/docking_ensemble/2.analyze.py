#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  docking_ensemble/2.analyze.py
## @brief this script is part of docking ensemble scientific test
## @author Sergey Lyskov
## @author Ameya Harmalkar

import os, sys, subprocess, math
import numpy as np
import benchmark
from benchmark.util import quality_measures as qm
from benchmark.util import bootstrap as bs

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into 	s
config = benchmark.config()

results = {}
scorefiles = []
#logfiles = []
cutoffs_rank_dict = {}
cutoffs_rmsd_dict = {}
failures = []

# inputs are header labels from the scorefile, for instance "total_score" and "rmsd"
# => it figures out the column numbers from there
x_label = "Irms"
y_label = "I_sc"
z_label = "CAPRI_rank"
outfile = "result.txt"
cutoffs = "cutoffs"

# scorefiles and logfiles
scorefiles.extend( [ f'{working_dir}/output/{t}/{t}.score' for t in targets ] )
#logfiles.extend( [ f'{working_dir}/hpc-logs/hpc.{testname}-{t}.*.log' for t in targets ] )

# get column numbers from labels, 1-indexed
x_index = str( subprocess.getoutput( "grep " + x_label + " " + scorefiles[0] ).split().index( x_label ) + 1 )
y_index = str( subprocess.getoutput( "grep " + y_label + " " + scorefiles[0] ).split().index( y_label ) + 1 )
z_index = str( subprocess.getoutput( "grep " + z_label + " " + scorefiles[0] ).split().index( z_label ) + 1 )

# get number of fields in scorefile
nfields = len( subprocess.getoutput( "grep -v SEQUENCE " + scorefiles[0] + " | grep " + y_label + " | head -n1" ).split())

# read cutoffs
protein = subprocess.getoutput( "grep -v '#' " + cutoffs + " | awk '{print $1}'" ).splitlines()
cutoffs_rank = subprocess.getoutput( "grep -v '#' " + cutoffs + " | awk '{print $2}'" ).splitlines()
cutoffs_rmsd = subprocess.getoutput( "grep -v '#' " + cutoffs + " | awk '{print $3}'" ).splitlines()
cutoffs_rank = map( float, cutoffs_rank )
cutoffs_rmsd = map( float, cutoffs_rmsd )
cutoffs_rank_dict.update( dict( zip ( protein, cutoffs_rank )))
cutoffs_rmsd_dict.update( dict( zip ( protein, cutoffs_rmsd )))

# open results output file
f = open( outfile, "w" )

# go through scorefiles of targets
for i in range( 0, len( scorefiles ) ):

	target_results = {}

	# read in score file, scores are sorted by I_sc, first one is lowest
	x = subprocess.getoutput( "awk '{if(NF==" + str(nfields) + ") print}' " + scorefiles[i] + " | grep -v SEQUENCE | grep -v " + y_label + " |  sort -nk6 | awk '{print $" + x_index + "}'" ).splitlines()
	y = subprocess.getoutput( "awk '{if(NF==" + str(nfields) + ") print}' " + scorefiles[i] + " | grep -v SEQUENCE | grep -v " + y_label + " |  sort -nk6 | awk '{print $" + y_index + "}'" ).splitlines()
	z = subprocess.getoutput( "awk '{if(NF==" + str(nfields) + ") print}' " + scorefiles[i] + " | grep -v SEQUENCE | grep -v " + y_label + " |  sort -nk6 | awk '{print $" + z_index + "}'" ).splitlines()

	# sort the randomly re-sampled subset of models on y_label, taking the top 5
	# then, check if the 5-top-scoring models have metric z_label >= 1 (True means use >= for the metric)
        # take len(y) random models with replacement from the original set of models
        # repeat the process of checking if the top-5 models have a CAPRI rank of 1 or better 1000 times
	nscore, nstd = bs.bootstrap_NX( scorefiles[i], y_label, 5, z_label, 1, True, sample_size=len(y), n_samples=1000)

	# map values to floats (were strings)
	x = list( map( float, x ))
	y = list( map( float, y ))
	z = list( map( float, z ))

	# check for CAPRI bins
	# bin of highest CAPRI score should be better (>= rank cutoff)
	rank = True
	if max(z) < cutoffs_rank_dict[targets[i]]:
		failures.append( targets[i] )
		rank = False
	f.write( targets[i] + "\t" + "highest CAPRI rank >= cutoff:\t" + str(rank) + "\n" )
		
	# check lowest Isc scoring model has low RMSD
	f.write( targets[i] + "\t" )
	val_topscoring = qm.check_rmsd_of_topscoring( x, cutoffs_rmsd_dict[targets[i]], f )
	target_results.update( val_topscoring )
	
	if val_topscoring["RMSD <= cutoff"] == False and targets[i] not in failures:
		failures.append( targets[i] )
	
	#check for N5 above cutoff
	f.write( targets[i] + "\t" )
	val_cutoff = qm.check_avgNX_above_cutoff(nscore, 3.0, "N5", f)
	target_results.update( val_cutoff )

	# check for N5 range
	f.write( targets[i] + "\t" )
	val_N5 = qm.get_N5_spread(nscore, nstd, "N5", f)
	target_results.update( val_N5 )

	# check for RMSD range
	f.write( targets[i] + "\t" )
	val_rms = qm.check_range( x, "Irms", f )
	target_results.update( val_rms )

	# check for score range
	f.write( targets[i] + "\t" )
	val_score = qm.check_range( y, "I_sc", f )
	target_results.update( val_score )

	results.update( {targets[i] : target_results} )
	f.write( "\n" )

f.close()

benchmark.save_variables('debug targets nstruct working_dir testname results scorefiles cutoffs_rmsd_dict cutoffs_score_dict failures')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
