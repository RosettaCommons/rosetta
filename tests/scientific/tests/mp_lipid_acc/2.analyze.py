#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  cartesian_relax/2.analyze.py
## @brief this script is part of the mp_lipid_acc scientific test
## @author JKLeman

import os, sys, subprocess, math
import numpy as np
import benchmark
from benchmark.util import quality_measures as qm

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into 	s
config = benchmark.config()

###################################
def accuracies_from_bfactors( protein, file_real, file_pred, filehandle ):
	
	# initialize variables
	tp = 0
	tn = 0
	fp = 0
	fn = 0
	acc = 0
	sens = 0
	spec = 0
    
	# checking for file existence
	if not os.path.exists( file_real ): 
	   print ("File " + file_real + " doesn't exist.")
	   return {"tp": 0, "tn": 0, "fp":0, "fn": 0, "acc": 0, "sens": 0, "spec": 0}
	   
	if not os.path.exists( file_pred ): 
	   print ("File " + file_pred + " doesn't exist.")
	   return {"tp": 0, "tn": 0, "fp":0, "fn": 0, "acc": 0, "sens": 0, "spec": 0}
	
	# get manually curated lipid acc
	real = (subprocess.getoutput("grep CA " + file_real + " | grep ATOM | awk '{print substr($0,62,2)}'") ).splitlines()
	resn_real = (subprocess.getoutput("grep CA " + file_real + " | grep ATOM | awk '{print substr($0,23,4)}'") ).splitlines()
	resn_real = list(map( int, resn_real ))

	# get predicted lipid acc
	pred = (subprocess.getoutput("grep CA " + file_pred + " | grep ATOM | awk '{print substr($0,62,2)}'") ).splitlines()
	resn_pred = (subprocess.getoutput("grep CA " + file_pred + " | grep ATOM | awk '{print substr($0,23,4)}'") ).splitlines()
	resn_pred = list(map( int, resn_pred ))

	# initialize empty arrays
	real_array = []
	pred_array = []
	real_array = [0 for i in range( 0, len( real )+1 )]
	pred_array = [0 for i in range( 0, len( real )+1 )]

	# order by residue number (this requires pose numbering)
	for i in range( 0, len( real ) ):

		real_array[ resn_real[ i ] ] = real[ i ]  
		pred_array[ resn_pred[ i ] ] = pred[ i ]  
	
	# compare real and prediction
	for i in range( 0, len( real ) ):
		
		if real_array[ i ] == "50" and pred_array[ i ] == "50": 
			tp += 1
		elif real_array[ i ] == " 0" and pred_array[ i ] == " 0":
			tn += 1
		elif real_array[ i ] == " 0" and pred_array[ i ] == "50":
			fp += 1
		elif real_array[ i ] == "50" and pred_array[ i ] == " 0":
			fn += 1
	
	# get accuracy, sensitivity and specificity
	acc = round( float(tp + tn) / float(tp + tn + fp + fn) * 100, 2 )
	sens = round( float(tp) / float(tp + fn) * 100, 2 )
	spec = round( float(tn) / float(tn + fp) * 100, 2 )

	# print results for protein
	filehandle.write(protein + "\t" + str(tp) + "\t" + str(tn) + "\t" + str(fp) + "\t" + str(fn) + "\t" + str(acc) + "\t" + str(sens) + "\t" + str(spec) + "\n")
	
	# return results in a dict
	return {"tp": tp, "tn": tn, "fp":fp, "fn": fn, "acc": acc, "sens": sens, "spec": spec}

###################################

rosetta_dir = config['rosetta_dir']

results = {}
cutoffs_rmsd_dict = {}
cutoffs_score_dict = {}
failures = []

# get filenames
pdb_truth = [ f'{rosetta_dir}/tests/scientific/data/{testname}/db_hand_curated/{t}__tr_0001.pdb' for t in targets] 
pdb_pred = [ f'{working_dir}/output/{t}__tr_0001.pdb' for t in targets ]

# open results output file
outfile = "result.txt"
f = open( outfile, "w" )
f.write( "protein\ttp\ttn\tfp\tfn\tacc\tsens\tspec\n" )

# inialize variables for averages
#	quality_measures = {}

# targets are already sorted by decreasing accuracies
targets_sorted = targets
    
acc_before = 100

# go through targets, get accuracies
for i in range(0, len( targets_sorted )):

	target = targets_sorted[i]
	target_results = {}
	
	# get accuracies and append to dict
	print (target, pdb_truth[i], pdb_pred[i])
	acc = accuracies_from_bfactors( target, pdb_truth[i], pdb_pred[i], f )
#		quality_measures[target] = acc['acc']
	target_results.update( acc )

	accuracy = acc['acc']

	# if new accuracies are not sorted any more, call this a failure
	# i.e. if accuracy is larger than previous one in the list (i.e. they 
	# are not sorted by decreasing values any more)
	if accuracy > acc_before:
		failures.append( target )
	
	# update overall results
	results.update( {targets[i] : target_results} )
	
	acc_before = accuracy

f.close()

benchmark.save_variables('targets nstruct working_dir testname results scorefiles cutoffs_rmsd_dict cutoffs_score_dict failures')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
