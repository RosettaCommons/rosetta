#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  mhc_epitope_energy/2.analyze.py
## @brief this script is part of mhc_epitope_energy scientific test
## @author Sergey Lyskov
## @author Brahm Yachnin (brahm.yachnin@rutgers.edu)

import os, sys, subprocess, math
import numpy as np
import benchmark
from benchmark.util import quality_measures as qm

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into 	s
config = benchmark.config()

results = {}
scorefiles = []
cutoffs_mhc_epitope_dict = {}
cutoffs_delta_mhc_epitope_dict = {}
cutoffs_mhc_pct_drop_dict = {}
failures = []
failures_dict = {}

# inputs are header labels from the scorefile, for instance "total_score" and "rmsd"
# => it figures out the column numbers from there
# We are going to test against the following scoreterms
label_mhc = "mhc_epitope"
label_delta_mhc = "delta_mhc_epitope"
label_base_total_score = "base_total_energy"
label_seqrec = "seqrec_seqrec"
label_coreseqrec = "seqrec_core_seqrec"
label_delta_packstat = "delta_packstat"
label_delta_buns = "delta_buried_unsat"
label_delta_netcharge = "delta_netcharge"

#One output filename per subtest
outfile_mhc_epitope = "result_mhc_epitope.txt"
outfile_delta_mhc_epitope = "result_delta_mhc_epitope.txt"
outfile_totalscore = "result_base_total_score_vs_mhc.txt"
outfile_seqrec = "result_sequence_rec_vs_mhc.txt"
outfile_coreseqrec = "result_core_sequence_rec_vs_mhc.txt"
outfile_delta_packstat = "result_delta_packstat_vs_mhc.txt"
outfile_delta_buns = "result_delta_buried_unsat_vs_mhc.txt"
outfile_delta_netcharge = "result_delta_netcharge_vs_mhc.txt"

cutoffs = "cutoffs"

# scorefiles
scorefiles.extend( [ f'{working_dir}/output/{t}/{t}.score' for t in targets ] )

# get column numbers from labels, 1-indexed
# Following label_ variables above
index_mhc = str( subprocess.getoutput( "grep " + label_mhc + " " + scorefiles[0] ).split().index( label_mhc ) + 1 )
index_delta_mhc = str( subprocess.getoutput( "grep " + label_delta_mhc + " " + scorefiles[0] ).split().index( label_delta_mhc ) + 1 )
index_base_total_score = str( subprocess.getoutput( "grep " + label_base_total_score + " " + scorefiles[0] ).split().index( label_base_total_score ) + 1 )
index_seqrec = str( subprocess.getoutput( "grep " + label_seqrec + " " + scorefiles[0] ).split().index( label_seqrec ) + 1 )
index_coreseqrec = str( subprocess.getoutput( "grep " + label_coreseqrec + " " + scorefiles[0] ).split().index( label_coreseqrec ) + 1 )
index_delta_packstat = str( subprocess.getoutput( "grep " + label_delta_packstat + " " + scorefiles[0] ).split().index( label_delta_packstat ) + 1 )
index_delta_buns = str( subprocess.getoutput( "grep " + label_delta_buns + " " + scorefiles[0] ).split().index( label_delta_buns ) + 1 )
index_delta_netcharge = str( subprocess.getoutput( "grep " + label_delta_netcharge + " " + scorefiles[0] ).split().index( label_delta_netcharge ) + 1 )

# read cutoffs
protein = subprocess.getoutput( "grep -v '#' " + cutoffs + " | awk '{print $1}'" ).splitlines()
cutoffs_mhc_epitope = subprocess.getoutput( "grep -v '#' " + cutoffs + " | awk '{print $2}'" ).splitlines()
cutoffs_delta_mhc_epitope = subprocess.getoutput( "grep -v '#' " + cutoffs + " | awk '{print $3}'" ).splitlines()
cutoffs_mhc_pct_drop = subprocess.getoutput( "grep -v '#' " + cutoffs + " | awk '{print $4}'" ).splitlines()

cutoffs_mhc_epitope = map( float, cutoffs_mhc_epitope )
cutoffs_delta_mhc_epitope = map( float, cutoffs_delta_mhc_epitope )
cutoffs_mhc_pct_drop = map( float, cutoffs_mhc_pct_drop )

cutoffs_mhc_epitope_dict.update( dict( zip ( protein, cutoffs_mhc_epitope )))
cutoffs_delta_mhc_epitope_dict.update( dict( zip ( protein, cutoffs_delta_mhc_epitope )))
cutoffs_mhc_pct_drop_dict.update( dict( zip ( protein, cutoffs_mhc_pct_drop )))

# open results output file -- one per subtest
f_mhc_epitope = open( outfile_mhc_epitope, "w" )
f_delta_mhc_epitope = open( outfile_delta_mhc_epitope, "w" )
f_totalscore = open( outfile_totalscore, "w" )
f_seqrec = open( outfile_seqrec, "w" )
f_coreseqrec = open( outfile_coreseqrec, "w" )
f_delta_packstat = open( outfile_delta_packstat, "w" )
f_delta_buns = open( outfile_delta_buns, "w" )
f_delta_netcharge = open( outfile_delta_netcharge, "w" )

# go through scorefiles of targets
for i in range( 0, len( scorefiles ) ):
	target_results = {}

	# read in score file, scores are sorted by total_score, first one is lowest
	# grep -v label_mhc is to remove the label line
	# Following label_ variables above
	mhc = subprocess.getoutput( "grep -v SEQUENCE " + scorefiles[i] + " | grep -v " + label_mhc + " | sort -nk2 | awk '{print $" + index_mhc + "}'" ).splitlines()
	delta_mhc = subprocess.getoutput( "grep -v SEQUENCE " + scorefiles[i] + " | grep -v " + label_mhc + " | sort -nk2 | awk '{print $" + index_delta_mhc + "}'" ).splitlines()
	base_total_score = subprocess.getoutput( "grep -v SEQUENCE " + scorefiles[i] + " | grep -v " + label_mhc + " | sort -nk2 | awk '{print $" + index_base_total_score + "}'" ).splitlines()
	seqrec = subprocess.getoutput( "grep -v SEQUENCE " + scorefiles[i] + " | grep -v " + label_mhc + " | sort -nk2 | awk '{print $" + index_seqrec + "}'" ).splitlines()
	coreseqrec = subprocess.getoutput( "grep -v SEQUENCE " + scorefiles[i] + " | grep -v " + label_mhc + " | sort -nk2 | awk '{print $" + index_coreseqrec + "}'" ).splitlines()
	delta_packstat = subprocess.getoutput( "grep -v SEQUENCE " + scorefiles[i] + " | grep -v " + label_mhc + " | sort -nk2 | awk '{print $" + index_delta_packstat + "}'" ).splitlines()
	delta_buns = subprocess.getoutput( "grep -v SEQUENCE " + scorefiles[i] + " | grep -v " + label_mhc + " | sort -nk2 | awk '{print $" + index_delta_buns + "}'" ).splitlines()
	delta_netcharge = subprocess.getoutput( "grep -v SEQUENCE " + scorefiles[i] + " | grep -v " + label_mhc + " | sort -nk2 | awk '{print $" + index_delta_netcharge + "}'" ).splitlines()

	# map values to floats (were strings)
	mhc = list( map( float, mhc ))
	delta_mhc = list( map( float, delta_mhc ))
	base_total_score = list( map( float, base_total_score ))
	seqrec = list( map( float, seqrec ))
	coreseqrec = list( map( float, coreseqrec ))
	delta_packstat = list( map( float, delta_packstat ))
	delta_buns = list( map( float, delta_buns ))
	delta_netcharge = list( map( float, delta_netcharge ))
	
	# Calculate pct_drop_mhc, the normalized drop in mhc_score
	# mhc - delta_mhc should give the native mhc_epitope score
	pct_drop_mhc = []
	for idx in range( 0, len(mhc) ):
		pct_drop_mhc.append( 100 * (mhc[idx] - (mhc[idx] - delta_mhc[idx]) ) / (mhc[idx] - delta_mhc[idx]) )

	###### Go through subtests using data from above here
	
	###Subtest 1: mhc_epitope cutoff
	# check that 90% of mhc_epitope scores are less than the cutoff for that pdb
	f_mhc_epitope.write( targets[i] + "\t" )
	result = qm.check_xpercent_values_below_cutoff( mhc, cutoffs_mhc_epitope_dict[targets[i]], "mhc_epitope", f_mhc_epitope, 90 )
	target_results.update( result )

    # add to failures if there are no decoys within both cutoffs
	if result['All mhc_epitopes < cutoff'] == False:
		failures.append( targets[i] )
		failures_dict.setdefault(targets[i], []).append('mhc_epitope')
		
	###Subtest 2: delta_mhc_epitope cutoff
	# check that 90% of  delta_mhc_epitope scores are less than the cutoff for that pdb
	f_delta_mhc_epitope.write( targets[i] + "\t" )
	result = qm.check_xpercent_values_below_cutoff( delta_mhc, cutoffs_delta_mhc_epitope_dict[targets[i]], "delta_mhc_epitope", f_delta_mhc_epitope, 90 )
	target_results.update( result )

    # add to failures if there are no decoys within both cutoffs
	if result['All delta_mhc_epitopes < cutoff'] == False:
		failures.append( targets[i] )
		failures_dict.setdefault(targets[i], []).append('delta_mhc_epitope')
	
	###Subtest 3: base_total_score vs. mhc_epitope
	# check for at least one decoy in top 50% of base_total_score and has decoys below the mhc_epitope percent drop cutoff
	f_totalscore.write( targets[i] + "\t" )
	result = qm.check_for_2d_top_values( pct_drop_mhc, base_total_score, "pct_drop_mhc", "base_total_score", cutoffs_mhc_pct_drop_dict[targets[i]], 50, f_totalscore, yabsolute = False)
	target_results.update( result )

    # add to failures if there are no decoys within both cutoffs (reading from position 0 of the result tuple)
	if result['Num decoys with pct_drop_mhc and base_total_score better than cutoffs'][0] == 0:
		failures.append( targets[i] )
		failures_dict.setdefault(targets[i], []).append('base_total_score_vs_mhc_epitope')
		
	###Subtest 4: sequence_recovery vs. mhc_epitope
	# check for at least one third of decoys with 80% sequence recovery and has decoys below the mhc_epitope percent drop cutoff
	f_seqrec.write( targets[i] + "\t" )
	result = qm.check_for_2d_top_values( pct_drop_mhc, seqrec, "pct_drop_mhc", "sequence_recovery", cutoffs_mhc_pct_drop_dict[targets[i]], 0.80, f_seqrec, yminimize = False)
	target_results.update( result )

    # add to failures if there are not enough decoys within both cutoffs (reading from position 0 of the result tuple)
	if result['Num decoys with pct_drop_mhc and sequence_recovery better than cutoffs'][0] < nstruct/3:
		failures.append( targets[i] )
		failures_dict.setdefault(targets[i], []).append('sequence_recovery_vs_mhc_epitope')

	results.update( {targets[i] : target_results} )
	
	###Subtest 5: core_sequence_recovery vs. mhc_epitope
	# check for at least one third of decoys with 87% core sequence recovery and has decoys below the mhc_epitope percent drop cutoff
	f_coreseqrec.write( targets[i] + "\t" )
	result = qm.check_for_2d_top_values( pct_drop_mhc, coreseqrec, "pct_drop_mhc", "core_sequence_recovery", cutoffs_mhc_pct_drop_dict[targets[i]], 0.87, f_coreseqrec, yminimize = False)
	target_results.update( result )

    # add to failures if there are not enough decoys within both cutoffs (reading from position 0 of the result tuple)
	if result['Num decoys with pct_drop_mhc and core_sequence_recovery better than cutoffs'][0] < nstruct/3:
		failures.append( targets[i] )
		failures_dict.setdefault(targets[i], []).append('core_sequence_recovery_vs_mhc_epitope')

	results.update( {targets[i] : target_results} )
	
	###Subtest 6: delta_packstat vs. mhc_epitope
	# check for at least one decoy with delta_packstat greater than 0 and has decoys below the mhc_epitope percent drop cutoff
	f_delta_packstat.write( targets[i] + "\t" )
	result = qm.check_for_2d_top_values( pct_drop_mhc, delta_packstat, "pct_drop_mhc", "delta_packstat", cutoffs_mhc_pct_drop_dict[targets[i]], 0, f_delta_packstat, yminimize = False)
	target_results.update( result )

    # add to failures if there are no decoys within both cutoffs (reading from position 0 of the result tuple)
	if result['Num decoys with pct_drop_mhc and delta_packstat better than cutoffs'][0] == 0:
		failures.append( targets[i] )
		failures_dict.setdefault(targets[i], []).append('delta_packstat_vs_mhc_epitope')

	results.update( {targets[i] : target_results} )
	
	###Subtest 7: delta_buried_unsat vs. mhc_epitope
	# check for at least one third of decoys with delta_buried_unsat less than 3 (no more than 3 new buried unsat) and has decoys below the mhc_epitope percent drop cutoff
	f_delta_buns.write( targets[i] + "\t" )
	result = qm.check_for_2d_top_values( pct_drop_mhc, delta_buns, "pct_drop_mhc", "delta_buried_unsat", cutoffs_mhc_pct_drop_dict[targets[i]], 3, f_delta_buns)
	target_results.update( result )

    # add to failures if there are not enough decoys within both cutoffs (reading from position 0 of the result tuple)
	if result['Num decoys with pct_drop_mhc and delta_buried_unsat better than cutoffs'][0] < nstruct/3:
		failures.append( targets[i] )
		failures_dict.setdefault(targets[i], []).append('delta_buried_unsat_vs_mhc_epitope')

	results.update( {targets[i] : target_results} )
	
	###Subtest 8: delta_netcharge vs. mhc_epitope
	# check that 5% of decoys with abs(delta_netcharge) less than 3 (no more than 3 units of charge drift) and has decoys below the mhc_epitope percent drop cutoff
	f_delta_netcharge.write( targets[i] + "\t" )
	result = qm.check_for_2d_top_values( pct_drop_mhc, list(map(abs,delta_netcharge)), "pct_drop_mhc", "abs(delta_netcharge)", cutoffs_mhc_pct_drop_dict[targets[i]], 3, f_delta_netcharge)
	target_results.update( result )

    # add to failures if there are not enough decoys within both cutoffs (reading from position 0 of the result tuple)
	if result['Num decoys with pct_drop_mhc and abs(delta_netcharge) better than cutoffs'][0] < nstruct/20:
		failures.append( targets[i] )
		failures_dict.setdefault(targets[i], []).append('delta_netcharge_vs_mhc_epitope')

	results.update( {targets[i] : target_results} )
	
	###End of subtests
	
f_totalscore.close()
f_seqrec.close()
f_coreseqrec.close()
f_delta_packstat.close()
f_delta_buns.close()
f_delta_netcharge.close()

benchmark.save_variables('targets nstruct working_dir testname results scorefiles cutoffs_mhc_epitope_dict cutoffs_delta_mhc_epitope_dict failures failures_dict')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
