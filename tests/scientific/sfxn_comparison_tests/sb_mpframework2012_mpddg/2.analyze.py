#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  mpddg/2.analyze.py
## @brief this script is part of mpddg scientific test
## @author Julia Koehler Leman, Hope Woods, Johanna Tiemann

import os, sys, subprocess, math
import numpy as np
import benchmark
from scipy.stats import linregress
from benchmark.util import quality_measures as qm

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into 	s
config = benchmark.config()

# filenames
ddg_out = "ddg.txt"
outfile = "result.txt"
cutoffs = "cutoffs"

# initializations
results = {}
cutoffs_corr_dict = {}
failures = []

# read in all data into ddg_data, a list of dicts
ddg_data = {}

# open results output file
f = open( outfile, "w" )

# write data into ddg out file
o = open( ddg_out, "w" )
o.write( "pdb\tresnum\tmut\tAA\tddGexp\tddGpred\n" )

# go through targets and get all data
for t in targets:
	
	# remove low and last pdbs as they're not needed
	low = f'{working_dir}/output/{t}/{t}_*low.pdb'
	last = f'{working_dir}/output/{t}/{t}_*last.pdb'
	os.system( "rm " + low )
	os.system( "rm " + last )
	
	ddg_dict = {}
	ddg_pred = []
	
	#exp_ddg = pd.read_table("%s_exp_ddg.txt" %pdb, sep=",").set_index(['res_num', 'mut'])
	file_ddg = f'{rosetta_dir}/tests/scientific/data/mp_ddg/{t}/{t}_exp_ddg.txt'
	res_num = subprocess.getoutput( "grep -v mut " + file_ddg + " | awk '{print $1}'" ).splitlines()
	mutant = subprocess.getoutput( "grep -v mut " + file_ddg + " | awk '{print $2}'" ).splitlines()
	residue = subprocess.getoutput( "grep -v mut " + file_ddg + " | awk '{print $3}'" ).splitlines()
	ddg_exp = subprocess.getoutput( "grep -v mut " + file_ddg + " | awk '{print $5}'" ).splitlines()

	# map to numbers
	res_num = list( map( int, res_num ) )
	ddg_exp = list( map( float, ddg_exp ) )

	# go through mutants
	for i in range(0, len(res_num)):
		
		#check that score files exist
		file_mut = f'{working_dir}/output/{t}/{t}_mut_%s_%s.sc' % (res_num[i], mutant[i])
		file_wt = f'{working_dir}/output/{t}/{t}_wt_%s_%s.sc' % (res_num[i], mutant[i])
		
		if not os.path.exists( file_mut ):
			f.write( file_mut + " missing!\n" )
			ddg_pred.append( 0 )
			continue

		if not os.path.exists( file_wt ):
			f.write( file_wt + " missing!\n" )
			ddg_pred.append( 0 )
			continue

		# read in scores, sorted with lowest at beginning
		score_mut = subprocess.getoutput( "grep -v SEQ " + file_mut + " | grep -v total_score | awk '{print $2}' | sort -nk1" ).splitlines()
		score_wt = subprocess.getoutput( "grep -v SEQ " + file_wt + " | grep -v total_score | awk '{print $2}' | sort -nk1" ).splitlines()
		score_mut = list( map( float, score_mut ) )
		score_wt = list( map( float, score_wt ) )

		#take average of top 3 scores from mutant and WT
		score_mut_avg = round( np.mean( score_mut[0:3] ), 3 )
		score_wt_avg = round( np.mean( score_wt[0:3] ), 3 )
		
		#subtract WT from mutant to get ddG
		ddg = round( score_mut_avg - score_wt_avg, 3 )
		ddg_pred.append( ddg )
		
		#write data points as table into ddG out file
		o.write( t + "\t" + str(res_num[i]) + "\t" + mutant[i] + "\t" + residue[i] + "\t" + str(ddg_exp[i]) + "\t" + str(ddg) + "\n")

	#append ddG info to dataframe 
	ddg_dict = {'res_num': res_num, 'mutant': mutant, 'residue': residue, 'ddg_exp': ddg_exp, 'ddg_pred': ddg_pred}

	# write data frame into dictionary	
	ddg_data[t] = ddg_dict

o.close()

# read cutoffs
protein = subprocess.getoutput( "grep -v '#' " + cutoffs + " | awk '{print $1}'" ).splitlines()
cutoffs_corr = subprocess.getoutput( "grep -v '#' " + cutoffs + " | awk '{print $2}'" ).splitlines()
cutoffs_corr = map( float, cutoffs_corr )
cutoffs_corr_dict.update( dict( zip ( protein, cutoffs_corr )))

# go through targets and compare correlation coefficient with cutoff
for t in targets:

	target_results = {}

	# run scipy linear regression
	linreg = linregress(ddg_data[t]['ddg_exp'], ddg_data[t]['ddg_pred'])
	corr = round( linreg.rvalue, 3 )
	slope = round( linreg.slope, 3 )
	intercept = round( linreg.intercept, 3 )

	# compare correlation coefficient with cutoff
	if corr < cutoffs_corr_dict[t]:
		failures.append( t )
		
	# write results and ddg data into json file
	target_results.update( {"corr" : corr, "slope" : linreg.slope, "intercept" : linreg.intercept} )
	target_results.update( ddg_data[t] )

	# write correlation coefficients into results file
	f.write( t + "\tcorr:\t" + str(corr) + "\tslope:\t" + str(slope) + "\tintercept:\t" + str(intercept) + "\n" )
		
	# update results in json
	results.update( {t : target_results} )
	
f.close()

benchmark.save_variables('debug targets mutants nstruct working_dir rosetta_dir testname results ddg_data scorefiles cutoffs_corr_dict failures')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
