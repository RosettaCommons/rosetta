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
from copy import deepcopy
import numpy as np
import benchmark
from scipy.stats import linregress, spearmanr
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
cutoffs_rho_dict = {}
failures = []

# read in all data into ddg_data, a list of dicts
ddg_data = {}

# open results output file
f = open( outfile, "w" )

# write data into ddg out file
o = open( ddg_out, "w" )
o.write( "pdb\tresnum\tmut\tAA\tddGexp\tddGpred\n" )

# experimental reference list for each target
exp_ref_list = {
	"1py6" : [['1py6_all', '1py6_exp_ddg.txt'], 
				['1py6_ref1', '1py6_exp_ddg_ref1.txt'], ['1py6_ref2', '1py6_exp_ddg_ref2.txt'], 
				['1py6_ref3', '1py6_exp_ddg_ref3.txt'], ['1py6_ref4', '1py6_exp_ddg_ref4.txt'], 
				['1py6_ref5', '1py6_exp_ddg_ref5.txt'], ['1py6_ref6', '1py6_exp_ddg_ref6.txt'], 
				['1py6_ref7', '1py6_exp_ddg_ref7.txt'], ['1py6_ref8', '1py6_exp_ddg_ref8.txt'],], 
	"2k74" : [['2k74_all', '2k74_exp_ddg.txt'], 
				['2k74_ref1', '2k74_exp_ddg_ref1.txt'], ['2k74_ref2', '2k74_exp_ddg_ref2.txt'],],
	"2xov" : [['2xov_all', '2xov_exp_ddg.txt'], ['2xov_ref1', '2xov_exp_ddg_ref1.txt'], 
				['2xov_ref2', '2xov_exp_ddg_ref2.txt'], ['2xov_ref3', '2xov_exp_ddg_ref3.txt']]
}

# go through targets and get all data
for t in targets:

	# remove low and last pdbs as they're not needed
	low = f'{working_dir}/output/{t}/{t}_*low.pdb'
	last = f'{working_dir}/output/{t}/{t}_*last.pdb'
	os.system( "rm " + low )
	os.system( "rm " + last )

	for ref in exp_ref_list[t]:

		ddg_dict = {}
		ddg_pred = []

		# exp_ddg = pd.read_table("%s_exp_ddg.txt" %pdb, sep=",").set_index(['res_num', 'mut'])
		file_ddg = f'{rosetta_dir}/tests/scientific/data/mp_ddg/{t}/{ref[1]}'
		res_num = subprocess.getoutput( "grep -v mut " + file_ddg + " | awk '{print $1}'" ).splitlines()
		mutant = subprocess.getoutput( "grep -v mut " + file_ddg + " | awk '{print $2}'" ).splitlines()
		residue = subprocess.getoutput( "grep -v mut " + file_ddg + " | awk '{print $3}'" ).splitlines()
		ddg_exp = subprocess.getoutput( "grep -v mut " + file_ddg + " | awk '{print $5}'" ).splitlines()

		# map to numbers
		res_num = list( map( int, res_num ) )
		ddg_exp = list( map( float, ddg_exp ) )

		# go through mutants
		for i in range(0, len(res_num)):

			# check that score files exist
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

			# take average of top 3 scores from mutant and WT
			score_mut_avg = round( np.mean( score_mut[0:3] ), 3 )
			score_wt_avg = round( np.mean( score_wt[0:3] ), 3 )

			# subtract WT from mutant to get ddG
			ddg = round( score_mut_avg - score_wt_avg, 3 )
			ddg_pred.append( ddg )

			# write data points as table into ddG out file
			o.write( ref[0] + "\t" + str(res_num[i]) + "\t" + mutant[i] + "\t" + residue[i] + "\t" + str(ddg_exp[i]) + "\t" + str(ddg) + "\n")

		# append ddG info to dataframe
		ddg_dict = {'res_num': res_num, 'mutant': mutant, 'residue': residue, 'ddg_exp': ddg_exp, 'ddg_pred': ddg_pred}

		# write data frame into dictionary	
		ddg_data[ref[0]] = ddg_dict

		# store values for combined plot/analysis
		if not ref[0].endswith('_all'):
			if 'all' in ddg_data.keys():
				for subkey in ddg_data['all'].keys():
					ddg_data['all'][subkey] = ddg_data['all'][subkey] + deepcopy(ddg_dict[subkey])
			else:
				ddg_data['all'] = deepcopy(ddg_dict)

o.close()

# read cutoffs
protein = subprocess.getoutput( "grep -v '#' " + cutoffs + " | awk '{print $1}'" ).splitlines()
cutoffs_corr = subprocess.getoutput( "grep -v '#' " + cutoffs + " | awk '{print $2}'" ).splitlines()
cutoffs_rho = subprocess.getoutput( "grep -v '#' " + cutoffs + " | awk '{print $3}'" ).splitlines()
cutoffs_corr = map( float, cutoffs_corr )
cutoffs_rho = map( float, cutoffs_rho )
cutoffs_corr_dict.update( dict( zip ( protein, cutoffs_corr )))
cutoffs_rho_dict.update( dict( zip ( protein, cutoffs_rho )))

# go through targets and compare correlation coefficient with cutoff
def correl_cal(file_pointer, result_dict, sub_ref, ddg_data):
	target_results = {}

	# run scipy linear regression
	linreg = linregress(ddg_data[sub_ref]['ddg_exp'], ddg_data[sub_ref]['ddg_pred'])
	corr = round( linreg.rvalue, 3 )
	slope = round( linreg.slope, 3 )
	intercept = round( linreg.intercept, 3 )

	# run scipy spearmanr
	rho, rho_pval = spearmanr(ddg_data[sub_ref]['ddg_exp'], ddg_data[sub_ref]['ddg_pred'])
	spearman_rho = round( rho, 3 )
	spearman_rho_pval = round( rho_pval, 3 )

	# compare correlation coefficient with cutoff
	if sub_ref in cutoffs_corr_dict.keys():
		if corr < cutoffs_corr_dict[sub_ref]:
			failures.append( f"{sub_ref} (pearson)" )

	# compare rho with cutoff
	if sub_ref in cutoffs_rho_dict.keys():
		if spearman_rho < cutoffs_rho_dict[sub_ref]:
			failures.append( f"{sub_ref} (spearman)" )

	# write results and ddg data into json file
	target_results.update( {"corr" : corr, "slope" : linreg.slope, "intercept" : linreg.intercept, "rho" : rho, "pval" : rho_pval} )
	target_results.update( ddg_data[sub_ref] )

	# write correlation coefficients into results file
	file_pointer.write( sub_ref + "\tcorr:\t" + str(corr) + "\tslope:\t" + str(slope) + "\tintercept:\t" + str(intercept) + 
		"\trho:\t" + str(spearman_rho) + "\tpval:\t" + str(spearman_rho_pval) + "\n" )

	# update results in json
	result_dict.update( { sub_ref: target_results} )

for t in targets:

	for ref in exp_ref_list[t]:

		correl_cal(f, results, ref[0], ddg_data)

correl_cal(f, results, 'all', ddg_data)

f.close()

benchmark.save_variables('debug targets mutants nstruct working_dir rosetta_dir testname results ddg_data scorefiles cutoffs_corr_dict exp_ref_list failures')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
