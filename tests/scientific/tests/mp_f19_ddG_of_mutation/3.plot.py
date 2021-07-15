#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  mp_f19_ddG_of_mutation/3.analyze.py
## @brief this script is part of mp_f19_ddG_of_mutation scientific test
## @author Sergey Lyskov, Rebecca F. Alford (ralford3@jhu.edu)

import os, sys, subprocess, math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import benchmark

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

# Output file for plotting results
outfile = "plot_results.png"

# Read information from ddG data files
exp_label = "experimental_ddG"
pred_label = "predicted_ddG"
aa_label = "Mut"

# get column numbers from labels, 1-indexed
exp_index = str( subprocess.getoutput( "grep " + exp_label + " " + ddG_files[0] ).split().index( exp_label ) + 1 )
pred_index = str( subprocess.getoutput( "grep " + pred_label + " " + ddG_files[0] ).split().index( pred_label ) + 1 )
aa_index = str( subprocess.getoutput( "grep " + aa_label + " " + ddG_files[0] ).split().index( aa_label ) + 1 )

#number of subplots
ncols = 3
nrows = 1

# figure size
width = 7.5 * ncols
height = 6 * nrows

plt.rc("font", size=20)
plt.rcParams['figure.figsize'] = width, height #width, height

# Correlation coefficients from Alford et al. 2019
#updating R values for the plots. In the paper the results were for R2. 
moon_fleming = 0.922
mcdonald_fleming = 0.158
marx_fleming = 0.919

lit_values = [ moon_fleming, mcdonald_fleming, marx_fleming ]

# go through scorefiles
for i in range( 0, len( ddG_files ) ):

	# read in score file, scores are sorted, first one is lowest
	exp = subprocess.getoutput( "grep -v SEQUENCE " + ddG_files[i] + " | grep -v " + exp_label + " | awk '{print $" + exp_index + "}'" ).splitlines()
	pred = subprocess.getoutput( "grep -v SEQUENCE " + ddG_files[i] + " | grep -v " + pred_label + " | awk '{print $" + pred_index + "}'" ).splitlines()
	aa = subprocess.getoutput( "grep -v SEQUENCE " + ddG_files[i] + " | grep -v " + aa_label + " | awk '{print $" + aa_index + "}'" ).splitlines()

	# map values to floats (were strings)
	exp = list( map( float, exp ))
	pred = list( map( float, pred ))
	aa = list( map( str, aa ))

	# Make numpy arrays
	exp_arr = np.asarray( exp )
	pred_arr = np.asarray( pred )
	aa_arr = np.asarray( aa )

	# Remove prolines
	exp_no_proline = exp_arr[ np.where( aa_arr != "P" ) ]
	pred_no_proline = pred_arr[ np.where( aa_arr != "P" ) ]
	aa_no_proline = aa_arr[ np.where( aa_arr != "P" ) ]

	# create subplot
	plt.subplot( nrows, ncols, i+1 )
	
	# x and y labels
	plt.xlabel( exp_label )
	plt.ylabel( pred_label )
	
	# set title
	plt.title( targets[i] )

	# Plot the cutoff/reference best line of fit (red, should not show)
	x = np.linspace(-10, 10, num=200)
	abline_values2 = [ cutoffs_slope_dict[targets[i]] * j + cutoffs_intercept_dict[targets[i]] for j in x ]
	plt.plot(x, abline_values2, 'lightgray', linewidth = 20.0 )

	# Plot the line of best fit
	abline_values = [ results[targets[i]][1] * j + results[targets[i]][2] for j in x]
	plt.plot(x, abline_values, 'b' )
	
	# scatterplot of the data
	plt.plot(exp_no_proline, pred_no_proline, 'ko')

	# Plot text labels
	corr_coeff = np.corrcoef( exp_no_proline, pred_no_proline )[1,0]
	plt.text( -5, 4, "corr=" + str(round(corr_coeff,3)) + "\n lit_coeff=" + str(lit_values[i]) , horizontalalignment='center', verticalalignment='center' )


#save figure
plt.tight_layout()
plt.savefig( outfile )

benchmark.save_variables('debug targets nstruct working_dir testname results ddG_files cutoffs_corr_dict cutoffs_slope_dict cutoffs_intercept_dict failures')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
