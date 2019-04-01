#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  cartesian_relax/3.plot.py
## @brief this script is part of the mp_lipid_acc scientific test
## @author JKLeman

import os, sys, subprocess, math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import benchmark

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

# inputs are header labels from the scorefile to plot, for instance "total_score" and "rmsd"
# => it figures out the column numbers from there
x_label = "protein"
y_label = ['acc', 'sens', 'spec']
outfile = "plot_results.png"

# get column numbers from labels, 1-indexed
x_index = str( subprocess.getoutput( "grep " + x_label + " result.txt" ).split().index( x_label ) + 1 )
y_index = [str( subprocess.getoutput( "grep " + i + " result.txt" ).split().index( i ) + 1 ) for i in y_label]

#number of subplots
ncols = 1
nrows = 3

# figure size
width = 40 * ncols
height = 6 * nrows

plt.rc("font", size=24)
plt.rcParams['figure.figsize'] = width, height #width, height

# go through scorefiles
for i in range( 0, 3 ):

	# read in score file
	tag = subprocess.getoutput( "grep -v protein result.txt | awk '{print $" + x_index + "}'" ).splitlines()
	y = subprocess.getoutput( "grep -v protein result.txt | awk '{print $" + y_index[i] + "}'" ).splitlines()
	
	# map all values to floats
	x = [i for i in range(1, len(tag)+1)]
	y = list( map( float, y ) )
	
	# create subplot
	plt.subplot( nrows, ncols, i+1 )
	
	# x and y labels
	plt.xticks( x, tag, rotation='vertical', fontsize=14 )
	plt.xlabel( x_label )
	plt.ylabel( y_label[i] )

	# plot limits
	plt.xlim( 0, 223 )
	
	# set title
	plt.title( y_label[i] )

	# scatterplot of the data
	plt.bar(x, y, width=0.8, color='silver')
	
#save figure
plt.tight_layout()
plt.savefig( outfile )

benchmark.save_variables('targets nstruct working_dir testname results outfile cutoffs_rmsd_dict cutoffs_score_dict failures')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
