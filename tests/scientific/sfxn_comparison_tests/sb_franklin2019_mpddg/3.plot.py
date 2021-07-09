#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  mpddg/3.plot.py
## @brief this script is part of mpddg scientific test
## @author Julia Koehler Leman, Hope Woods, Johanna Tiemann

import os, sys, subprocess, math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import benchmark

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

# file names
x_label = "Experimental ddG (kcal/mol)"
y_label = "Predicted ddG (REU)"
outfile = "plot_results.png"

# number of subplots
ncols = 4
nrows = 1
n_plots = ['all'] + [ref[0] for t in targets for ref in exp_ref_list[t]]
if len( n_plots ) < 5:
	ncols = len( n_plots )
else:
	nrows = math.ceil( len( n_plots ) / 4 )

# figure size
width = 6 * ncols
height = 6 * nrows

plt.rc("font", size=20)
plt.rcParams['figure.figsize'] = width, height #width, height

# go through scorefiles
for i in range(0, len(n_plots)):

	ref = n_plots[i]
		
	# get x and y
	x = list( map( float, ddg_data[ref]['ddg_exp'] ) )
	y = list( map( float, ddg_data[ref]['ddg_pred'] ) )

	# create subplot
	if i ==0:
		ax1 = plt.subplot( nrows, ncols, i+1 )
	else:
		ax1 = plt.subplot( nrows, ncols, i+1, aspect=1 )

	# x and y labels
	plt.xlabel( x_label )
	plt.ylabel( y_label )

	# set title
	plt.title( ref + ': Pearson = %0.3f' % results[ref]['corr'] + '\n' + 
		ref + ': Spearman = %0.3f' % results[ref]['rho'] )

	# scatterplot of the data
	plt.plot( x, y, 'ko' )

	# add linear regression line: slope + intercept
	xlist = list(np.linspace( min(x), max(x), 1000 ))
	ylist = [(results[ref]['slope'] * j + results[ref]['intercept']) for j in xlist ]
	plt.plot( xlist, ylist, c='b' )

	# add correlation coefficient as text
#	plt.text( min(x), max(y)-0.5, 'pearson = %0.3f' % results[t]['corr'] )

	# add horizontal and vertical lines for axis at 0
	plt.axvline(x=0.0, color='lightgray', linestyle='-')
	plt.axhline(y=0.0, color='lightgray', linestyle='-')

	# get same x,y limits as for overall plot
	if ref != 'all':
		ref_ref = f"{ref.split('_')[0]}_all"
		plt.xlim(min(ddg_data[ref_ref]['ddg_exp'])-0.5, max(ddg_data[ref_ref]['ddg_exp'])+0.5)
		plt.ylim(min(ddg_data[ref_ref]['ddg_pred'])-0.5, max(ddg_data[ref_ref]['ddg_pred'])+0.5)

	# make plots square
	x0,x1 = ax1.get_xlim()
	y0,y1 = ax1.get_ylim()
	ax1.set_aspect((x1-x0)/(y1-y0))

	# x axis limits
#	if targets[i] in failures:
#		plt.xlim( left=0 )
#	else:
#		plt.xlim( 0, 10 )

# save figure
plt.tight_layout()
plt.savefig( outfile )

benchmark.save_variables('debug targets mutants nstruct working_dir testname results scorefiles cutoffs_corr_dict ddg_data exp_ref_list failures')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
