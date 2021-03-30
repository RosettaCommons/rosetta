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

#number of subplots
ncols = 3
nrows = 1
if len( targets ) < 4:
	ncols = len( targets )
else:
	nrows = math.ceil( len( targets ) / 4 )

# figure size
width = 6 * ncols
height = 6 * nrows

plt.rc("font", size=20)
plt.rcParams['figure.figsize'] = width, height #width, height

# go through scorefiles
for i in range(0, len(targets)):
	
	t = targets[i]
	
	# get x and y
	x = list( map( float, ddg_data[t]['ddg_exp'] ) )
	y = list( map( float, ddg_data[t]['ddg_pred'] ) )
	
	# create subplot
	plt.subplot( nrows, ncols, i+1 )
	
	# x and y labels
	plt.xlabel( x_label )
	plt.ylabel( y_label )
	
	# set title
	plt.title( t + ': Pearson = %0.3f' % results[t]['corr'] )

	# scatterplot of the data
	plt.plot( x, y, 'ko' )
	
	# add linear regression line: slope + intercept
	xlist = list(np.linspace( min(x), max(x), 1000 ))
	ylist = [(results[t]['slope'] * j + results[t]['intercept']) for j in xlist ]
	plt.plot( xlist, ylist, c='b' )
	
	# add correlation coefficient as text
#	plt.text( min(x), max(y)-0.5, 'pearson = %0.3f' % results[t]['corr'] )
	
	# add horizontal and vertical lines for axis at 0
	plt.axvline(x=0.0, color='lightgray', linestyle='-')
	plt.axhline(y=0.0, color='lightgray', linestyle='-')
	
	# x axis limits
#	if targets[i] in failures:
#		plt.xlim( left=0 )
#	else:
#		plt.xlim( 0, 10 )
	
#save figure
plt.tight_layout()
plt.savefig( outfile )


benchmark.save_variables('debug targets mutants nstruct working_dir testname results scorefiles cutoffs_corr_dict ddg_data failures')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
