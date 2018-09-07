#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  cartesian_relax/3.plot.py
## @brief this script is part of cartesian_relax scientific test
## @author Sergey Lyskov

import os, sys, subprocess, math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import benchmark

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

# inputs are header labels from the scorefile to plot, for instance "total_score" and "rmsd"
# => it figures out the column numbers from there
x_label = "rmsd"
y_label = "total_score"
outfile = "score_vs_rmsd.png"

# get column numbers from labels, 1-indexed
x_index = str( subprocess.getoutput( "grep " + x_label + " " + scorefiles[0] ).split().index( x_label ) + 1 )
y_index = str( subprocess.getoutput( "grep " + y_label + " " + scorefiles[0] ).split().index( y_label ) + 1 )

#number of subplots
ncols = 4
nrows = 1
if len( targets ) < 4:
    ncols = len( targets )
else:
    nrows = math.ceil( len( targets ) / 4 )

# figure size
width = 7.5 * ncols
height = 6 * nrows

plt.rc("font", size=20)
plt.rcParams['figure.figsize'] = width, height #width, height

# go through scorefiles
for i in range( 0, len( scorefiles ) ):

	# read in score file
	x = subprocess.getoutput( "grep -v SEQUENCE " + scorefiles[i] + " | grep -v " + y_label + " | awk '{print $" + x_index + "}'" ).splitlines()
	y = subprocess.getoutput( "grep -v SEQUENCE " + scorefiles[i] + " | grep -v " + y_label + " | awk '{print $" + y_index + "}'" ).splitlines()
	
    # map all values to floats
	x = list( map( float, x ) )
	y = list( map( float, y ) )
	
	# create subplot
	plt.subplot( nrows, ncols, i+1 )
	
	# x and y labels
	plt.xlabel( x_label )
	plt.ylabel( y_label )
	
	# set title
	plt.title( targets[i] )

	# scatterplot of the data
	plt.plot(x, y, 'ko')
	
	# x axis limits
	plt.xlim( 0, 10 )
	
#save figure
plt.tight_layout()
plt.savefig( outfile )

benchmark.save_variables('targets nstruct working_dir testname results outfile')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
	
	