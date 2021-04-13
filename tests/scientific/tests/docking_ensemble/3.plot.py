#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  docking_ensemble/3.plot.py
## @brief this script is part of docking ensemble scientific test
## @author Sergey Lyskov
## @author Ameya Harmalkar

import os, sys, subprocess, math
import matplotlib
import numpy as np
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import benchmark

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

# inputs are header labels from the scorefile to plot, for instance "total_score" and "rmsd"
# => it figures out the column numbers from there
x_label = "Irms"
y_label = "I_sc"
z_label = "CAPRI_rank"
outfile = "plot_results.png"

# get column numbers from labels, 1-indexed
x_index = str( subprocess.getoutput( "grep " + x_label + " " + scorefiles[0] ).split().index( x_label ) + 1 )
y_index = str( subprocess.getoutput( "grep " + y_label + " " + scorefiles[0] ).split().index( y_label ) + 1 )
z_index = str( subprocess.getoutput( "grep " + z_label + " " + scorefiles[0] ).split().index( z_label ) + 1 )

#number of subplots
ncols = 4
nrows = 3

# figure size
width = 7.5 * ncols
height = 6 * nrows

plt.rc("font", size=20)
plt.rcParams['figure.figsize'] = width, height #width, height

# get number of fields in scorefile
nfields = len( subprocess.getoutput( "grep -v SEQUENCE " + scorefiles[0] + " | grep " + y_label + " | head -n1" ).split())

# go through scorefiles
for i in range( 0, len( scorefiles ) ):
	
	# read in score file
	x = subprocess.getoutput( "awk '{if(NF==" + str(nfields) + ") print}' " + scorefiles[i] + " | grep -v SEQUENCE | grep -v " + y_label + " |  sort -nk6 | awk '{print $" + x_index + "}'" ).splitlines()
	y = subprocess.getoutput( "awk '{if(NF==" + str(nfields) + ") print}' " + scorefiles[i] + " | grep -v SEQUENCE | grep -v " + y_label + " |  sort -nk6 | awk '{print $" + y_index + "}'" ).splitlines()
	z = subprocess.getoutput( "awk '{if(NF==" + str(nfields) + ") print}' " + scorefiles[i] + " | grep -v SEQUENCE | grep -v " + y_label + " |  sort -nk6 | awk '{print $" + z_index + "}'" ).splitlines()
    
	# map all values to floats
	x = list( map( float, x ) )
	y = list( map( float, y ) )
	z = list( map( float, z ) )

	# containers to store the type of target
	incorrect = []
	acceptable = []
	medium = []
	high = []
	
	for j in range(len(x)):
		score_terms = [x[j], y[j]]
		if z[j] == 0.0:
			incorrect.append( score_terms )
		elif z[j] == 1.0:
			acceptable.append( score_terms )
		elif z[j] == 2.0:
			medium.append( score_terms )
		elif z[j] == 3.0:
			high.append( score_terms )

	incorrect = list(map( list, zip(*incorrect) ))
	if len(acceptable) > 0:
		acceptable = list( map( list, zip(*acceptable) ))
	if len(medium) > 0:
		medium = list(map( list, zip(*medium) ))
	if len(high) > 0:
		high = list(map( list, zip(*high) ))
		
	# common settings for the subplots
	s = 25

	plt.subplot( nrows, ncols, i+1 )
	
	if( len(incorrect) > 0 ):
		plt.scatter( incorrect[0], incorrect[1], c='black', zorder=3, s=s )
	if( len(acceptable) > 0 ):
		plt.scatter( acceptable[0], acceptable[1], c='orange', zorder=3, s=s )
	if( len(medium) > 0 ):
		plt.scatter( medium[0], medium[1], c='red', zorder=3, s=s )
	if( len(high) > 0 ):
		plt.scatter( high[0], high[1], c='green', zorder=3, s=s )

	plt.xlim([0,20])
	plt.ylim([min(y)- 0.5, 0.5])
	plt.title( targets[i] )
	plt.xlabel("Irms")
	plt.ylabel("Interface Score (REU)")

	# add horizontal and vertical lines for cutoff
	plt.axvline(x=float(cutoffs_rmsd_dict[targets[i]]), color='b', linestyle='-')
#	plt.axhline(y=float(cutoffs_score_dict[targets[i]]), color='b', linestyle='-')
	
#save figure
plt.tight_layout()
plt.savefig( outfile )

benchmark.save_variables('debug targets nstruct working_dir testname results outfile cutoffs_rmsd_dict cutoffs_score_dict failures')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
