#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  mhc_epitope_energy/3.plot.py
## @brief this script is part of mhc_epitope_energy scientific test
## @author Sergey Lyskov
## @author Brahm Yachnin (brahm.yachnin@rutgers.edu)

import os, sys, subprocess, math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import benchmark

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

###Distribution plot section
plotlist = [
	("plot_results_mhc_epitope.png", "mhc_epitope", cutoffs_mhc_epitope_dict),
	("plot_results_delta_mhc_epitope.png", "delta_mhc_epitope", cutoffs_delta_mhc_epitope_dict)
	]
	
for plot in plotlist:
	# inputs are header labels from the scorefile to plot, for instance "total_score" and "rmsd"
	# take these from plot
	outfile = plot[0]
	var_label = plot[1]

	# get column numbers from labels, 1-indexed
	var_index = str( subprocess.getoutput( "grep " + var_label + " " + scorefiles[0] ).split().index( var_label ) + 1 )
	
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
		var = subprocess.getoutput( "grep -v SEQUENCE " + scorefiles[i] + " | grep -v " + var_label + " | awk '{print $" + var_index + "}'" ).splitlines()
		
		# map all values to floats
		var = list( map( float, var ) )
		
		# create subplot
		plt.subplot( nrows, ncols, i+1 )
		
		# labels
		plt.xlabel( var_label )
		#plt.ylabel( y_label )
		
		# set title
		plt.title( targets[i] )

		# scatterplot of the data
		plt.hist(var)
		
		# add vertical line for cutoff
		plt.axvline(x=float( plot[2].get(targets[i]) ), color='tab:orange', linestyle='-')
		
	#save figure
	plt.tight_layout()
	plt.savefig( outfile )
	plt.close()

###2D plot section
plotlist = [ #(outfile, x_label, y_label, result_key)
	("plot_results_base_total_score_vs_mhc_epitope.png", "pct_drop_mhc", "base_total_energy", 'Num decoys with pct_drop_mhc and base_total_score better than cutoffs'),
	("plot_results_sequence_recovery_vs_mhc_epitope.png", "pct_drop_mhc", "seqrec_seqrec", 'Num decoys with pct_drop_mhc and sequence_recovery better than cutoffs'),
	("plot_results_core_sequence_recovery_vs_mhc_epitope.png", "pct_drop_mhc", "seqrec_core_seqrec", 'Num decoys with pct_drop_mhc and core_sequence_recovery better than cutoffs'),
	("plot_results_delta_packstat_vs_mhc_epitope.png", "pct_drop_mhc", "delta_packstat", 'Num decoys with pct_drop_mhc and delta_packstat better than cutoffs'),
	("plot_results_delta_buried_unsat_vs_mhc_epitope.png", "pct_drop_mhc", "delta_buried_unsat", 'Num decoys with pct_drop_mhc and delta_buried_unsat better than cutoffs'),
	("plot_results_delta_netcharge_vs_mhc_epitope.png", "pct_drop_mhc", "delta_netcharge", 'Num decoys with pct_drop_mhc and abs(delta_netcharge) better than cutoffs')
	]
	
for plot in plotlist:
	# inputs are header labels from the scorefile to plot, for instance "total_score" and "rmsd"
	# take these from plot
	outfile = plot[0]
	x_label = plot[1]
	y_label = plot[2]

	# get column numbers from labels, 1-indexed
	if x_label == "pct_drop_mhc":
		x_index = x_label
		mhc_index = str( subprocess.getoutput( "grep mhc_epitope " + scorefiles[0] ).split().index( "mhc_epitope" ) + 1 )
		delta_mhc_index = str( subprocess.getoutput( "grep delta_mhc_epitope " + scorefiles[0] ).split().index( "delta_mhc_epitope" ) + 1 )
	else:
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
		if x_label == "pct_drop_mhc":
			mhc = subprocess.getoutput( "grep -v SEQUENCE " + scorefiles[i] + " | grep -v " + y_label + " | awk '{print $" + mhc_index + "}'" ).splitlines()
			mhc = list( map( float, mhc ) )
			delta_mhc = subprocess.getoutput( "grep -v SEQUENCE " + scorefiles[i] + " | grep -v " + y_label + " | awk '{print $" + delta_mhc_index + "}'" ).splitlines()
			delta_mhc = list( map( float, delta_mhc ) )
			
			x = []
			for idx in range( 0, len(mhc) ):
				x.append( 100 * (mhc[idx] - (mhc[idx] - delta_mhc[idx]) ) / (mhc[idx] - delta_mhc[idx]) )

		else:
			x = subprocess.getoutput( "grep -v SEQUENCE " + scorefiles[i] + " | grep -v " + y_label + " | awk '{print $" + x_index + "}'" ).splitlines()
			x = list( map( float, x ) )
			
		y = subprocess.getoutput( "grep -v SEQUENCE " + scorefiles[i] + " | grep -v " + y_label + " | awk '{print $" + y_index + "}'" ).splitlines()
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
		
		# add horizontal and vertical lines for cutoff
		plt.axvline(x=float(results.get(targets[i]).get(plot[3])[1]), color='r', linestyle='-')
		plt.axhline(y=float(results.get(targets[i]).get(plot[3])[2]), color='g', linestyle='-')
		if plot[0] == "plot_results_delta_netcharge_vs_mhc_epitope.png":	#Add a line at -1*cutoff too for delta netcharge, since we look at the absolute value
			plt.axhline(y=-float(results.get(targets[i]).get(plot[3])[2]), color='g', linestyle='-')
		
	#save figure
	plt.tight_layout()
	plt.savefig( outfile )
	plt.close()

benchmark.save_variables('targets nstruct working_dir testname results outfile cutoffs_mhc_epitope_dict cutoffs_delta_mhc_epitope_dict failures failures_dict')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
