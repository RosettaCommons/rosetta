#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   benchmark.py
## @brief  Helper functions for quality measures
## @author Sergey Lyskov

import os, os.path, sys, imp, shutil, json, math
import numpy as np

#=======================================
def check_all_values_below_cutoff( col, cutoff, tag, filehandle ):
	
	return check_xpercent_values_below_cutoff( col, cutoff, tag, filehandle, 100 )

#=======================================
def check_xpercent_values_below_cutoff( col, cutoff, tag, filehandle, percentage ):

	out = "All " + tag + "s < cutoff"
	filehandle.write( out + " " + str(cutoff) + "\t" )
	
	# sort the values from smallest to largest, then take the first x records
	col = sorted(col)
	partial_col = col[:int( math.ceil(percentage * len( col )/100.0) )]

	if all( i <= cutoff for i in partial_col ):
		value = True
	else:
		value = False

	filehandle.write( str( value ) + "\n" )
	return {out : value}

#=======================================
def check_rmsd_of_topscoring( rmsd_col_sorted, cutoff, filehandle ):

	out = "rmsd of topscoring model < " + str( cutoff )
	filehandle.write( out + "\t" )

	if rmsd_col_sorted[0] <= cutoff:
		value = True
	else:
		value = False

	filehandle.write( str( value ) + "\n" )
	return {out : value}

#=======================================
def check_range( col, tag, filehandle ):

	filehandle.write( tag + "\tmin, max, avg, std:" + '% 12.3f % 12.3f % 12.3f % 12.3f\n' % (min( col ), max( col ), np.mean( col ), np.std( col )) )
	value = { "min" : round(min( col ), 4), "max" : round(max( col ), 4), "avg" : round(np.mean( col ), 4), "std" : round(np.std( col ), 4) }
	return {tag : value}

#=======================================
def check_for_2d_top_values( xcol, ycol, xtag, ytag, xcutoff, ycutoff, filehandle, xminimize = True, yminimize = True, xabsolute = True, yabsolute = True ):
	"""
	Check if there are decoys that are in both the top xcutoff% and top ycutoff%.
	Returns a tuple, containing (number of passes, x_absolute_cutoff, y_absolute_cutoff).
	xcol, ycol = Arrays of scoreterm values to be checked
	xtag, ytag = Used to genereate the out tag.  Should usually be the scoreterm names being checked (but is arbitrary)
	xcutoff, ycutoff = The absolute or percentile cutoff for xcol and ycol, respectively.  If percentile, should be between 0 and 100.
	filehandle = The result file to write to.
	xminimize, yminimize = By default (True), the "top" values are those less than the cutoff.  If false, the "top" values are those greater than the cutoff.
	xabsolute, yabsolute = Should the threshold be using an absolute cutoff (True) or a percentile cutoff (False).
	"""
	
	assert(len(xcol) == len(ycol))
	
	out = "Num decoys with " + xtag + " and " + ytag + " better than cutoffs"
	filehandle.write( out + " " + str(xcutoff) + " and " + str(ycutoff) + ":\t" )
	
	#Set the cutoff values in x and y
	x_val = xcutoff
	y_val = ycutoff
	#Get the percentile values for xcol and ycol if we ask for it
	if not xabsolute: x_val = np.percentile( xcol, xcutoff )
	if not yabsolute: y_val = np.percentile( ycol, ycutoff )
	
	#Check each value in xcol to see if it is greater/less than xperc_val
	#Then check if the corresponding value in ycol is greater/less than yperc_val
	#If both passing conditions are met, append True to passed.  If either condition fails, append False to passed.
	passed = []
	for idx in range(len(xcol)):
		if xminimize:
			if xcol[idx] <= x_val:
				if yminimize:
					if ycol[idx] <= y_val:
						passed.append(True)
					else:
						passed.append(False)
				else: #ymaximize
					if ycol[idx] >= y_val:
						passed.append(True)
					else:
						passed.append(False)
			else:
				passed.append(False)
		else: #xmaximize
			if xcol[idx] >= x_val:
				if yminimize:
					if ycol[idx] <= y_val:
						passed.append(True)
					else:
						passed.append(False)
				else: #ymaximize
					if ycol[idx] >= y_val:
						passed.append(True)
					else:
						passed.append(False)
			else:
				passed.append(False)
	assert(len(passed) == len(xcol))
	assert(len(passed) == len(ycol))
				
	#The number of passes is the sum of passed (true = 1, false = 0)
	npasses = sum(passed)
	filehandle.write( str( npasses ) + " (" + str(x_val) + ", " + str(y_val) + ")\n" )
	return {out : (npasses,x_val,y_val)}