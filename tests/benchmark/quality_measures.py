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

import os, os.path, sys, imp, shutil, json
import numpy as np

#=======================================
def check_all_values_below_cutoff( col, cutoff, tag, filehandle ):

	out = "All " + tag + "s < cutoff"
	filehandle.write( out + " " + str(cutoff) + "\t" )

	if all( i <= cutoff for i in col ):
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

