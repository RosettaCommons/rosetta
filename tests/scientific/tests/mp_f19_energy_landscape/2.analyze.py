#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  mp_f19_energy_landscape/3.plot.py
## @brief Plot energy landscape data for the mp_19_energy_landscape scientific benchmark test
## @author Sergey Lyskov, Rebecca F. Alford (ralford3@jhu.edu)

import os, sys, subprocess, math
import numpy as np
import benchmark
import energy_landscape_metrics
from benchmark.util import quality_measures as qm

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into 	s
config = benchmark.config()

results = {}
energy_landscapes = []
cutoffs_ddG_ins_dict = {}
cutoffs_zmin_dict = {}
cutoffs_anglemin_dict = {}
failures = []

# Define output and cutoffs filename
outfile = "result.txt"
cutoffs = "cutoffs"

# Get a list of files containing energy landscapes for each target
energy_landscapes.extend( [ f'{working_dir}/output/{t}/{t}_franklin2019_landscape.dat' for t in targets ] )

# read cutoffs
protein = subprocess.getoutput( "grep -v '#' " + cutoffs + " | awk '{print $1}'" ).splitlines()
cutoffs_ddG_ins = subprocess.getoutput( "grep -v '#' " + cutoffs + " | awk '{print $2}'" ).splitlines()
cutoffs_zmin = subprocess.getoutput( "grep -v '#' " + cutoffs + " | awk '{print $3}'" ).splitlines()
cutoffs_anglemin = subprocess.getoutput( "grep -v '#' " + cutoffs + " | awk '{print $4}'" ).splitlines()

cutoffs_ddG_ins = map( float, cutoffs_ddG_ins ) 
cutoffs_zmin = map( float, cutoffs_zmin )
cutoffs_anglemin = map( float, cutoffs_anglemin )

cutoffs_ddG_ins_dict.update( dict( zip ( protein, cutoffs_ddG_ins )))
cutoffs_zmin_dict.update( dict( zip ( protein, cutoffs_zmin )))
cutoffs_anglemin_dict.update( dict( zip( protein, cutoffs_anglemin )))

# open results output file
f = open( outfile, "w" )
f.write( "target\tdG\tzmin\tanglemin" )

# go through the energy landscape files for each target
for i in range( 0, len( energy_landscapes ) ): 

	target_results = {}

	# Read in the energy landscape file
	if ( not os.path.isfile( energy_landscapes[i] ) ): 
		sys.exit( "Output data file " + energy_landscapes[i] + " not found!" )

	with open( energy_landscapes[i], 'rt' ) as g: 
		contents = g.readlines()
		contents = [ x.strip() for x in contents ]
		contents = [ x.split(" ") for x in contents ]

	zcoords = []
	angles = []
	total_scores = []

	for x in range(1, len(contents)): 
		zcoords.append( float( contents[x][0] ) )
		angles.append( float( contents[x][1] ))
		total_scores.append( float( contents[x][2] ))

	zcoords_arr = np.asarray( zcoords )
	angles_arr = np.asarray( angles )
	total_scores_arr = np.asarray( total_scores )

	# Calculate the ddG of insertion and the minimum energy orientations
	dG_transfer = energy_landscape_metrics.compute_dG_transfer_energy( zcoords_arr, angles_arr, total_scores_arr )
	zmin, anglemin = energy_landscape_metrics.compute_minimum_energy_orientation( zcoords_arr, angles_arr, total_scores_arr )

	# Is the ddG_ins within +/- 1 REU of the previously established value? 
	f.write( "\n" + targets[i] + "\t" )
	f.write( str(dG_transfer) + "\t" )
	val_cutoff = cutoffs_ddG_ins_dict[targets[i]]
	target_results.update( dG_transfer = dG_transfer )

	if ( not ( dG_transfer <= (val_cutoff+1) and dG_transfer >= (val_cutoff-1) ) ): 
		f.write( "fails in dG: " + str(val_cutoff) + "\t" )
		if ( targets[i] not in failures ):
			failures.append( targets[i] )

	# Is the mimimum energy orinetation within +/- 1angstrom and +/- 5 degrees of the previously estabished value? 
	f.write( str(round(zmin,3)) + "\t" )
	val_cutoff = cutoffs_zmin_dict[targets[i]]
	target_results.update( zmin = zmin )

	if ( not ( zmin <= (val_cutoff+2) and zmin >= (val_cutoff-2) ) ): 
		f.write(  "fails in zmin: " + str(val_cutoff) + "\t" )
		if ( targets[i] not in failures ):
			failures.append( targets[i] )

	# test for angle_min
	f.write( str(anglemin) + "\t" )
	val_cutoff = cutoffs_anglemin_dict[targets[i]]
	target_results.update( anglemin = anglemin )
	
	if ( not ( anglemin <= (val_cutoff+10) and anglemin >= (val_cutoff-10) ) ): 
		f.write(  "fails in anglemin: " + str(val_cutoff) + "\t" )
		if ( targets[i] not in failures ):
			failures.append( targets[i] )

	results.update( {targets[i] : target_results} )

f.close()

benchmark.save_variables('debug targets working_dir testname results energy_landscapes cutoffs_ddG_ins_dict cutoffs_zmin_dict cutoffs_anglemin_dict failures')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
