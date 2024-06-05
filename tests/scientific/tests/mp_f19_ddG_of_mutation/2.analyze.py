#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  mp_f19_ddG_of_mutation/2.analyze.py
## @brief this script is part of mp_f19_ddG_of_mutation scientific test
## @author Sergey Lyskov, Rebecca F. Alford (ralford3@jhu.edu)

import os, os.path, sys, subprocess, math
import numpy as np
import benchmark
from benchmark.util import quality_measures as qm

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into 	s
config = benchmark.config()

results = {}
ddG_files = []
outfile = "result.txt"
cutoffs = "cutoffs"
failures = []

# scorefiles and logfiles
ddG_files.extend( [ f'{working_dir}/output/{t}/ddG_franklin2019.dat' for t in targets ] )

exp_label = "experimental_ddG"
pred_label = "predicted_ddG"
aa_label = "Mut"

cutoffs_corr_dict = {}
cutoffs_slope_dict = {}
cutoffs_intercept_dict = {}

# read cutoffs
protein = subprocess.getoutput( "grep -v '#' " + cutoffs + " | awk '{print $1}'" ).splitlines()
cutoffs_corr = subprocess.getoutput( "grep -v '#' " + cutoffs + " | awk '{print $2}'" ).splitlines()
cutoffs_corr = map( float, cutoffs_corr )
cutoffs_corr_dict.update( dict( zip ( protein, cutoffs_corr )))

cutoffs_slope = subprocess.getoutput( "grep -v '#' " + cutoffs + " | awk '{print $3}'" ).splitlines()
cutoffs_slope = map( float, cutoffs_slope )
cutoffs_slope_dict.update( dict( zip ( protein, cutoffs_slope )))

cutoffs_intercept = subprocess.getoutput( "grep -v '#' " + cutoffs + " | awk '{print $4}'" ).splitlines()
cutoffs_intercept = map( float, cutoffs_intercept )
cutoffs_intercept_dict.update( dict( zip ( protein, cutoffs_intercept )))

# open results output file
f = open( outfile, "w" )
f.write( "target\tcorrcoeff\tslope\tintercept\n" )

# go through scorefiles of targets
for i in range( 0, len( ddG_files ) ):
    if not os.path.isfile(ddG_files[i]):
        raise RuntimeError("Cannot find datafile: "+ddG_files[i])

    try:
        # get column numbers from labels, 1-indexed
        exp_index = str( subprocess.getoutput( "grep " + exp_label + " " + ddG_files[0] ).split().index( exp_label ) + 1 )
        pred_index = str( subprocess.getoutput( "grep " + pred_label + " " + ddG_files[0] ).split().index( pred_label ) + 1 )
        aa_index = str( subprocess.getoutput( "grep " + aa_label + " " + ddG_files[0] ).split().index( aa_label ) + 1 )

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

        # Calculate a correlation coefficient
        pearson_correl = np.corrcoef( exp_no_proline, pred_no_proline ).item((0,1))
        best_fit = np.polyfit( exp_no_proline, pred_no_proline, 1 )

        slope = best_fit[0]
        intercept = best_fit[1]

        # check for correlation above cutoff
        #Its a stochastic value, exact value may not match. giving a range of -0.05
        f.write( targets[i] + "\t" )
        f.write( str(round(pearson_correl, 3)) + "\t" )
        if ( cutoffs_corr_dict[targets[i]]-0.05 > pearson_correl ):
            failures.append( "correl " + targets[i] )

        # check for slope within range
        f.write( str(round(slope, 3)) + "\t" )
        if ( cutoffs_slope_dict[targets[i]]+2.0 < slope or cutoffs_slope_dict[targets[i]]-2.0 > slope ):
            failures.append( "slope " + targets[i] )

        # check for intercept within range
        f.write( str(round(intercept, 3)) + "\t" )
        if ( cutoffs_intercept_dict[targets[i]]+3.0 < intercept or cutoffs_intercept_dict[targets[i]]-3.0 > intercept ):
            failures.append( "intercept " + targets[i] )

        results.update( {targets[i] : [pearson_correl, slope, intercept] } )
        f.write( "\n" )
    except:
        print("ERROR when processing", ddG_files[i])

f.close()

benchmark.save_variables('debug targets nstruct working_dir testname results ddG_files cutoffs_corr_dict cutoffs_slope_dict cutoffs_intercept_dict failures')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
