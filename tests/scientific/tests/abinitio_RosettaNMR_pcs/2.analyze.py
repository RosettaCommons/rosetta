#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  abinitio_RosettaNMR_pcs/2.analyze.py
## @brief this script is part of abinitio_RosettaNMR_pcs scientific test
## @author Sergey Lyskov

import os, sys, subprocess, math
import numpy as np
import benchmark
from benchmark.util import quality_measures as qm

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into     s
config = benchmark.config()

results = {}
scorefiles = []
#logfiles = []
cutoffs_rmsd_dict = {}
cutoffs_score_dict = {}
cutoffs_pnear_dict = {}
failures = []

# inputs are header labels from the scorefile, for instance "total_score" and "rmsd"
# => it figures out the column numbers from there
x_label = "rms"
y_label = "total_score"
outfile = "result.txt"
cutoffs = "cutoffs"

# scorefiles and logfiles
scorefiles.extend( [ f'{working_dir}/output/{t}/{t}.score' for t in targets ] )
#logfiles.extend( [ f'{working_dir}/hpc-logs/hpc.{testname}-{t}.*.log' for t in targets ] )

# get column numbers from labels, 1-indexed
x_index = str( subprocess.getoutput( "grep " + x_label + " " + scorefiles[0] ).split().index( x_label ) + 1 )
y_index = str( subprocess.getoutput( "grep " + y_label + " " + scorefiles[0] ).split().index( y_label ) + 1 )

# get number of fields in scorefile
nfields = len( subprocess.getoutput( "grep -v SEQUENCE " + scorefiles[0] + " | grep " + y_label + " | head -n1" ).split())

# read cutoffs
protein = subprocess.getoutput( "grep -v '#' " + cutoffs + " | awk '{print $1}'" ).splitlines()
cutoffs_rmsd = subprocess.getoutput( "grep -v '#' " + cutoffs + " | awk '{print $2}'" ).splitlines()
cutoffs_score = subprocess.getoutput( "grep -v '#' " + cutoffs + " | awk '{print $3}'" ).splitlines()
cutoffs_pnear = subprocess.getoutput( "grep -v '#' " + cutoffs + " | awk '{print $4}'" ).splitlines()
cutoffs_rmsd = map( float, cutoffs_rmsd )
cutoffs_score = map( float, cutoffs_score )
cutoffs_pnear = map( float, cutoffs_pnear )
cutoffs_rmsd_dict.update( dict( zip ( protein, cutoffs_rmsd )))
cutoffs_score_dict.update( dict( zip ( protein, cutoffs_score )))
cutoffs_pnear_dict.update( dict( zip ( protein, cutoffs_pnear )))
Pnear = {}

lambdas = {"1x0n": 8.0, "2ae9": 4.0, "2k61":6.0}

# open results output file
f = open( outfile, "w" )

# go through scorefiles of targets
for i in range( 0, len( scorefiles ) ):

    target_results = {}

    # read in score file, scores are sorted, first one is lowest
    x = subprocess.getoutput( "awk '{if(NF==" + str(nfields) + ") print}' " + scorefiles[i] + " | grep -v SEQUENCE | grep -v " + y_label + " | sort -nk2 | awk '{print $" + x_index + "}'" ).splitlines()
    y = subprocess.getoutput( "awk '{if(NF==" + str(nfields) + ") print}' " + scorefiles[i] + " | grep -v SEQUENCE | grep -v " + y_label + " | sort -nk2 | awk '{print $" + y_index + "}'" ).splitlines()

    # map values to floats (were strings)
    x = list( map( float, x ))
    y = list( map( float, y ))

    # check for RMSDs below cutoff
    f.write( targets[i] + "\t" )
    val_cutoff = qm.check_xpercent_values_below_cutoff( x, cutoffs_rmsd_dict[targets[i]], "rmsd", f, 1 )
    target_results.update( val_cutoff )

    # add to failues
    if val_cutoff['All rmsds < cutoff'] == False:
        failures.append( targets[i] )

    # check for scores below cutoff
    f.write( targets[i] + "\t" )
    val_cutoff = qm.check_xpercent_values_below_cutoff( y, cutoffs_score_dict[targets[i]], "score", f, 1 )
    target_results.update( val_cutoff )

    # add to failures
    if val_cutoff['All scores < cutoff'] == False:
        failures.append( targets[i] )

    # compute Pnear
    f.write( targets[i] + "\t" )
    val_cutoff = round( float(qm.calculate_pnear( y, x, lambda_val=lambdas[targets[i]] )), 4)
    target_results.update( {"Pnear <= cutoff": val_cutoff} )
    Pnear[targets[i]] = val_cutoff    

    # add to failures
    if val_cutoff < cutoffs_pnear_dict[targets[i]]:
        failures.append( targets[i] )

    # check lowest scoring model has low RMSD
#        f.write( targets[i] + "\t" )
#        val_topscoring = check_rmsd_of_topscoring( x, cutoffs_rmsd_dict[targets[i]], f )
#        target_results.update( val_topscoring )

    # check for RMSD range
    f.write( targets[i] + "\t" )
    val_rms = qm.check_range( x, "rmsd", f )
    target_results.update( val_rms )

    # check for score range
    f.write( targets[i] + "\t" )
    val_score = qm.check_range( y, "score", f )
    target_results.update( val_score )

    results.update( {targets[i] : target_results} )
    f.write( "\n" )

f.close()

benchmark.save_variables('debug targets nstruct working_dir testname results scorefiles cutoffs_rmsd_dict cutoffs_score_dict failures Pnear')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
