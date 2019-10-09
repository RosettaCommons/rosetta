#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  stepwise_RNA_favorites/2.analyze.py
## @brief this script is part of stepwise_RNA_favorites scientific test
## @author Andy Watkins

import os, sys, subprocess, math
import numpy as np
import benchmark

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into     s
config = benchmark.config()

results = {}
scorefiles = []
#logfiles = []

# TODO: populate me for each target.
min_e_cutoffs = {'gagua_pentaloop': -1, 'gcaa_tetraloop': -3, 'gg_mismatch': 0, 'j44a_p4p6': -10, 'r2_4x4': -12, 'srl_fixed': -9, 'srl_free_bulgedG': -9, 'srp_domainIV': -10, 'srp_domainIV_fixed': -9, 'tandem_ga_imino': -10, 'tandem_ga_sheared': -12, 'uucg_tetraloop': 5 }
min_e_rmsd_cutoffs = { 'gagua_pentaloop': 3, 'gcaa_tetraloop': 2, 'gg_mismatch': 2, 'j44a_p4p6': 2, 'r2_4x4': 4, 'srl_fixed': 5.5, 'srl_free_bulgedG': 6, 'srp_domainIV': 4, 'srp_domainIV_fixed': 3, 'tandem_ga_imino': 3, 'tandem_ga_sheared': 1, 'uucg_tetraloop': 5 }

def main(args):
    # inputs are header labels from the scorefile, for instance "score" and "rmsd"
    # => it figures out the column numbers from there
    x_label = "rms"
    y_label = "score"

    # scorefiles and logfiles
    scorefiles.extend( [ f'{working_dir}/output/{t}/{t}.out' for t in targets ] )
    #logfiles.extend( [ f'{working_dir}/hpc-logs/hpc.{testname}-{t}.*.log' for t in targets ] )

    # go through scorefiles of targets
    for i in range( 0, len( scorefiles ) ):

        target_results = {}
        
        # get column numbers from labels, 1-indexed
        # We do this in the inner loop because it's fast, and also because it's possible for the index 
        # to change because there is a bizarre instability in how the numbers get set.
        x_index = str( subprocess.getoutput( "grep " + x_label + " " + scorefiles[0] ).split().index( x_label ) + 1 )
        y_index = str( subprocess.getoutput( "grep " + y_label + " " + scorefiles[0] ).split().index( y_label ) + 1 )


        # read in score file, scores are sorted, first one is lowest
        x = subprocess.getoutput( "grep \"^SCORE:\" " + scorefiles[i] + " | grep -v SEQUENCE | grep -v " + x_label + " | sort -nk2 | awk '{print $" + x_index + "}'" ).splitlines()
        y = subprocess.getoutput( "grep \"^SCORE:\" " + scorefiles[i] + " | grep -v SEQUENCE | grep -v " + y_label + " | sort -nk2 | awk '{print $" + y_index + "}'" ).splitlines()

        # map values to floats (were strings)
        # skip any failures??
        newx, newy = [], []
        for x_, y_ in zip(x,y):
            try:
                fx_, fy_ = float(x_), float(y_)
                newx.append(fx_)
                newy.append(fy_)
            except:
                continue
        x, y = newx, newy

        # check for lowest score below cutoff
        print (targets[i], "\t", end=""),
        val_cutoff = check_min_e( y, min_e_cutoffs[targets[i]] )
        target_results.update( val_cutoff )
        target_results.update( {"Min E:": min(y)})

        # check lowest scoring model has low RMSD
        print (targets[i], "\t", end=""),
        val_topscoring = check_rmsd_of_topscoring( x, min_e_rmsd_cutoffs[targets[i]] )
        target_results.update( val_topscoring )
        target_results.update( {"Min RMSD of 400 topscoring:": min(x[:400])})

        results.update( {targets[i] : target_results} )
        print ("\n")

    benchmark.save_variables('targets nstruct working_dir testname results scorefiles min_e_cutoffs min_e_rmsd_cutoffs')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)

#=======================================
def check_rmsd_of_topscoring( rmsd_col_sorted, cutoff ):
    """
    400 top scoring is 5%
    """

    out = "min rmsd of 400 topscoring models below " + str( cutoff )
    print (out, "\t", end="")

    if min(rmsd_col_sorted[:400]) <= cutoff:
        value = "TRUE"
    else:
        value = "FALSE"

    print (value)
    return {out : value}

#=======================================
def check_min_e( energy_col_sorted, cutoff ):

    out = "energy of low scoring model below " + str( cutoff )
    print (out, "\t", end="")

    if energy_col_sorted[0] <= cutoff:
        value = "TRUE"
    else:
        value = "FALSE"

    print (value)
    return {out : value}


if __name__ == "__main__": main(sys.argv)
