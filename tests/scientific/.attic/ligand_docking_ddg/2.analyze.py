#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  cartesian_relax/2.analyze.py
## @brief this script is part of cartesian_relax scientific test
## @author Sergey Lyskov
## @author Shannon Smith

import os, sys, subprocess, math
import numpy as np
import pandas as pd
import scipy.stats
import benchmark

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into 	s
config = benchmark.config()

#scorefiles = []
#logfiles = []
exp_file = f'{rosetta_dir}/tests/scientific/data/{testname}/experimental_values.txt'

#=======================================
def import_scorefile(input_scorefile):
	selection = pd.read_csv(input_scorefile, delim_whitespace=True, skiprows=1)[['ligand_rms_no_super_X', 'interface_delta_X']]
	return selection

#=======================================
def get_exp_val(input_file, pdb_name, ligand_name ):
    if os.path.exists(input_file):
        out = "log of experimental binding affinity "
        selection = pd.read_csv(input_file, delim_whitespace=True)[['PDB','LIG','pEXP']]
        for index, line in selection.iterrows():
            if line[0] == str(pdb_name) and line[1]== (str(ligand_name)):
                value = '%0.2f' % line[2]
                return { out : value }
            else:
                continue
    else:
        print("file not here")
        exit(1)

#=======================================
def check_all_values_below_cutoff( rmsd_col, cutoff, tag ):

    out = "All " + tag + "s < " + str( cutoff )
    #print (out, end=""),

    if all( i <= cutoff for i in rmsd_col ):
        value = "TRUE"
    else:
        value = "FALSE"

    #print (value)
    return {out : value}

#=======================================
def check_rmsd_of_topscoring( rmsd_col_sorted, cutoff ):

    out = "ligand_rms_no_super_X of topscoring model below " + str( cutoff )
    #print (out, "\t", end="")

    if rmsd_col_sorted[0] <= cutoff:
        value = "TRUE"
    else:
        value = "FALSE"

    #print (value)

    return {out : value}

#=======================================
def check_range( col, tag ):

    #print (tag, "\tmin, max, avg, std:", '% 12.3f % 12.3f % 12.3f % 12.3f' % (min( col ), max( col ), np.mean( col ), np.std( col )))
    value = { "min" : round(min( col ), 4), "max" : round(max( col ), 4), "avg" : round(np.mean( col ), 4), "std" : round(np.std( col ), 4) }
    return {tag : value}

#=======================================
# experimental score versus predicted score
def spearman_rank_correlation( exp , pred ):
    exp_ranks = pd.Series(exp).rank()
    pred_ranks = pd.Series(pred).rank()
    rho = float(str(scipy.stats.spearmanr( exp_ranks, pred_ranks )).split("correlation=")[1].split(',')[0]) ## scipy.stats.spearmanr automatically prints out as SpearmanR(rho, p-value) to 984498 digits. sorry this is annoying
    p_value = float(str(scipy.stats.spearmanr( exp_ranks, pred_ranks )).split("pvalue=")[1].split(')')[0])
    return '%0.2f %0.2f' % (rho , p_value)
#=======================================

results = {}
# inputs are header labels from the scorefile, for instance "total_score" and "rmsd"
# => it figures out the column numbers from there
for t in targets:
    target_results = {}
    for l in ligands:
        target_ligand_results = {}
        scorefile = f'{working_dir}/output/{t}/{t}_{l}.sc'
        #scorefiles.append()
        if not os.path.exists(scorefile):
            continue
        else:
            # Can't figure out a way to separate out ligands for the same target--probably something simple, but oh well. Just have to iterate through possible target-ligand pairs and if it exists, then run this.
            # Having multiple ligands for the same target made nesting these more complicated than I was expecting....
            target_ligand_results.update({'ligand' : l })
            # print(results)
            # get column numbers from labels, 1-indexed
            x_label = "ligand_rms_no_super_X"
            y_label = "interface_delta_X"
            x_index = str(subprocess.getoutput("grep " + x_label + " " + scorefile).split().index(x_label) + 1)
            y_index = str(subprocess.getoutput("grep " + y_label + " " + scorefile).split().index(y_label) + 1)

            # go through scorefiles of targets
            # read in score file, scores are sorted, first one is lowest
            x = subprocess.getoutput( "grep -v SEQUENCE " + scorefile + " | grep -v " + x_label + " | sort -nk2 | awk '{print $" + x_index + "}'" ).splitlines()
            y = subprocess.getoutput( "grep -v SEQUENCE " + scorefile + " | grep -v " + y_label + " | sort -nk2 | awk '{print $" + y_index + "}'" ).splitlines()

            # map values to floats (were strings)
            x = list( map( float, x ))
            y = list( map( float, y ))

            # check for RMSDs below cutoff (native aligned pose in this case)
            #print (targets[i], "\t", end=""),
            val_cutoff = check_all_values_below_cutoff( x, 2, "ligand_rms_no_super_X" )
            target_ligand_results.update( val_cutoff )

            # check lowest scoring model has low RMSD
            #print (targets[i], "\t", end=""),
            val_topscoring = check_rmsd_of_topscoring( x, 2 )
            target_ligand_results.update( val_topscoring )

            # check for RMSD range
            #print (targets[i], "\t", end=""),
            val_rms = check_range( x, "ligand_rms_no_super_X" )
            target_ligand_results.update( val_rms )

            # check for score range
            #print (targets[i], "\t", end=""),
            val_score = check_range( y, "interface_delta_X" )
            target_ligand_results.update( val_score )

            # get predicted value for pdb-ligand pair
            #print (targets[i], "\t", end=""),
            exp_val = get_exp_val(exp_file, t, l)
            target_ligand_results.update( exp_val )

            #print(target_ligand_results)
            #target_results.update({ l : target_ligand_results })
        target_results.update({l : target_ligand_results})

        target_spearman = spearman_rank_correlation(list({t: l['log of experimental binding affinity ']}) , list({t : l['interface_delta_X']['avg'] }) )
    results.update({t: target_results})
    #print(results)
    benchmark.save_variables('working_dir testname results spearman')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)