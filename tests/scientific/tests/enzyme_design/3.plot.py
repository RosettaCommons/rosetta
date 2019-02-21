#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  enzyme_design/3.plot.py
## @brief this script is part of enzyme_design scientific test
## @author Sergey Lyskov
## @author Rocco Moretti (rmorettiase@gmail.com)

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
plot_outfile = "plot_results.png"
plot_metrics = "seqrec pssm_seqrec pssm_delta_seqrec".split()

# What's the canonical sorted order for the targets (descending by median performance of ~10 reference runs.)
target_order = {
 "seqrec":            "1FZQ 1y52 2HZQ 1sw1 1n4h 1fby 2h6b 2Q89 1XT8 1USK 2qo4 1OPB 2FR3 2ifb 1wdn 1y3n 1hsl 1nq7 2RDE 1db1 2FQX 1XZX 1y2u 1hmr 2e2r 2GM1 2Q2Y 1A99 1x7r 1urg 1H6H 2f5t 2p0d 2UYI 1ZHX 1z17 1uw1 1RBP 2PFY 3B50 1l8b 1J6Z 2ioy 2DRI 1LKE 1TYR 1nl5 1POT 2b3b 2rct 2FME".split(),
 "pssm_seqrec":       "2h6b 1y52 1XT8 1y3n 2Q89 2qo4 1OPB 2PFY 1n4h 1db1 2p0d 1A99 1RBP 2e2r 1fby 2HZQ 1z17 2FR3 1hsl 1USK 2ifb 1nq7 1sw1 1J6Z 1FZQ 2RDE 1TYR 2rct 1H6H 2f5t 1wdn 1x7r 1y2u 1XZX 1LKE 3B50 1nl5 1ZHX 2Q2Y 1uw1 2GM1 1l8b 1hmr 2UYI 2FQX 1urg 1POT 2b3b 2ioy 2DRI 2FME".split(),
 "pssm_delta_seqrec": "2h6b 2qo4 2Q89 1y3n 2ifb 1y52 1LKE 1n4h 1USK 1RBP 1OPB 2e2r 1fby 1XZX 1db1 1XT8 2HZQ 1TYR 2FQX 2FR3 1H6H 2f5t 2PFY 1wdn 1nq7 2p0d 1FZQ 3B50 1z17 1sw1 1hsl 2rct 2RDE 1J6Z 1A99 1ZHX 1y2u 2GM1 1x7r 1hmr 2Q2Y 1uw1 1l8b 1urg 2UYI 1nl5 2ioy 1POT 2DRI 2b3b 2FME".split(),
}

#number of subplots
ncols = 4
nrows = 1
if len( plot_metrics ) < ncols:
	ncols = len( plot_metrics )
else:
	nrows = math.ceil( len( plot_metrics ) / ncols )

# figure size
width = 7.5 * ncols
height = 6 * nrows

plt.rc("font", size=20)
plt.rcParams['figure.figsize'] = width, height #width, height


## get column numbers from labels, 1-indexed
#x_index = str( subprocess.getoutput( "grep " + x_label + " " + scorefiles[0] ).split().index( x_label ) + 1 )
#y_index = str( subprocess.getoutput( "grep " + y_label + " " + scorefiles[0] ).split().index( y_label ) + 1 )

# go through scorefiles
for i, metric in enumerate( plot_metrics ):
    values = []
    targ_labels = []
    for targ in target_order.get( metric, targets ):
        if targ not in results:
            continue
        if metric not in results[targ]:
            continue
        values.append( results[targ][metric] )
        targ_labels.append( targ )

    # map all values to floats
    x = range(len(targ_labels))
    y = list( map( float, values ) )

    # create subplot
    plt.subplot( nrows, ncols, i+1 )

    # y labels
    plt.ylabel( metric )

    # set title
    plt.title( metric )

    # bar plot of the data
    plt.bar(x, y)

    plt.xlim( [ min(x)-1, max(x)+1 ] ) # Defaults give a bit too much space

    # label X-axis
    plt.xticks(x, targ_labels, rotation=90, fontsize="9") # Need a smaller font to fit all the proteins on the chart

    # add horizontal lines for cutoff
    if "ALL" in thresholds and metric in thresholds["ALL"]:
        plt.axhline(y=float(thresholds["ALL"][metric][1]), color='b', linestyle='-')

#save figure
plt.tight_layout()
plt.savefig( plot_outfile )

#benchmark.save_variables('targets nstruct working_dir testname results outfile cutoffs_rmsd_dict cutoffs_score_dict failures')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
benchmark.save_variables('targets nstruct working_dir testname results thresholds plot_outfile')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
