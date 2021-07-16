#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available
# (c) under license.
# (c) The Rosetta software is developed by the contributing members of the
# (c) Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org.
# (c) Questions about this can be addressed to
# (c) University of Washington CoMotion, email: license@uw.edu.

## @file   glycan_dock/3.plot.py
## @brief  This script is part of the GlycanDock scientific test.
## @author Sergey Lyskov
## @author Morgan Nance (@mlnance; revised for GlycanDock sci test)

import os, sys, subprocess, math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import benchmark
from benchmark.util import scorefile_io_utils as sciu


benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

# inputs are header labels from the scorefile to plot
# for instance "total_score" and "rmsd"
x_label = "ring_Lrmsd"
y_label = "interaction_energy"
plot_cutoffs = "plot_cutoffs"
outfile = "plot_results.png"

# read plot cutoffs (ymin)
protein = subprocess.getoutput( "grep -v '#' " + plot_cutoffs + " | awk '{print $1}'" ).splitlines()
y_mins = subprocess.getoutput( "grep -v '#' " + plot_cutoffs + " | awk '{print $2}'" ).splitlines()
plot_cutoffs_dict = dict( zip ( protein, [float(y) for y in y_mins] ) )

# number of subplots
ncols = 5
nrows = 1
if len( targets ) < ncols:
        ncols = len( targets )
else:
        nrows = math.ceil( len( targets ) / ncols )

# figure size
width = 7.5 * ncols
height = 6 * nrows

plt.rc("font", size=20)
plt.rcParams['figure.figsize'] = width, height #width, height

# go through scorefiles
for sc_f, ii, target in zip(scorefiles, range(len(scorefiles)), targets):

        # import scorefile as a pandas dataframe
        df = sciu.scorefile_to_dataframe(sc_f).sort_values(y_label)

        # create subplot
        plt.subplot( nrows, ncols, ii+1 )

        # x and y labels
        plt.xlabel( x_label + " (Ang)" )
        plt.ylabel( y_label + " (REU)" )

        # title
        tag = "<N5>"
        plt.title( target + " " + tag + ": " + str(results[target][tag]["avg"]) + " (" + str(results[target][tag]["std"]) + ")" )

        # scatterplot of the data
        plt.plot(df[x_label], df[y_label], 'ko')

        # add veritical line for the near-native cutoff of 2.0
        plt.axvline(x=2.0, color='r', linestyle=':', alpha=0.7)

        # add veritical line for the target's ring_Lrmsd cutoff
        plt.axvline(x=float(cutoffs_ring_Lrmsd_dict[target]), color='b', linestyle='-')
        # and add horizontal line for the target's top-5 interaction_energy
        plt.axhline(y=float(df[y_label].iloc[5]), color='b', linestyle='-')

        plt.xlim( 0, 15 )
        plt.ylim( plot_cutoffs_dict[target], 0 )
	
#save figure
plt.tight_layout()
plt.savefig( outfile )

benchmark.save_variables('debug targets nstruct working_dir testname results outfile cutoffs_ring_Lrmsd_dict failures')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
