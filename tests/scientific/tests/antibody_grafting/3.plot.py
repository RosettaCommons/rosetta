#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  cartesian_relax/3.plot.py
## @brief this script is part of cartesian_relax scientific test
## @author Sergey Lyskov

import os, sys, subprocess, math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import benchmark

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

# inputs are header labels from the scorefile to plot, for instance "total_score" and "rmsd"
# => it figures out the column numbers from there
x_label = "backbone rmsd"
y_label = "count"
outfile = "plot_results.png"

# number of subplots is the number of regions
# we have: frh, h1, h2, h3, frl, l1, l2, l3
ncols = 4
nrows = 2

# figure size
width = 7.5 * ncols
height = 6 * nrows

plt.rc("font", size=20)
plt.rcParams['figure.figsize'] = width, height #width, height
plt.rcParams['axes.linewidth'] = 1.5 # thicker borders

# for each region in the results dict, plot the density
i = 0
for k,v in target_results.items():

    # create subplot
    i+=1
    plt.subplot( nrows, ncols, i )

    # label
    plt.xlabel( x_label )
    plt.ylabel( y_label)
    plt.title( k )

    # plot data
    data = target_results[k]
    n_points = len(data)
    plt.hist(data, bins=math.floor(math.sqrt(n_points)), density=False, fill=False, histtype='step', color='k', linewidth=1.5)
    # add points for unbinned data
    plt.plot(data, len(data)*[0.5], '|k', markersize=18, markeredgewidth=2)

    # add vertical line for cutoff (1.0 A)
    plt.axvline(x=cutoffs_rms_dict[k], color='r', linestyle='-', linewidth=1.5)

#save figure
plt.tight_layout()
plt.savefig( outfile )

benchmark.save_variables('targets working_dir testname results outfile target_results cutoffs_rms_dict cutoffs_fraction_dict failures off_by_more_than_X')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
