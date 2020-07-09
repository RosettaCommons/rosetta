#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  simple_cycpep_predict/3.plot.py
## @brief This script is part of simple_cycpep_predict scientific test
## @author Sergey Lyskov.
## @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

import os, sys, subprocess, math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import benchmark
from matplotlib.ticker import MultipleLocator

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

# inputs are header labels from the scorefile to plot, for instance "total_score" and "rmsd"
# => it figures out the column numbers from there
x_label = "RMSD from known native state (A)"
x_label2 = "RMSD from lowest energy state found (A)"
y_label = "Rosetta energy (kcal/mol)"
outfile = "plot_results.png"
outfile2 = "plot_results2.png"
plottitle = "Computed energy landscape"

# Read in data:
if( os.path.exists( f'{working_dir}/hpc-logs/.hpc.{testname}.output.0.log' ) == True ):
	logfile = f'{working_dir}/hpc-logs/.hpc.{testname}.output.0.log'
else:
	logfile = f'{working_dir}/hpc-logs/.hpc.{testname}.log'
print ( "Reading data from " + logfile + "." )
rmsd_vals = [ float(i) for i in str( subprocess.getoutput( "grep MPI_worker " + logfile + " -A 1000000 | tail -n+2 | awk '{if( NF == 13 ) {print $3} }'" ) ).split() ]
rmsd_vals_to_lowest = [ float(i) for i in str( subprocess.getoutput( "grep MPI_worker " + logfile + " -A 1000000 | tail -n+2 | awk '{if( NF == 13 ) {print $4} }'" ) ).split() ]
energy_vals = [ float(i) for i in str( subprocess.getoutput( "grep MPI_worker " + logfile + " -A 1000000 | tail -n+2 | awk '{if( NF == 13 ) {print $5} }'" ) ).split() ]
minenergy = min(energy_vals)
maxrms = max( rmsd_vals )
maxrms_to_lowest = max( rmsd_vals_to_lowest )
max_overall_rms = max( maxrms, maxrms_to_lowest )

# Plot the graph:
print ( "Plotting " + str(len(rmsd_vals)) + " points." )

for f in range(2):
    # figure size
    width = 7.5
    height = 6
    plt.rc("font", size=12)
    plt.rcParams['figure.figsize'] = width, height #width, height

    # add gridlines
    fig, ax = plt.subplots()
    ax.xaxis.set_major_locator( MultipleLocator(0.5) )
    ax.xaxis.set_minor_locator( MultipleLocator(0.1) )
    ax.yaxis.set_major_locator( MultipleLocator(5) )
    ax.yaxis.set_minor_locator( MultipleLocator(1) )
    plt.grid( c='blue', linestyle='-', linewidth=1.0, alpha=0.15, which='major' )
    plt.grid( c='blue', linestyle='-', linewidth=0.5, alpha=0.1, which='minor' )

    # x and y labels
    if( f == 0 ) :
        plt.xlabel( x_label )
    else :
        plt.xlabel( x_label2 )
    plt.ylabel( y_label )
	
    # set title
    plt.title( plottitle )

    # scatterplot of the data
    if( f == 0 ) :
        plt.scatter(rmsd_vals, energy_vals, c='blue', s=15, alpha=0.1 )
    else :
        plt.scatter(rmsd_vals_to_lowest, energy_vals, c='purple', s=15, alpha=0.1 )

    # Y-axis range
    plt.ylim( minenergy - 2.5, minenergy +27.5  )
    plt.xlim( 0, max_overall_rms + 0.1  )
	
    #save figure
    plt.tight_layout()
    if f == 0 :
        plt.savefig( outfile )
    else :
        plt.savefig( outfile2 )

benchmark.save_variables('working_dir testname outfile enough_sampling pnear_good pnear_to_lowest_good lowest_E_close_enough sampling_under_0_25_A sampling_beyond_1_5_A sampling_beyond_2_6_A big_energy_gap overall_pass')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
