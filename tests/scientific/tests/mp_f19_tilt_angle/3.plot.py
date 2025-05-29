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
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
import energy_landscape_metrics
import numpy as np
from scipy import interpolate
import benchmark

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

# Output file for results
outfile = "plot_results.png"

#number of subplots
ncols = 5
nrows = 1

# figure size
width = 7.5 * ncols
height = 6 * nrows

plt.rc("font", size=20)
plt.rcParams['figure.figsize'] = width, height #width, height

# read cutoffs
cutoffs = "cutoffs"
protein = subprocess.getoutput( "grep -v '#' " + cutoffs + " | awk '{print $1}'" ).splitlines()
cutoffs_zmin = subprocess.getoutput( "grep -v '#' " + cutoffs + " | awk '{print $3}'" ).splitlines()
cutoffs_anglemin = subprocess.getoutput( "grep -v '#' " + cutoffs + " | awk '{print $4}'" ).splitlines()
cutoffs_zexp = subprocess.getoutput( "grep -v '#' " + cutoffs + " | awk '{print $5}'" ).splitlines()
cutoffs_angleexp = subprocess.getoutput( "grep -v '#' " + cutoffs + " | awk '{print $6}'" ).splitlines()

cutoffs_zmin = map( float, cutoffs_zmin )
cutoffs_anglemin = map( float, cutoffs_anglemin )
cutoffs_zexp = map( float, cutoffs_zexp )
cutoffs_angleexp = map( float, cutoffs_angleexp )

cutoffs_zmin_dict.update( dict( zip ( protein, cutoffs_zmin )))
cutoffs_anglemin_dict.update( dict( zip( protein, cutoffs_anglemin )))
cutoffs_zexp_dict.update( dict( zip( protein, cutoffs_zexp )))
cutoffs_angleexp_dict.update( dict( zip( protein, cutoffs_angleexp )))

# go through energy landscapes
for i in range( 0, len( energy_landscapes ) ):

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
        angles.append( float( contents[x][1] ) )
        total_scores.append( float( contents[x][3] ) )

    zcoords_arr = np.asarray( zcoords )
    angles_arr = np.asarray( angles )
    total_scores_arr = np.asarray( total_scores )

    # create subplot
    ax = plt.subplot( nrows, ncols, i+1 )

    # x and y labels
    plt.xlabel( "Angle (degrees)" )
    plt.ylabel( "Membrane Depth (Angstroms)" )
    plt.title( targets[i] )

    # energy landscape plot of the data
    X, Y = np.mgrid[ angles_arr.min():angles_arr.max(), zcoords_arr.min():zcoords_arr.max()]
    points = np.c_[X.ravel(), Y.ravel()]
    try:
        Z = interpolate.griddata(np.c_[angles_arr, zcoords_arr], total_scores_arr, points).reshape(X.shape)
    except:
        print("ERROR - unable to interpolate data for", energy_landscapes[i])
        continue

    #RS: adding contours to have a sense of error bars in the energy landscape.
    contours = plt.contour( X, Y, Z, 10, colors='black')
    plt.clabel(contours, inline=True, fontsize=8, colors='red')

    cmap = cm.get_cmap('viridis', 256)
    im = ax.pcolormesh( X,Y,Z, cmap=cmap)
    cbar = plt.colorbar(im)

    # Calculate the minimum energy point
    zmin, anglemin = energy_landscape_metrics.compute_minimum_energy_orientation( zcoords_arr, angles_arr, total_scores_arr )
    #changing the minimum angle to first phase; results are symmetrical about the z-axis which is 0 and 180 degrees.

    if( anglemin>90 and anglemin<=180 ):
        anglemin = 180 - anglemin
    elif( anglemin>180 and anglemin<=270 ):
        anglemin = anglemin - 180
    elif( anglemin>270 and anglemin<=360 ):
        anglemin = 360 - anglemin
    plt.plot( anglemin, zmin, 'wo', markersize=12 )

    # Plot experimental values as red open squares
    plt.plot( cutoffs_angleexp_dict[targets[i]], cutoffs_zexp_dict[targets[i]], '^r', fillstyle='none', markeredgewidth=5, markersize=20 )

    #changing the value of ticks in the x-axis
    plt.xticks(np.arange(0, 180, 45))

#save figure
plt.tight_layout()
plt.savefig( outfile )

benchmark.save_variables('debug targets nstruct working_dir testname results outfile cutoffs_rmsd_dict cutoffs_score_dict failures')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
