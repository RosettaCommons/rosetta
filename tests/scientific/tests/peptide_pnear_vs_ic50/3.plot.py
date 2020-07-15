#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  peptide_pnear_vs_ic50/3.plot.py
## @brief This script is part of peptide_pnear_vs_ic50 scientific test
## @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

import os, sys, subprocess, math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import benchmark
from matplotlib.ticker import MultipleLocator
import numpy as np

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

# Some global settings for MatPlotLib:
font = {'family' : 'Helvetica',
        'weight' : 'regular',
        'size'   : 6
        }
matplotlib.rc( 'font', **font )
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['axes.linewidth'] = 0.25
matplotlib.rcParams['mathtext.default'] = 'regular'

# Draw one energy vs. RMSD plot:
def plot_graph( evals, rmsdvals, hbonds, pnear, delta_g_folding, fig, axs, i, j ) :
    fontbold = {'family' : 'Helvetica', 'weight' : 'bold', 'size' : 6 }
    fontsmall = {'family' : 'Helvetica', 'weight' : 'regular', 'size' : 5 }
    plotmarker = matplotlib.markers.MarkerStyle( marker="o", fillstyle="full" )
    scalefactor_x = 1.0
    scalefactor_y = 1.0

    cdict = {'red':  ((0.0, 0.5, 0.5),
                  (0.25, 0.0, 0.0),
                  (0.5, 0.0, 0.0),
                  (0.75, 1.0, 1.0),
                  (1.0, 1.0, 1.0)),

        'green': ((0.0, 0.0, 0.0),
                  (0.25, 0.0, 0.0),
                  (0.35, 0.75, 0.75),
                  (0.5, 1.0, 1.0),
                  (0.75, 0.7, 0.7),
                  (1.0, 0.0, 0.0)),

        'blue':  ((0.0, 0.5, 0.5),
                  (0.25, 1.0, 1.0),
                  (0.35, 1.0, 1.0),
                  (0.5, 0.0, 0.0),
                  (1.0, 0.0, 0.0))
       }
    rainbow = matplotlib.colors.LinearSegmentedColormap( "rainbow", cdict )

    ax = axs[i,j]
    if( j == 0 ) :
        ax.set_xlabel( "RMSD to design (Å)", labelpad=0.65, fontdict=fontbold )
    else :
        ax.set_xlabel( "RMSD to lowest-energy state (Å)", labelpad=0.65, fontdict=fontbold )
    ax.set_ylabel( "Energy (kcal/mol)", labelpad=0.65, fontdict=fontbold )
    ax.tick_params( axis='x', pad=0.65 )
    ax.tick_params( axis='y', pad=0.15 )
    ax.tick_params( width=0.25 )
    ax.tick_params( length=2 )
    sc = ax.scatter( rmsdvals, evals, marker=plotmarker, alpha=1.0, s=1.0, edgecolors="none", c=hbonds, cmap=rainbow, antialiased=True, rasterized=True )
    ax.set_xlim(0,3.0)
    minval = min( evals ) - 6.0
    maxval = minval + 30.0
    ax.set_ylim(minval, maxval)

    cb = fig.colorbar( sc, ax=ax, ticks=[0,1,2,3,4,5,6,7,8] )
    cb.ax.tick_params( axis='y', pad=0.15 )
    cb.ax.tick_params( width=0.25 )
    cb.ax.tick_params( length=2 )
    cb.ax.set_ylabel( "Hydrogen bonds", rotation=-90, va="bottom", labelpad=0.8, fontdict=fontbold )

    box2 = cb.ax.get_position()
    orig_height2 = box2.y1 - box2.y0
    orig_width2 = box2.x1 - box2.x0
    final_height2 = orig_height2 * scalefactor_y
    delta_y0_2 = (orig_height2 - final_height2)/2.0
    cb.ax.set_position( [box2.x0, box2.y0+delta_y0_2, orig_width2, final_height2 ] )
    
    box1 = ax.get_position()
    orig_height = box1.y1 - box1.y0
    orig_width = box1.x1 - box1.x0
    final_height = orig_height * scalefactor_y
    final_width = orig_width * scalefactor_x
    delta_y0 = (orig_height - final_height)/2.0
    delta_x0 = (orig_width - final_width)/2.0
    ax.set_position( [box1.x0 + delta_x0, box1.y0 + delta_y0, final_width, final_height ] )

    # Write PNear on chart
    pnearstring = "={pnear:.4f}, "
    dgfoldingstring = "={delta_g_folding:.4f}"
    textstring  = r"$P_{Near}$" + pnearstring.format(pnear=pnear) + r"$\Delta$" + r"$G_{folding}$" + dgfoldingstring.format(delta_g_folding=delta_g_folding)
    plt.text( .05, .05, textstring, transform=ax.transAxes, fontdict=fontsmall )

def plot_one_dataset( fig, axs, suffix='1A', row=1 ) :

    # inputs are header labels from the scorefile to plot, for instance "total_score" and "rmsd"
    # => it figures out the column numbers from there
    x_label = "RMSD from known native state (A)"
    x_label2 = "RMSD from lowest energy state found (A)"
    y_label = "Rosetta energy (kcal/mol)"
    outfile = f"plot_results_{suffix}.png"
    outfile2 = f"plot_results2_{suffix}.png"
    plottitle = f"Computed energy landscape for NDM1i_{suffix}"

    # Read in data:
    if( os.path.exists( f'{working_dir}/hpc-logs_{suffix}/.hpc.{testname}_NDM1i_{suffix}.output.0.log' ) == True ):
        logfile = f'{working_dir}/hpc-logs_{suffix}/.hpc.{testname}_NDM1i_{suffix}.output.0.log'
    else:
        logfile = f'{working_dir}/hpc-logs_{suffix}/.hpc.{testname}_NDM1i_{suffix}.log'
    print ( "Reading data from " + logfile + "." )
    rmsd_vals = [ float(i) for i in str( subprocess.getoutput( "grep MPI_worker " + logfile + " -A 1000000 | tail -n+2 | awk '{if( NF == 13 ) {print $3} }'" ) ).split() ]
    rmsd_vals_to_lowest = [ float(i) for i in str( subprocess.getoutput( "grep MPI_worker " + logfile + " -A 1000000 | tail -n+2 | awk '{if( NF == 13 ) {print $4} }'" ) ).split() ]
    hbonds = [ float(i) for i in str( subprocess.getoutput( "grep MPI_worker " + logfile + " -A 1000000 | tail -n+2 | awk '{if( NF == 13 ) {print $6} }'" ) ).split() ]
    energy_vals = [ float(i) for i in str( subprocess.getoutput( "grep MPI_worker " + logfile + " -A 1000000 | tail -n+2 | awk '{if( NF == 13 ) {print $5} }'" ) ).split() ]
    pnear_val = float( subprocess.getoutput( "grep 'PNear:' " + logfile + " | awk '{print $2}'" ) )
    pnear_val_to_lowest = float( subprocess.getoutput( "grep 'PNearLowest:' " + logfile + " | awk '{print $2}'" ) )
    delta_g_folding_val = float( str(subprocess.getoutput( "grep 'kB\*T\*ln(Keq):' " + logfile + " | awk '{print $2}'" ) ) )
    delta_g_folding_val_to_lowest = float( str(subprocess.getoutput( "grep 'kB\*T\*ln(KeqLowest):' " + logfile + " | awk '{print $2}'" ) ) )
    minenergy = min(energy_vals)
    maxrms = max( rmsd_vals )
    maxrms_to_lowest = max( rmsd_vals_to_lowest )
    max_overall_rms = max( maxrms, maxrms_to_lowest )

    # Plot the graph:
    print ( "Plotting " + str(len(rmsd_vals)) + " points." )
    plot_graph( energy_vals, rmsd_vals, hbonds, pnear_val, delta_g_folding_val, fig, axs, row, 0 )
    plot_graph( energy_vals, rmsd_vals_to_lowest, hbonds, pnear_val_to_lowest, delta_g_folding_val_to_lowest, fig, axs, row, 1 )

     # Write peptide name at left:
    ax = axs[row, 0]
    fontbold = {'family' : 'Helvetica', 'weight' : 'bold', 'size' : 8 }
    plt.text( -0.25, .5, "NDM1i_" + suffix, transform=ax.transAxes, fontdict=fontbold, ha="center", va="center", rotation=90 )

def plot_fitcurve_figure( fig, ax, Afitval, kfitval, Rsqfitval, yvals, yerr, all_DG_foldings, suffices ) :
    fontbold = {'family' : 'Helvetica', 'weight' : 'bold', 'size' : 6 }
    fontsmall = {'family' : 'Helvetica', 'weight' : 'regular', 'size' : 5 }

    xvals = np.zeros( len(all_DG_foldings) )
    for i in range( len( all_DG_foldings) ) :
        xvals[i] = all_DG_foldings[i][1]
    plotmarker = matplotlib.markers.MarkerStyle( marker="o", fillstyle="full" )

    ax.set_xlabel( r'$\Delta$' + r'$G_{folding}$' + ' (μM)', labelpad=0.65, fontdict=fontbold )
    ax.set_ylabel( r'$IC_{50}$' + ' (μM)', labelpad=0.65, fontdict=fontbold )
    ax.set_yscale('log')
    sc = ax.scatter( xvals, yvals, marker=plotmarker, edgecolors="none", s=12, c='blue' )
    eb = ax.errorbar( xvals, yvals, yerr=yerr, linestyle="None", elinewidth=0.25, c='black', capsize=1.5, capthick=0.25  )

    #Construct the best fit curve
    xprime = []
    interpolation_points = 100
    for j in range(len(xvals) - 1) :
        for k in range(interpolation_points) :
            xprime.append( float(k)/float(interpolation_points) * ( xvals[j+1] - xvals[j] ) + xvals[j] )
    xprime.append( xvals[len(xvals)-1] )
    yprime = []
    for x in xprime :
            yprime.append( math.exp( kfitval*x + Afitval ) )
    #print( "xprime:\t", xprime)
    #print( "yprime:\t", yprime)

    ax.plot( xprime, yprime, 'blue', linewidth=0.5 )

    # Write R^2 and equation on chart:
    eqnstring = "ln(Y) = {kval:.4f}X + {Aval:.4f}"
    r2string = " = {r2val:.4f}"
    textstring = eqnstring.format( kval=kfitval, Aval=Afitval ) + "\n" + r"$R^2$" + r2string.format( r2val=Rsqfitval )
    plt.text( .05, .975, textstring, transform=ax.transAxes, fontdict=fontsmall, ha='left', va='top' )

    #Labels for points:
    fontsmall_blue = {'family' : 'Arial',
        'weight' : 'regular',
        'size'   : 5,
        'color' : 'blue'
        }
    for i in range(len(all_DG_foldings)) :
        label = "NDM1i_" + suffices[i]
        x = (all_DG_foldings[i][1] - ax.get_xlim()[0] ) / ( ax.get_xlim()[1] - ax.get_xlim()[0] )
        y = (math.log(yvals[i]) - math.log( ax.get_ylim()[0] ) ) / ( math.log(ax.get_ylim()[1]) - math.log(ax.get_ylim()[0]) )
        #print( x, y, label )
        # this method is called for each point
        plt.text( x, y - 0.025, label, transform=ax.transAxes, fontdict=fontsmall_blue, ha='center', va='top')

# Information for each dataset:
suffices = ["1A", "1B", "1C", "1D", "1E", "1F", "1G"]

#Make the binding plots figure:
fig, axs = plt.subplots( nrows=len(suffices), ncols=2, gridspec_kw={ 'width_ratios': [1,1] } )
fig.set_size_inches( 6.5, len(suffices)*1.75 )
fig.subplots_adjust(bottom=0.025, top=0.975, hspace=0.3, left=0.1, right=0.975, wspace=0.2)

for i in range( len(suffices) ) :
    plot_one_dataset( fig, axs, suffix=suffices[i], row=i )

plt.savefig( "e_vs_rmsd_plots.png", dpi=600 )

# Make the fit curve figure:
fig2, axs2 = plt.subplots( nrows=1, ncols=1 )
fitted_A = bundle[0]
fitted_k = bundle[1]
fitted_Rsq = bundle[2]
print( "A=", fitted_A, ", k=", fitted_k, ", R^2=", fitted_Rsq)
fig.set_size_inches( 2.5, 2.5 )
fig.subplots_adjust(bottom=0.025, top=0.975, hspace=0.3, left=0.025, right=0.975, wspace=0.2)
plot_fitcurve_figure( fig2, axs2, fitted_A, fitted_k, fitted_Rsq, fitting_ic50_vals, fitting_ic50_errs, all_DG_foldings, suffices )
plt.savefig( "fitcurve.png", dpi=600 )

benchmark.save_variables('working_dir testname bundle outfile enough_sampling pnear_good pnear_to_lowest_good lowest_E_close_enough sampling_under_0_25_A sampling_beyond_1_5_A sampling_beyond_2_6_A big_energy_gap overall_pass')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
