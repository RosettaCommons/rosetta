#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  peptide_pnear_vs_ic50/2.analyze.py
## @brief This script is part of peptide_pnear_vs_ic50 scientific test.
## @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).

import os, sys, subprocess, math
import re
import benchmark
import numpy as np
from scipy.optimize import least_squares
from scipy.optimize import curve_fit
from scipy.stats import linregress

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into     s
config = benchmark.config()

# If bool is true, print YES; otherwise print NO
def bool_to_string( bool_in ):
    assert( bool_in == True or bool_in == False )
    if( bool_in == True ):
        return "YES"
    return "NO"

logfile_names = {}

def analyze_one_funnel( suffix='1A', expected_rmsd_of_lowest=0.3, expected_min_rmsd=0.25, expected_max_rmsd=2.5, expected_min_pnear=0.9, expected_energy_gap=6 ):
    # Things we'll check
    pnear_good = False
    pnear_to_lowest_good = False
    enough_sampling = False
    lowest_E_is_first = False
    lowest_E_close_enough = False
    sampling_under_0_25_A = False
    sampling_beyond_1_5_A = False
    sampling_beyond_2_6_A = False
    big_energy_gap = False
    overall_pass = False

    # Expect at least a certain fraction of attempts to produce a sample:
    samples_expected = 450 if debug else 18000

    # Relevant filenames
    outfile = f'{working_dir}/result_{suffix}.txt'

    for logfile in f'{working_dir}/hpc-logs_{suffix}/.hpc.{testname}_NDM1i_{suffix}.output.0.log {working_dir}/hpc-logs_{suffix}/.hpc.{testname}_NDM1i_{suffix}.log {working_dir}/hpc-logs_{suffix}/.hpc.{testname}_NDM1i_{suffix}.output'.split():
        if os.path.isfile(logfile):
            logfile_names[suffix] = logfile
            break

    # Check that it exists and has contents.
    if not os.path.exists( logfile ):
        raise ValueError( "Logfile `" + logfile + "` does not exist, but should." )
    with open( logfile ) as f:
        logfile_contents = f.readlines()
    if len( ''.join( l.strip() for l in logfile_contents[:10] ) ) == 0:
        #Empty, or at least nothing in the first 10 lines ... problem
        raise ValueError( "Logfile `" + logfile + "` is empty, but it shouldn't be." )

    # read relevant data:

    rmsd_vals = []
    rmsd_vals_to_lowest = []
    energy_vals = []
    pnear = None
    pnear_to_lowest = None
    DG_folding = None
    DG_folding_to_lowest = None

    in_table = False
    for line in logfile_contents:
        try:
            sline = line.split()
            if in_table:
                if line.startswith("End summary"):
                    in_table = False
                    continue
                rmsd_vals.append( float(sline[2]) )
                rmsd_vals_to_lowest.append( float(sline[3]) )
                energy_vals.append( float(sline[4]) )
            elif sline[0] == "MPI_worker_node":
                in_table = True
                continue
            elif sline[0] == "PNear:":
                pnear = float( sline[1] )
            elif sline[0] == "PNearLowest:":
                pnear_to_lowest = float( sline[1] )
            elif sline[0] == "-kB*T*ln(Keq):":
                DG_folding = float( sline[1] )
            elif sline[0] == "-kB*T*ln(KeqLowest):":
                DG_folding_to_lowest = float( sline[1] )
        except:
            # Provide a slightly more useful debugging info for outputs
            print("ISSUE IN", logfile, "can't parse", line)
            raise

    print( "Read PNear=" + str(pnear) + " from " + logfile + "." )
    print( "Read PNearLowest=" + str(pnear_to_lowest) + " from " + logfile + "." )
    print( "Read DG_folding=" + str(DG_folding) + " from " + logfile + "." )
    print( "Read DG_folding_to_lowest=" + str(DG_folding_to_lowest) + " from " + logfile + "." )

    assert None not in [ pnear, pnear_to_lowest, DG_folding, DG_folding_to_lowest ]

    sample_match = re.search("[1-90]+ jobs returned structures", str( subprocess.getoutput( 'grep "application completed" ' + logfile ) ) )
    if sample_match is not None:
        total_samples = int( sample_match.group().split()[0] )
    else:
        print("ERROR: Could not find the number of structures!")
        total_samples = 0

    print ( "Determined that " + str( total_samples ) + " samples were performed." )

    # Check that total samples is high enough.
    if( total_samples > samples_expected ) :
        enough_sampling = True

    # Check that PNear is reasonable:
    if ( pnear > expected_min_pnear ):
        pnear_good = True

    # Check that PNear to loweset energy is reasonable:
    if ( pnear_to_lowest > expected_min_pnear ):
        pnear_to_lowest_good = True

    # Cheack that lowest energy is first entry
    lowestE = min( energy_vals )
    if ( abs( energy_vals[0] - lowestE ) < 1e-8 ) :
        lowest_E_is_first = True

    # Check that lowest energy is below threshold.
    if( rmsd_vals[0] < expected_rmsd_of_lowest ) :
        lowest_E_close_enough = True

    # Check that we're sampling diverse structures
    max_rmsd = max( rmsd_vals )
    print( "Max RMSD = " + str( max_rmsd ) )
    if( max_rmsd > 1.5 ) :
        sampling_beyond_1_5_A = True
    if( max_rmsd > expected_max_rmsd ) :
        sampling_beyond_2_6_A = True

    # Check that we're sampling under 0.25 A
    min_rmsd = min( rmsd_vals )
    print( "Min RMSD = " + str( min_rmsd ) )
    if( min_rmsd < expected_min_rmsd ) :
        sampling_under_0_25_A = True

    # Construct the vector of entries with RMSD > 1.5 A
    energies_over_1_5 = []
    rmsds_over_1_5 = []
    for i in range(0, len(rmsd_vals)):
        if( rmsd_vals[i] > 1.5 ):
            energies_over_1_5.append( energy_vals[i] )
            rmsds_over_1_5.append( rmsd_vals[i] )

    if( len( energies_over_1_5 ) > 0 ) :
        min_energy_over_1_5 = min( energies_over_1_5 )
        energy_gap_val = min_energy_over_1_5 - energy_vals[0]
        print( "Energy gap = " + str(energy_gap_val) )
        if( energy_gap_val > expected_energy_gap ):
            big_energy_gap = True

    # Determine if we passed overall
    if( pnear_good == True and pnear_to_lowest_good == True and enough_sampling and lowest_E_is_first == True and lowest_E_close_enough == True and sampling_under_0_25_A == True and  sampling_beyond_1_5_A == True and sampling_beyond_2_6_A == True and big_energy_gap == True ):
        overall_pass = True

    # Write out results
    with open( outfile, "w" ) as f:
        f.write( "NDM1i-" + suffix + ":\n" )
        f.write( "Total samples =\t" + str(total_samples) + "\n" )
        f.write( "Computed PNear =\t" + "%.4f" % round(pnear,4) + "\n" )
        f.write( "Computed PNear to lowest E =\t" + "%.4f" % round(pnear_to_lowest,4) + "\n" )
        f.write( "Computed DG_folding =\t" + "%.4f" % round(DG_folding,4) + "\n" )
        f.write( "Computed DG_folding to lowest =\t" + "%.4f" % round(DG_folding_to_lowest,4) + "\n" )
        f.write( "Lowest energy =\t" + "%.4f" % round(lowestE,4) + " kcal/mol\n" )
        f.write( "RMSD of lowest energy =\t" + "%.3f" % round(rmsd_vals[0],3) + " Angstroms\n" )
        f.write( "Lowest RMSD =\t" + "%.3f" % round(min_rmsd,3) + " Angstroms\n" )
        f.write( "Highest RMSD =\t" + "%.3f" % round(max_rmsd,3) + " Angstroms\n" )
        f.write( "Energy gap (minE>1.5A - minE) =\t" + "%.4f" % round(energy_gap_val,4) + " kcal/mol\n" )
        f.write( "\n" )

        if( debug == True ):
            f.write( "More than 450 samples?\t" + bool_to_string( enough_sampling ) + "\n" )
        else:
            f.write( "More than 18,000 samples?\t" + bool_to_string( enough_sampling ) + "\n" )
        f.write( "PNear value over " + str(expected_min_pnear) + "?\t" + bool_to_string( pnear_good ) + "\n" )
        f.write( "PNear value to lowest E over " + str(expected_min_pnear) + "?\t" + bool_to_string( pnear_to_lowest_good ) + "\n" )
        f.write( "Lowest energy under " + str( expected_rmsd_of_lowest ) + " A RMSD?\t" + bool_to_string( lowest_E_close_enough ) + "\n" )
        f.write( "Sampling below expected lower threshold RMSD(" + str(expected_min_rmsd) + " A)?\t" + bool_to_string( sampling_under_0_25_A ) + "\n" )
        f.write( "Sampling beyond 1.5 A RMSD?\t" + bool_to_string( sampling_beyond_1_5_A ) + "\n" )
        f.write( "Sampling beyond " + str(expected_max_rmsd) + " A RMSD?\t" + bool_to_string( sampling_beyond_2_6_A ) + "\n" )
        f.write( str(expected_energy_gap) + "+ kcal/mol energy gap?\t" + bool_to_string( big_energy_gap ) + "\n" )
        f.write( "OVERALL PASS?\t" + bool_to_string( overall_pass ) + "\n" )

    print( "Finished analyzing NDM1i-" + suffix + ".  Overall pass = " + bool_to_string( overall_pass ) )
    return ( suffix, DG_folding, DG_folding_to_lowest )

# Read the experimentally-measured binding data from a data file:
def load_ic50_data( filename ) :
    vals = []
    errs = []

    with open( filename, "r" ) as f:
        lines = f.readlines()
    first = True
    for line in lines :
        if first :
            first = False
            continue
        linesplit = line.split()
        assert( len(linesplit) == 3 )
        vals.append( float(linesplit[1]) )
        errs.append( float(linesplit[2]) )
    return vals, errs

# Fit the data to y = A exp(k*x)
def do_data_fitting( all_DG_foldings, fitting_ic50_vals, fitting_ic50_errs ) :
    x_data = np.zeros( len(all_DG_foldings), dtype=np.float64 )
    for  i in range( len(all_DG_foldings) ) :
        x_data[i] = all_DG_foldings[i][1]
    y_data = np.array( np.log(fitting_ic50_vals), dtype=np.float64 )
    y_err = np.array( np.log(fitting_ic50_errs), dtype=np.float64 )

    print( "x_data\t", x_data )
    print( "y_data\t", y_data )
    print( "y_err\t", y_err )

    k, A, rval, pval, stderr = linregress(x_data, y_data )

    return A, k, rval**2


# Analyze the folding funnels.
# Suffices and expected values:
suffices = ["1A", "1B", "1C", "1D", "1E", "1F", "1G"]
expected_min_pnears = [ 0.83, 0.68, 0.57, 0.80, 0.75, 0.88, 0.90 ]
expected_energy_gaps = [ 3, 3, 3, 3, 3, 3, 5 ]
expected_min_rmsds = [ 0.45, 0.7, 0.55, 0.25, 0.3, 0.22, 0.25 ]
if debug == True :
    for i in range( len(expected_min_rmsds) ) :
        expected_min_rmsds[i] = expected_min_rmsds[i] + 0.2 #More permissive thresholds in debug mode, since there's less sampling.
expected_max_rmsds = [ 2.5, 2.4, 2.9, 2.4, 2.2, 2.5, 2.6 ]
if debug == True :
    for i in range( len(expected_max_rmsds) ) :
        expected_max_rmsds[i] = expected_max_rmsds[i] - 0.2 #More permissive thresholds in debug mode, since there's less sampling.
expected_rmsds_of_lowest = [ 0.45, 0.85, 1.05, 0.55, 0.65, 0.32, 0.32 ]
# Will be a list of tuples of (suffix, DG_folding, DG_folding_to_lowest):
all_DG_foldings = []
for i in range( len(suffices) ) :
    all_DG_foldings.append( analyze_one_funnel( suffix=suffices[i], expected_rmsd_of_lowest=expected_rmsds_of_lowest[i], expected_min_rmsd=expected_min_rmsds[i], expected_max_rmsd=expected_max_rmsds[i], expected_min_pnear=expected_min_pnears[i], expected_energy_gap=expected_energy_gaps[i] ) )

print( all_DG_foldings )

# Do the data fitting:
fitting_ic50_vals, fitting_ic50_errs = load_ic50_data( f'{working_dir}/inputs/ic50_data.txt' )
print( "IC50_vals:\t", fitting_ic50_vals )
print( "IC50_errors:\t", fitting_ic50_errs )

fitted_A, fitted_k, fitted_Rsq = do_data_fitting( all_DG_foldings, fitting_ic50_vals, fitting_ic50_errs )
print( "Fitted A:\t", fitted_A )
print( "Fitted k:\t", fitted_k )
print( "R^2:\t", fitted_Rsq )

bundle = [fitted_A, fitted_k, fitted_Rsq]

# Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
benchmark.save_variables('''fitting_ic50_vals fitting_ic50_errs all_DG_foldings bundle
working_dir testname enough_sampling pnear_good pnear_to_lowest_good lowest_E_close_enough
sampling_under_0_25_A sampling_beyond_1_5_A sampling_beyond_2_6_A big_energy_gap overall_pass
logfile_names
''')
