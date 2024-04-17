#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  simple_cycpep_predict/2.analyze.py
## @brief This script is part of simple_cycpep_predict scientific test.
## @author Sergey Lyskov
## @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).

import os, sys, subprocess, math
import re
import benchmark

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into 	s
config = benchmark.config()

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

samples_expected = 6000 if debug else 230000

# If bool is true, print YES; otherwise print NO
def bool_to_string( bool_in ):
	assert( bool_in == True or bool_in == False )
	if( bool_in == True ):
		return "YES"
	return "NO"

# Relevant filenames
outfile = f'{working_dir}/result.txt'


for f in f'.hpc.{testname}.output.0.log .hpc.{testname}.output.log .hpc.{testname}.output'.split():
    logfile = f'{working_dir}/hpc-logs/{f}'
    if os.path.isfile(logfile): break


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
    except:
        # Provide a slightly more useful debugging info for outputs
        print("ISSUE IN", logfile, "can't parse", line)
        raise

print( "Read PNear=" + str(pnear) + " from " + logfile + "." )
print( "Read PNearLowest=" + str(pnear_to_lowest) + " from " + logfile + "." )

assert None not in [ pnear, pnear_to_lowest ]

sample_match = re.search("simple_cycpep_predict application completed [1-90]+", str( subprocess.getoutput( 'grep "application completed" ' + logfile ) ) )
if sample_match is not None:
    total_samples = int( sample_match.group().split()[-1] )
else:
    print("ERROR: Could not find the number of structures!")
    total_samples = 0

print ( "Determined that " + str( total_samples ) + " samples were performed." )

# Check that total samples is high enough.
if( total_samples > samples_expected ) :
	enough_sampling = True

# Check that PNear is reasonable:
if ( pnear >0.94 ):
	pnear_good = True

# Check that PNear to loweset energy is reasonable:
if ( pnear_to_lowest >0.97 ):
	pnear_to_lowest_good = True

# Cheack that lowest energy is first entry
lowestE = min( energy_vals )
if ( abs( energy_vals[0] - lowestE ) < 1e-8 ) :
	lowest_E_is_first = True

# Check that lowest energy is below 0.3 A.
if( rmsd_vals[0] < 0.3 ) :
	lowest_E_close_enough = True

# Check that we're sampling diverse structures
max_rmsd = max( rmsd_vals )
print( "Max RMSD = " + str( max_rmsd ) )
if( max_rmsd > 1.5 ) :
	sampling_beyond_1_5_A = True
if( max_rmsd > 2.6 ) :
	sampling_beyond_2_6_A = True

# Check that we're sampling under 0.25 A
min_rmsd = min( rmsd_vals )
print( "Min RMSD = " + str( min_rmsd ) )
if( min_rmsd < 0.25 ) :
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
	if( energy_gap_val > 6.0 ):
		big_energy_gap = True

# Determine if we passed overall
if( pnear_good == True and pnear_to_lowest_good == True and enough_sampling and lowest_E_is_first == True and lowest_E_close_enough == True and sampling_under_0_25_A == True and  sampling_beyond_1_5_A == True and sampling_beyond_2_6_A == True and big_energy_gap == True ):
	overall_pass = True

# Write out results
with open( outfile, "w" ) as f:
	f.write( "Total samples =\t" + str(total_samples) + "\n" )
	f.write( "Computed PNear =\t" + str(pnear) + "\n" )
	f.write( "Computed PNear to lowest E =\t" + str(pnear_to_lowest) + "\n" )
	f.write( "Lowest energy =\t" + str(lowestE) + " kcal/mol\n" )
	f.write( "RMSD of lowest energy =\t" + str(rmsd_vals[0]) + " Angstroms\n" )
	f.write( "Lowest RMSD =\t" + str(min_rmsd) + " Angstroms\n" )
	f.write( "Highest RMSD =\t" + str(max_rmsd) + " Angstroms\n" )
	f.write( "Energy gap (minE>1.5A - minE) =\t" + str(energy_gap_val) + " kcal/mol\n" )
	f.write( "\n" )

	if( debug == True ):
		f.write( "More than 6,000 samples?\t" + bool_to_string( enough_sampling ) + "\n" )
	else:
		f.write( "More than 230,000 samples?\t" + bool_to_string( enough_sampling ) + "\n" )
	f.write( "PNear value over 0.94?\t" + bool_to_string( pnear_good ) + "\n" )
	f.write( "PNear value to lowest E over 0.97?\t" + bool_to_string( pnear_to_lowest_good ) + "\n" )
	f.write( "Lowest energy under 0.3 A RMSD?\t" + bool_to_string( lowest_E_close_enough ) + "\n" )
	f.write( "Sampling below 0.25 A RMSD?\t" + bool_to_string( sampling_under_0_25_A ) + "\n" )
	f.write( "Sampling beyond 1.5 A RMSD?\t" + bool_to_string( sampling_beyond_1_5_A ) + "\n" )
	f.write( "Sampling beyond 2.6 A RMSD?\t" + bool_to_string( sampling_beyond_2_6_A ) + "\n" )
	f.write( "6+ kcal/mol energy gap?\t" + bool_to_string( big_energy_gap ) + "\n" )
	f.write( "OVERALL PASS?\t" + bool_to_string( overall_pass ) + "\n" )

benchmark.save_variables('working_dir testname enough_sampling pnear_good pnear_to_lowest_good lowest_E_close_enough sampling_under_0_25_A sampling_beyond_1_5_A sampling_beyond_2_6_A big_energy_gap overall_pass logfile')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
#benchmark.save_variables('targets nstruct working_dir testname results scorefiles cutoffs_rmsd_dict cutoffs_score_dict failures')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
