#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  mp_f19_sequence_recovery/1.submit.py
## @brief this script is part of the franklin2019 sequence recovery test
## @author Sergey Lyskov, Rebecca F. Alford (ralford3@jhu.edu)

import os, sys, subprocess, math
import numpy as np
import benchmark
import seqrecov_metrics
from benchmark.util import quality_measures as qm

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into 	s
config = benchmark.config()

results = {}
failures = []
outfile = "result.txt"
cutoffs = "cutoffs"

# Make a list of the designed PDB files
design_list = subprocess.check_output(f'ls {working_dir}/output/*/*_0001.pdb', shell=True).decode("utf-8").split('\n')
targets_completed = [ x.split("/")[-1].split("_")[0] for x in design_list[:-1] ]
data_dir = f'{rosetta_dir}/tests/scientific/data/{testname}/'
natives_list = []
for t in range(0, len(targets_completed)): 
	curr_target = targets_completed[t]
	natives_list.append( f'{data_dir}{curr_target}/{curr_target}_tr_ignorechain.pdb' )

# Write the output lists
design_list_file = f'{working_dir}/output/design_pdbs.txt'
with open( design_list_file, 'w' ) as f:
	for designed_pdb in design_list[:-1]: 
		f.write( designed_pdb + "\n" )

native_list_file = f'{working_dir}/output/native_pdbs.txt'
with open( native_list_file, 'w' ) as g: 
	for native_pdbs in natives_list: 
		g.write( native_pdbs + "\n" )

# Calculate the raw sequence recovery statistics
sequence_recovery_output = "mp_f19_design_info.txt"
executable = f'{rosetta_dir}/source/bin/mp_seqrecov.{extension}'
arguments = f'-database {rosetta_dir}/database -native_pdb_list {native_list_file} -redesign_pdb_list {design_list_file} -seq_recov_filename {working_dir}/output/{sequence_recovery_output} -read_only_ATOM_entries -in:ignore_unrecognized_res'
cmd = executable + " " + arguments; print(f'{cmd!r}')
benchmark.execute('Running Rosetta', cmd )

# Using the sequence recovery statistics, calculate Naa, Dkl and Rseq
seqdata = seqrecov_metrics.read_sequence_data( f'output/{sequence_recovery_output}' )
non_random_recov = seqrecov_metrics.compute_non_random_recovered_fraction( seqdata )
recov = seqrecov_metrics.compute_sequence_recovery( seqdata )
kl_divergence = seqrecov_metrics.compute_kl_divergence( seqdata )

# read cutoffs
cutoff_info = subprocess.getoutput( "grep -v '#' " + cutoffs + " | awk '{print $0}'" ).splitlines()[0].split(' ')
cutoffs_recov = float( cutoff_info[0] )
cutoffs_non_random = float( cutoff_info[1] )
cutoffs_kl_div_subset = float( cutoff_info[2] )
cutoffs_kl_div_all = float( cutoff_info[3] )

# open results output file
f = open( outfile, "w" )

subset = "all lipid aqueous".split()

# check for non_random above cutoff
for v in range( 0, len(non_random_recov)): 
	f.write( "non_random\t" + subset[v] + "\t" + str(non_random_recov[v]) + "\n")
	if non_random_recov[v] < cutoffs_non_random: 
		failures.append( "non_random\n" ) 

# check for recov above cutoff
for v in range( 0, len(recov)): 
	f.write( "recovery\t" + subset[v] + "\t" + str(round(recov[v],3)) + "\n")
	if recov[v] < cutoffs_recov: 
		failures.append( "recovery\n" )

# check for kl values below cutoff
for v in range( 0, len(recov)): 
	f.write( "kl_divergence\t" + subset[v] + "\t" + str(round(kl_divergence[v],3)) + "\n")
if kl_divergence[0] > cutoffs_kl_div_all: 
	failures.append( "kl_div_all\n" )

# check for kl values below cutoff
if kl_divergence[1] > cutoffs_kl_div_subset and kl_divergence[2] > cutoffs_kl_div_subset: 
	failures.append( "kl_div_subset\n" )

f.close()

benchmark.save_variables('debug targets nstruct working_dir testname results scorefiles cutoffs_non_random cutoffs_recov cutoffs_kl_div_subset cutoffs_kl_div_all non_random_recov recov kl_divergence failures')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
