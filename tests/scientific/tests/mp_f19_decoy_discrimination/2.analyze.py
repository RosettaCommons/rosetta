#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  mp_f19_decoy_discrimination/1.submit.py
## @brief this script is part of the franklin19 decoy discrimination scientific test
## @author Sergey Lyskov, Rebecca F. Alford (ralford3@jhu.edu)

import os, sys, subprocess, math
import numpy as np
import benchmark
from benchmark.util import quality_measures as qm

benchmark.load_variables()  # Python black magic: load all variables saved by previous script
config = benchmark.config()

results = {}
scorefiles = []
cutoffs_sampled_rms_dict = {}
cutoffs_weighted_rms_dict = {}
cutoffs_pnear_dict = {}
failures = []

# scorefiles and logfiles
scorefiles.extend( [ f'{working_dir}/output/{t}/{t}.rescore' for t in targets ] )

outfile = "result.txt"
cutoffs = "cutoffs"

# read cutoffs
protein = subprocess.getoutput( "grep -v '#' " + cutoffs + " | awk '{print $1}'" ).splitlines()
cutoffs_sampled_rms = subprocess.getoutput( "grep -v '#' " + cutoffs + " | awk '{print $2}'" ).splitlines()
cutoffs_pnear = subprocess.getoutput( "grep -v '#' " + cutoffs + " | awk '{print $3}'" ).splitlines()
cutoffs_weighted_rms = subprocess.getoutput( "grep -v '#' " + cutoffs + " | awk '{print $4}'" ).splitlines()
cutoffs_sampled_rms = map( float, cutoffs_sampled_rms )
cutoffs_pnear = map( float, cutoffs_pnear )
cutoffs_weighted_rms = map( float, cutoffs_weighted_rms )

cutoffs_sampled_rms_dict.update( dict( zip ( protein, cutoffs_sampled_rms )))
cutoffs_pnear_dict.update( dict( zip ( protein, cutoffs_pnear )))
cutoffs_weighted_rms_dict.update( dict( zip ( protein, cutoffs_weighted_rms )))

# open results output file
with open( outfile, "w" ) as f:
	f.write( "#pdb\tsampled_rms\tpnear\tweighted_rms\n" )

	for target in targets: 
		
		# Rescore the pdbs to compute RMS values
		get_refined = subprocess.check_output(f'ls {working_dir}/output/{target}/*_0001.pdb', shell=True).decode("utf-8").split('\n')
		rescore_list_file = f'{working_dir}/output/{target}/refined_pdbs.txt'

		with open( rescore_list_file, 'w' ) as g:
			for relaxed_pdb in get_refined[:-1]: 
				g.write( relaxed_pdb + "\n" )

		executable = f'{rosetta_dir}/source/bin/score_jd2.{extension}'
		arguments = f"-database {rosetta_dir}/database -in:file:native {rosetta_dir}/tests/scientific/data/{testname}/{target}/{target}_native.pdb -mp:setup:spanfiles {rosetta_dir}/tests/scientific/data/{testname}/{target}/{target}.span -in:membrane -mp:lipids:has_pore false -in:file:l {rescore_list_file} -score:weights franklin2019 -out:file:scorefile {working_dir}/output/{target}/{target}.rescore"
		cmd = executable + " " + arguments; print(f'{cmd!r}')
		benchmark.execute('Running Rosetta', cmd )

		# Calculate the discrimination score
		decoy_disc_data = subprocess.check_output( f'python3 score_energy_landscape.py -terms rms total_score -abinitio_scorefile {working_dir}/output/{target}/{target}.rescore', shell=True ).decode("utf-8").replace('\n', '\t').split('\t')

		if ( len(decoy_disc_data) != 8 ): 
			sys.exit( "Score energy landscape script failed!" )

		sampled_rms = round( float(decoy_disc_data[3]), 3 )
		pnear = round( float(decoy_disc_data[4]), 3 )
		weighted_rms = round( float(decoy_disc_data[5]), 3)

		target_results = {}

		# check for sampled_rms below cutoff
		f.write( target + "\t" )
		f.write( str(sampled_rms) + "\t" )
		if ( sampled_rms > cutoffs_sampled_rms_dict[ target ] ): 
			failures.append( "sampled_rms " + target + " value: " + str(sampled_rms) )

		# check for pnear better than cutoff
		# NOTE: Pnear isn't actually used here because it requires better sampling
		f.write( str(pnear) + "\t" )
		# if ( pnear < cutoffs_pnear_dict[ target ] ): 
		# 	failures.append( target )

		# check for weighted rms better than cutoff
		f.write( str(weighted_rms) )
		if ( weighted_rms > cutoffs_weighted_rms_dict[ target ] ): 
			failures.append( "weighted_rms " + target + " value: " + str(weighted_rms) )

		results.update( {target : [sampled_rms, pnear, weighted_rms]} )
		f.write( "\n" )

benchmark.save_variables('debug targets nstruct working_dir testname results scorefiles cutoffs_sampled_rms_dict cutoffs_pnear_dict cutoffs_weighted_rms_dict failures')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
