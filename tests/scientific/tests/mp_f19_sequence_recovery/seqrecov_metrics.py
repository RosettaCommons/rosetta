#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  mp_f19_sequence_recovery/seqrecov_metrics.py
## @brief this script is part of the franklin2019 sequence recovery test
## @author Rebecca F. Alford (ralford3@jhu.edu)

import sys, os, csv
import numpy as np 
from collections import defaultdict

def read_sequence_data( seqrecov_file ): 

	with open( seqrecov_file, 'rt' ) as f: 
		seqdata = f.readlines()
		seqdata = [ x.strip() for x in seqdata ]
		seqdata = [ x.split('\t') for x in seqdata ]

	temp_keys = []
	for key in seqdata[0]: 
		temp_keys.append( key )

	temp_data = [[] for i in range(len(seqdata[0]))]
	for k in range(0, len(seqdata[0])): 
		for d in range(1, len(seqdata)-1):
			if ( k == 0 ): 
				temp_data[k].append( seqdata[d][k] )
			else: 
				temp_data[k].append( float(seqdata[d][k]) )

	seqdict = defaultdict(list)
	for k in range(0, len(seqdata[0])): 
		seqdict[ temp_keys[k] ] = temp_data[k]
	return seqdict

def skip_when_native_is_zero( array1, array2 ): 
	return array1[np.nonzero(array2)], array2[np.nonzero(array2)]

def compute_non_random_recovered_fraction( seqdict ): 
	"""
	For each category, calculate the number of residues with an overall recovery
	rate above 0.05 (or 5%, the background probability).
	"""

	# all positions
	n_correct_all_arr = np.asarray( seqdict[ 'No.correct' ] )
	n_native_all_arr = np.asarray( seqdict[ 'No.native' ] )
	n_correct_all, n_native_all = skip_when_native_is_zero( n_correct_all_arr, n_native_all_arr )
	all_recovered = np.divide( n_correct_all, n_native_all )
	all_non_random = len(np.where( all_recovered >= 0.05 )[0])

	# lipid-facing positons
	n_correct_lipid_arr = np.asarray( seqdict[ 'No.correct.lipid.facing' ] )
	n_native_lipid_arr = np.asarray( seqdict[ 'No.native.lipid.facing' ] )
	n_correct_lipid, n_native_lipid = skip_when_native_is_zero( n_correct_lipid_arr, n_native_lipid_arr )
	lipid_recovered = np.divide( n_correct_lipid, n_native_lipid )
	lipid_non_random = len(np.where( lipid_recovered >= 0.05 )[0])

	# aqueous-facing positions
	n_correct_aqueous_arr = np.asarray( seqdict[ 'No.correct.aqueous' ] )
	n_native_aqueous_arr = np.asarray( seqdict[ 'No.native.aqueous' ] )
	n_correct_aqueous, n_native_aqueous = skip_when_native_is_zero( n_correct_aqueous_arr, n_native_aqueous_arr )
	aqueous_recovered = np.divide( n_correct_aqueous, n_native_aqueous )
	aqueous_non_random = len(np.where( aqueous_recovered >= 0.05 )[0])

	return [ all_non_random, lipid_non_random, aqueous_non_random ]

def compute_sequence_recovery( seqdict ): 
	"""
	Calculate the rate that a position is designed correctly according to the native sequence
	"""

	# all positions
	n_correct_all_arr = np.asarray( seqdict[ 'No.correct' ] )
	n_native_all_arr = np.asarray( seqdict[ 'No.native' ] )
	n_correct_all, n_native_all = skip_when_native_is_zero( n_correct_all_arr, n_native_all_arr )
	all_recovered = np.divide( np.sum(n_correct_all), np.sum(n_native_all) )
	
	# lipid-facing positons
	n_correct_lipid_arr = np.asarray( seqdict[ 'No.correct.lipid.facing' ] )
	n_native_lipid_arr = np.asarray( seqdict[ 'No.native.lipid.facing' ] )
	n_correct_lipid, n_native_lipid = skip_when_native_is_zero( n_correct_lipid_arr, n_native_lipid_arr )
	lipid_recovered = np.divide( np.sum(n_correct_lipid), np.sum(n_native_lipid) )
	
	# aqueous-facing positions
	n_correct_aqueous_arr = np.asarray( seqdict[ 'No.correct.aqueous' ] )
	n_native_aqueous_arr = np.asarray( seqdict[ 'No.native.aqueous' ] )
	n_correct_aqueous, n_native_aqueous = skip_when_native_is_zero( n_correct_aqueous_arr, n_native_aqueous_arr )
	aqueous_recovered = np.divide( np.sum(n_correct_aqueous), np.sum(n_native_aqueous) )
	
	return [all_recovered, lipid_recovered, aqueous_recovered]

def compute_kl_divergence( seqdict ): 
	"""
	Compute the sequence divergence between the designed and natives
	"""

	# all positions
	n_designed_all_arr = np.asarray( seqdict[ 'No.designed' ] )
	n_native_all_arr = np.asarray( seqdict[ 'No.native' ] )	
	n_designed_all, n_native_all = skip_when_native_is_zero( n_designed_all_arr, n_native_all_arr )
	#eliminates those with n_designed_all;
	n_designed_all[ np.where(n_designed_all == 0)[0] ] = n_native_all[ np.where(n_designed_all == 0)[0] ] 
	#n_designed_all[ n_designed_all == 0 ] = 0.0001
	all_KL = -np.sum( np.log( np.divide( n_designed_all, n_native_all ) ) )

	# lipid-facing positons
	n_designed_lipid_arr = np.asarray( seqdict[ 'No.designed.lipid.facing' ] )
	n_native_lipid_arr = np.asarray( seqdict[ 'No.native.lipid.facing' ] )
	n_designed_lipid, n_native_lipid = skip_when_native_is_zero( n_designed_lipid_arr, n_native_lipid_arr )
	n_designed_lipid[ np.where(n_designed_lipid == 0)[0] ] = n_native_lipid[ np.where(n_designed_lipid == 0)[0] ]
	#n_designed_lipid[ n_designed_lipid == 0 ] = 0.0001
	lipid_KL = -np.sum( np.log( np.divide( n_designed_lipid, n_native_lipid ) ) )

	# aqueous-facing positions
	n_designed_aqueous_arr = np.asarray( seqdict[ 'No.designed.aqueous' ] )
	n_native_aqueous_arr = np.asarray( seqdict[ 'No.native.aqueous' ] )
	n_designed_aqueous, n_native_aqueous = skip_when_native_is_zero( n_designed_aqueous_arr, n_native_aqueous_arr )
	n_designed_aqueous[ np.where(n_designed_aqueous == 0)[0] ] = n_native_aqueous[ np.where(n_designed_aqueous == 0)[0] ]
	#n_designed_aqueous[ n_designed_aqueous == 0 ] = 0.0001
	aqueous_KL = -np.sum( np.log( np.divide( n_designed_aqueous, n_native_aqueous ) ) )

	return [all_KL, lipid_KL, aqueous_KL]
