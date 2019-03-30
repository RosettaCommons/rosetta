#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  cofactor_binding_sites/2.analyze.py
## @brief this script is part of cofactor_binding_sites scientific test
## @author Amanda Loshbaugh

import os
import sys
import math
#import operator
import subprocess
import glob

import benchmark
from benchmark import quality_measures as qm

background = { # background amino acid probability in all proteins. These numbers are taken from the analysis scripts used in Ollikainen, N., de Jong, R. M., and Kortemme, T. Coupling protein side-chain and backbone flexibility improves the re-design of protein-ligand specificity. PLoS Comput Biol, 2015.
'A' : 0.0853130414059,
'C' : 0.0145091808885,
'E' : 0.0697042211031,
'D' : 0.0576610517405,
'G' : 0.0677683836625,
'F' : 0.0368651011894,
'I' : 0.0658157819481,
'H' : 0.0211289643495,
'K' : 0.0581917850968,
'M' : 0.0190262038642,
'L' : 0.0958638899794,
'N' : 0.0369395202374,
'Q' : 0.036293485414,
'P' : 0.0391082335344,
'S' : 0.0594265039867,
'R' : 0.0562652852139,
'T' : 0.0541996845528,
'W' : 0.0108604669712,
'V' : 0.0866667775459,
'Y' : 0.0283924373158
}

def get_RankTop( systems_Q, systems_P ):
	# Predicted rank of the most frequent experimentally observed amino acid.
	# In case of ties in prediction, the maximum rank was used.
	#print("get_RankTop")
	design_RankTop = {}
	for system in sorted(systems):
		pdb = system.split()[0]
		system_RankTop = {}
		design_count_dict = systems_Q[ pdb ][0]
		length_Q = systems_Q[ pdb ][2]
		natural_frequency_dict = systems_P[ pdb ][0]
		
		for position in range( 0, length_Q ):
			max_natural_freq = 0
			max_natural_aa = ''
			for aa in sorted( natural_frequency_dict[position] ):
				if natural_frequency_dict[position][aa] > max_natural_freq:
					max_natural_freq = natural_frequency_dict[position][aa]
					max_natural_aa = aa
			rank = 1
			aa_design_freq = design_count_dict[position][max_natural_aa]
			for aa in sorted( design_count_dict[position] ):
				if design_count_dict[position][aa] > aa_design_freq:
					rank += 1
			if aa_design_freq == 0:
				rank = 20
			system_RankTop[position] = rank
		design_RankTop[ pdb ] = system_RankTop
	return design_RankTop

def get_freq( sequences_dict, natural_or_design ):
	#print("get_freq")
	systems_Q = {}
	for system in systems:
		pdb = system.split()[0]
		#print "system="+pdb
		Q = {}
		sequences = sequences_dict[ pdb ][0]
		seq_length = sequences_dict[ pdb ][1]
		# Count up probability of each aa occurring at a position
		aa = 'ACDEFGHIKLMNPQRSTVWY'
		for i in range( 0, seq_length ):
			Q[i] = single.copy()
		count_Q = 0.0
		# Count
		for sequence in range(0, len(sequences)):
			if len( sequences[sequence] ) != seq_length:
				print("{0} {1} sequence ({2}) is not the correct length! Expected length = {3}\n".format(natural_or_design, system, sequences[sequence], seq_length ))
				continue
			count_Q += 1.0
			for i in range( 0, seq_length ):
				result_aa = sequences[sequence][i]
				Q[i][result_aa] += 1.0
		# End count
		# Convert count to frequency
		for i in range( 0, seq_length ):
			for aa in Q[i]:
				Q[i][aa] = float(Q[i][aa]) / float(count_Q)
		systems_Q[ pdb ] = [ Q, count_Q, seq_length ]
	return systems_Q

def get_freq_null( design_sequences_dict ):
	#print("get_freq")
	systems_Q = {}
	for system in systems:
		pdb = system.split()[0]
		Q = {}
		seq_length = design_sequences_dict[ pdb ][1]
		aa = 'ACDEFGHIKLMNPQRSTVWY'
		for i in range( 0, seq_length ):
			Q[i] = single.copy()
		for i in range( 0, seq_length ):
			for aa in Q[i]:
				Q[i][aa] = float(1) / float(20)
		count_Q = 1.0
		systems_Q[ pdb ] = [ Q, count_Q, seq_length ]
	return systems_Q

def get_similarity( systems_Q, systems_P ):
	#print("get_similarity")
	#design_divergences = {}
	design_similarity = {}
	
	for system in systems:
		pdb = system.split()[0]
		#print "system="+pdb
		system_similarity = {}
		designed_frequency_dict = systems_Q[ pdb ][0]
		count_Q = systems_Q[ pdb ][1]
		length_Q = systems_Q[ pdb ][2]

		natural_frequency_dict = systems_P[ pdb ][0]
		count_P = systems_P[ pdb ][1]
		length_P = systems_P[ pdb ][2]
		
		assert (  length_Q == length_P ), "\nlength_Q != length_P\nlength_Q={0}\nlength_P={1}".format( length_Q, length_P )

		for position in range( 0, length_Q ):
			similarity = {}
			sum_PR = 0.0
			sum_QR = 0.0
			sum_bg_RR = 0.0
			sum_R_RR = 0.0
			for a in natural_frequency_dict[position]: #for each possible amino acid at that position
				# Discrete probability distribution of that aa in sequence
				P_i = natural_frequency_dict[position][a]# float(natural_frequency_dict[position][a]) / float(count_P)
				Q_i = designed_frequency_dict[position][a] #float(designed_frequency_dict[position][a]) / float(count_Q)
				R_i = P_i * 0.5 + Q_i * 0.5 # M
				bg_i = background[a]
				RR_i = bg_i * 0.5 + R_i * 0.5 # bg_i = background aa probability in all proteins; R_i = P_i * 0.5 + Q_i * 0.5
				if R_i != 0 and Q_i != 0:
					sum_PR += math.log( (Q_i / R_i), 2 ) * Q_i
				if R_i != 0 and P_i != 0:
					sum_QR += math.log( (P_i / R_i), 2) * P_i
				if RR_i != 0 and bg_i != 0:
					sum_bg_RR += math.log( (bg_i / RR_i), 2) * bg_i
				if RR_i != 0 and R_i != 0:
					sum_R_RR += math.log( (R_i / RR_i), 2) * R_i
			divergence = float(sum_PR + sum_QR) * 0.5
			significance = float(sum_bg_RR + sum_R_RR) * 0.5
			#final_similarity = ((1 - divergence) * (1 + significance)) * 0.5
			final_similarity = (1 - divergence)
			system_similarity[ position ] = final_similarity
		design_similarity[ pdb ] = system_similarity
	return design_similarity

def get_cm_fasta( fasta_files, pdb ):
	#print('pdb={0}'.format(pdb))
	seqs = []
	for fasta_file in fasta_files:
		f = open(fasta_file)
		for line in f:
			if ">" not in line:
				seqs.append(line.strip().replace('\n',''))
	to_skip = {
		'1ZK4': set([21,22,23,24,25,26,27]), #same as to_skip.txt
		'3DLC': set([0,1,2,3]), #same as to_skip.txt
		'3DK9': set([30,31,32,33,34,35]), # use for Coupled Moves analysis, reproduces Fig 6 paper result
		'3R2Q': set([11,12]), #same as to_skip.txt
	}
	for i in range( 0, len(seqs) ):
		seq = seqs[i]
		edited = ''
		for j in range(0, len(seq)):
			if pdb in to_skip:
				if j in to_skip[pdb]:
					continue
			edited += seq[j]
		seqs[i] = edited
	return seqs

def pdb_to_fasta( ):
	#if sys.argv[2] == "f":
	#	print('pdb to fasta...')
	#	os.system( f'python {rosetta_dir}/tests/scientific/data/{testname}/pdb_to_fasta.py {working_dir}/' )
	# Gather result fasta files
	system_fastas = {}
	for system in systems:
		#print('1-{0}'.format(system))
		pdb = system.split()[0]
		system_fastas[ pdb ] = []
		# dorrrk = f'{working_dir}/output/*{pdb}*fasta'
		# print(dorrrk)
		fasta_files = glob.glob( f'{working_dir}/output/*{pdb}*fasta' )
		for fasta_file in fasta_files:
			#print(fasta_file)
			system_fastas[ pdb ].append( fasta_file )
	return system_fastas

def slice_fasta( ):
	#print("Making sliced pdb for analysis")
	system_fastas = pdb_to_fasta()
	design_sequences_dict = {}
	natural_sequences_dict = {}
	for system in sorted(systems):
		pdb = system.split()[0]
		
		# vvv COMMENT THIS REGION WHEN ANALYZING COUPLED MOVES FIXBB WITH 0 TRIALS vvv
		#     OR IF SLICED FASTAS ARE ALREADY PREPARED
		fasta_files = system_fastas[ pdb ] #[:n_fasta]
		#print '{0} {1} fasta files'.format( pdb, len(fasta_files) )
		#if sys.argv[2] == "cm":
		seq_list = get_cm_fasta( fasta_files, pdb )
		#if sys.argv[2] == "f":
		#	seq_list = get_f_seq( fasta_files, pdb )
		
		# Comment to prevent overwriting existing slices:
		slice_file = f'{working_dir}/output/{pdb}.fasta'
		if os.path.isfile(slice_file) == False or os.stat(slice_file).st_size == 0:
			output = open( ( slice_file ), 'w')
			for i in range(0, len(seq_list)):
				output.write(">"+str(i)+"\n")
				output.write(seq_list[i]+"\n")
			output.close()
		# ^^^ COMMENT THIS REGION WHEN ANALYZING COUPLED MOVES FIXBB WITH 0 TRIALS ^^^
		
		# vvv DON'T COMMENT
		designed_sequences = []
		design_length = 0
		fasta_file = f'{working_dir}/output/{pdb}.fasta' #'{0}/{1}.fasta'.format( dir, pdb )
		
		f = open(fasta_file)
		for line in f:
			if ">" not in line:
				line = line.strip()
				designed_sequences.append(line)
				design_length = len(line)
				#print( 'design sequence = {0}</>'.format( line.strip() ) )
		#print('{0} {1} fasta files --> {2} design sequences'.format( pdb, len(fasta_files), len(designed_sequences) ) )
		design_sequences_dict[ pdb ] = [ designed_sequences, design_length ]
		### ^^^ DON'T COMMENT
		
		# Extract natural sequences
		# Natural sequences are pre-sliced 
		natural_sequences = []
		natural_length = 0
		f = f'{rosetta_dir}/tests/scientific/data/{testname}/natural_seqs/{pdb}_binding_site.txt'
		f = open(f)
		for line in f:
			if ">" not in line:
				if "X" not in line and "U" not in line:
					natural_sequences.append(line.strip().upper())
					natural_length = len(line.strip())
		#print(natural_sequences[0:4])
		natural_sequences_dict[ pdb ] = [ natural_sequences, natural_length ]
		# Done extracting natural sequences
		assert (  design_length == natural_length ), "\ndesign_length != natural_length\ndesign_length={0}\nnatural_length={1}".format( design_length, natural_length )
	return design_sequences_dict, natural_sequences_dict

def write( design_natural_similarity, null_natural_similarity, design_RankTop ):
	results = {}
	out = open( f"{working_dir}/output/{testname}_results.txt", "w")
	out.write("pdb aa design_natural_similarity null_natural_similarity RankTop\n")
	for system in sorted(design_natural_similarity):
		for position in sorted(design_natural_similarity[ system ]):
			key = '{0}_{1}'.format( system, position )
			results[ key ] = {}
			dn_similarity = design_natural_similarity[ system ][position]
			nn_similarity = null_natural_similarity[ system ][position]
			RankTop = design_RankTop[ system ][position]
			results[ key ][ 'design PPS' ] = dn_similarity
			results[ key ][ 'null PPS' ] = nn_similarity
			results[ key ][ 'design RT' ] = RankTop
			out.write( "{0} {1} {2} {3} {4}\n".format( system, position, dn_similarity, nn_similarity, RankTop ) )
	out.close()
	return results

benchmark.load_variables()  # Python black magic: load all variables saved by previous script
config = benchmark.config()
systems = open( f'{rosetta_dir}/tests/scientific/data/{testname}/input/systems.txt' ).readlines()
aa = 'ACDEFGHIKLMNPQRSTVWY'
single = {}
for a in aa:
	single[a] = 0.0

design_sequences_dict, natural_sequences_dict = slice_fasta( )
design_freq = systems_Q = get_freq( design_sequences_dict, 'design' )
natural_freq = systems_P = get_freq( natural_sequences_dict, 'natural' )
null_freq = get_freq_null( design_sequences_dict )
design_natural_similarity = get_similarity( design_freq, natural_freq )
null_natural_similarity = get_similarity( null_freq, natural_freq )
design_RankTop = get_RankTop( design_freq, natural_freq )

results = write( design_natural_similarity, null_natural_similarity, design_RankTop )

benchmark.save_variables( 'results targets working_dir testname')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)