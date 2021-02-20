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
from collections import defaultdict
import numpy as np
import benchmark
import pandas, json
from typing import *

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into 	s
config = benchmark.config()


#=================================
def read_res_size(f='pdb_roots_five.txt') -> Dict[str, int]:
	"""
	Reads a delimited text file of each PDB, Glycan start point, glycan tree length, and resolution.
	Ex:

	#PDB Branches Residues resolution
	1jnd	200A	4	1.3

	Returns a dictionary with the branch name (pdb_start) and length  ex: 3uue_45A
	"""
	branch_lengths = defaultdict()
	lines = open(f, 'r').readlines()
	for line in lines:
		line = line.strip()
		if not line: continue
		if line.startswith('#'): continue
		lineSP = line.split()

		pdb = lineSP[0]
		branches = lineSP[1].split(',')
		lengths = lineSP[2].split(',')

		for x in zip(branches, lengths):
			branch_lengths[pdb + '_' + x[0]] = int(x[1])

	return branch_lengths


branch_lengths = read_res_size()
#print(branch_lengths)

#=================================
def load_json_scorefile(local_lines: List[str]) -> pandas.DataFrame:
	"""
	Read scorefile lines as a dataframe, sorted by total_score with Nan's correctly replaced.
	"""
	decoys=[]
	for line in local_lines:
			o = json.loads(line.replace("nan", "NaN"))
			# print o[self.decoy_field_name]
			# print repr(o)
			decoys.append(o)
	local_df = pandas.DataFrame.from_dict(decoys)
	local_df = local_df.infer_objects()
	# df.to_csv("debugging.csv", sep=",")

	local_df = local_df.sort_values("total_score", ascending=True)

	return local_df

#=================================
def get_pdbID(name) -> AnyStr:
	"""
	Get the PDB ID from a whole name.
	:param name:
	:return:
	"""
	return os.path.basename(name).split('_')[-3]

#=================================
def get_decoy_name(name) -> str:
	"""
	Get the basename of the decoy.
	"""
	return os.path.basename(name)

#=================================
def get_branchID(name) -> AnyStr:
	"""
	Get the branch ID such as 54A
	:param name:
	:return:
	"""
	return os.path.basename(name).split('_')[-4]

#=================================
def order_by_position(d: pandas.DataFrame):
	"""
	Orders the dataframe by pdb_branch and decoyname.
	:param d:
	:return:
	"""
	d.reset_index(inplace=True)
	d.set_index(['pdb_branch', 'decoy_name'], inplace=True)
	d.sort_index(inplace=True)
	d.reset_index(inplace=True)

#=================================
def get_branch_length(name) -> int:
	"""
	Get the length of a particular branch
	:param name:
	:return:
	"""
	return branch_lengths[name]

#=================================
#Concatonate all the scores
os.chdir("decoys")
os.system("cat score.sc* > scores.json")
os.chdir(working_dir)

lines = open("decoys/scores.json").readlines()
df = load_json_scorefile(lines)

df['pdbID'] = df['decoy'].apply(get_pdbID)
df['decoy_name'] = df['decoy'].apply(get_decoy_name)
df['pdb_branch'] = df['pdbID'] + '_' + df['branch_selection']
df['branch_length'] = df['pdb_branch'].apply(get_branch_length)
df['pdb-branch'] = df['pdb_branch'].apply(lambda x: x.replace('_', ' ').upper())
df['pdb-branch-low'] = df['pdb_branch'].apply(lambda x: x.replace('_', ' '))
df['pdb-branch-low-size'] = df['pdb-branch-low'] + ", " + df['branch_length'].astype(str)
df['experiment'] = df['opt_job_tag']
df['experiment-neat'] = df['experiment'].apply(lambda x: x.replace('-', ' ').title())


order_by_position(df)


#from benchmark import quality_measures as qm

# Basic Enrichments

full_enrichments=[]
totals = defaultdict(int)
for name, f in df.reset_index().groupby(['experiment']):
	#print('\t', name)
	#print('\t', len(f))
	branches = f['pdb_branch'].unique()
	#print('\t', 'branches', len(branches))

	# Mean
	n_rmsd = []
	for name2, f2 in f.groupby(['pdb_branch']):
		#print(name2)
		#print(len(f2))
		n_rmsd.append(len(f2))


		for rmsd in [1.0, 2.5, 5.0]:
			#print('\n')
			#print(n)
			enrichments = defaultdict()
			enrichments["pdb_branch"] = name2
			enrichments["n"] = 0
			enrichments["experiment"] = name
			enrichments["fit6_rmsd"] = rmsd
			enrichments["experiment"] = name

			filt_df = f2[f2['fit6_rmsd'] <= rmsd]

			enriched_totals = len(filt_df)
			enrichments["n"] = enriched_totals
			totals[rmsd]+= enriched_totals

			full_enrichments.append(enrichments)


		enrichments = defaultdict()


	# print('\t','mean/sd', np.mean(n_rmsd), " +/-",np.std(n_rmsd))

# Write basic enrichment results and the full df, which we will load back in.
enrichment_df = pandas.DataFrame.from_dict(full_enrichments)

df.to_csv("data.csv")
enrichment_df.to_csv("sampling_enrichments.csv", index=False, columns=['experiment', 'pdb_branch', 'fit6_rmsd', 'n'])


# Fail if sampling fails to enrich for low-rmsd models.  We don't really care about the score function here - only sampling.
#  More specific cutoffs will come after testing, and after the actual benchmarking for the paper is complete.
#  1k outputs isn't quite enough to get a good funnel for some of these glycans, but is enough to start enriching for low-RMSD models.

#  These were determined based on results from benchmarking for the paper at 1k nstruct
any_failures = []
branch_failures = []
specific_failures = []

must_pass_one_angstrum = ["1jnd_200A", "4nyq_35A"]
must_pass_five_angstrum= ["1jnd_200A", "4nyq_35A", "3pfx_267A", "1gai_171A"]
for d in full_enrichments:
	if d["n"] > 0: continue

	#print(repr(d))
	rmsd = d["fit6_rmsd"]
	pdb_branch = d["pdb_branch"]


	#All should pass 5 A enrichment, even at 1k nstruct.
	if rmsd == 5.0:
		branch_failures.append(d["pdb_branch"]+":\tALL Required - No Enrichment at 5.0 RMSD ")
	elif rmsd == 2.5 and pdb_branch in must_pass_five_angstrum:
		specific_failures.append(d["pdb_branch"] + ": \tInput Specific - No Enrichment at 2.5 RMSD ")
	elif rmsd == 1.0 and pdb_branch in must_pass_one_angstrum:
		specific_failures.append(d["pdb_branch"]+":\tInput Specific - No Enrichment at 1.0 RMSD ")

for rmsd in totals:
	if totals[rmsd] == 0:
		any_failures.append("ANY: No Enrichment at "+str(rmsd))

failures = []
failures.extend(any_failures)
failures.extend(branch_failures)
failures.extend(specific_failures)

with open("failures.txt", 'w') as failure_file:
	for failure in failures:
		failure_file.write(failure+"\n")


benchmark.save_variables('working_dir testname failures debug')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
#benchmark.save_variables('targets nstruct working_dir testname results scorefiles cutoffs_rmsd_dict cutoffs_score_dict failures')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
